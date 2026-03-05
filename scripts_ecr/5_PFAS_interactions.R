##########################################
# Interaction: PFAS × sex/age
##########################################

rm(list=ls())
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)

input_file <- "./input/2023_2024_2025_AllData_clean.csv"
output_dir <- "./output/5_PFAS_interactions"
best_shape_igy_file <- "./output/1_igy_pfas_relationship/best_shape_summary.csv"
best_shape_clysis_file <- "./output/2_clysis_pfas_relationship/clysis_best_shape_summary.csv"
key_file <- "./input/PFAS_subtypes.csv"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    log_igy = log(igy + 1e-6),
    log_clysis = log(clysis + 1e-6)
  )

# Build PFAS category sums (chain_, func_, age_) by bird_id
key <- read.csv(key_file, stringsAsFactors = FALSE)
pfas_present <- intersect(key$pfas_analyte, colnames(df))
sum_cols <- character(0)

if ("bird_id" %in% names(df) && length(pfas_present) > 0) {
  long_pf <- df %>%
    select(bird_id, all_of(pfas_present)) %>%
    pivot_longer(cols = all_of(pfas_present),
                 names_to = "pfas_analyte",
                 values_to = "conc") %>%
    filter(!is.na(conc)) %>%
    left_join(key, by = "pfas_analyte")
  
  make_sum_table <- function(df_long, group_var, prefix) {
    df_long %>%
      filter(!is.na(.data[[group_var]])) %>%
      group_by(bird_id, .data[[group_var]]) %>%
      summarize(sum_val = sum(conc), .groups = "drop") %>%
      pivot_wider(names_from = all_of(group_var),
                  values_from = sum_val,
                  values_fill = 0,
                  names_prefix = prefix)
  }
  
  chain_sums <- make_sum_table(long_pf, "chain_cat", "chain_")
  func_sums  <- make_sum_table(long_pf, "func_cat",  "func_")
  age_sums   <- make_sum_table(long_pf, "pfas_age",  "age_")
  
  df <- df %>%
    left_join(chain_sums, by = "bird_id") %>%
    left_join(func_sums,  by = "bird_id") %>%
    left_join(age_sums,   by = "bird_id")
  
  sum_cols <- grep("^(chain_|func_|age_)", colnames(df), value = TRUE)
  df[sum_cols][is.na(df[sum_cols])] <- 0
}

# Ensure TOTAL_PFAS exists for modeling
pfas_analytes <- names(df %>% select(starts_with("PF")))
if (!"TOTAL_PFAS" %in% names(df) && length(pfas_analytes) > 0) {
  pfas_matrix <- df[, pfas_analytes, drop = FALSE]
  df$TOTAL_PFAS <- rowSums(pfas_matrix, na.rm = TRUE)
  df$TOTAL_PFAS[rowSums(!is.na(pfas_matrix)) == 0] <- NA_real_
}

# Load best-shape lookups per immune endpoint
best_shape_igy_lookup <- setNames(character(0), character(0))
if (file.exists(best_shape_igy_file)) {
  best_shape_igy_tbl <- read.csv(best_shape_igy_file, stringsAsFactors = FALSE) %>%
    filter(Response == "log_IgY") %>%
    distinct(PFAS, .keep_all = TRUE)
  best_shape_igy_lookup <- setNames(best_shape_igy_tbl$Best_Model, best_shape_igy_tbl$PFAS)
}

best_shape_clysis_lookup <- setNames(character(0), character(0))
if (file.exists(best_shape_clysis_file)) {
  best_shape_clysis_tbl <- read.csv(best_shape_clysis_file, stringsAsFactors = FALSE) %>%
    filter(Response == "log_clysis") %>%
    distinct(PFAS, .keep_all = TRUE)
  best_shape_clysis_lookup <- setNames(best_shape_clysis_tbl$Best_Model, best_shape_clysis_tbl$PFAS)
}

pfas_cols <- unique(c(
  names(df %>% select(starts_with("PF"))),
  intersect(c("TOTAL_PFAS"), names(df)),
  sum_cols
))

pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal", "DNA_sex")]

nest_var <- "nest_no"

immune_vars <- c("log_igy","log_clysis")

results_list <- list()

for(immune in immune_vars){
  for(pfas in pfas_cols){
    
    if(!is.numeric(df[[pfas]]) || sum(!is.na(df[[pfas]])) < 8) next
    
    df_model <- df %>%
      filter(!is.na(.data[[pfas]]),
             !is.na(.data[[immune]]),
             !is.na(DNA_sex),
             !is.na(age),
             !is.na(.data[[nest_var]])) %>%
      mutate(
        PFAS_raw = .data[[pfas]],
        PFAS_scaled = as.numeric(scale(PFAS_raw)),
        log_PFAS_scaled = as.numeric(scale(log(PFAS_raw + 1e-6))),
        PFAS_sq_scaled = as.numeric(scale(PFAS_raw^2)),
        exp_PFAS_scaled = as.numeric(scale(exp(PFAS_scaled)))
      )
    
    best_model_name <- if (immune == "log_igy") {
      ifelse(pfas %in% names(best_shape_igy_lookup), best_shape_igy_lookup[[pfas]], "linear")
    } else {
      ifelse(pfas %in% names(best_shape_clysis_lookup), best_shape_clysis_lookup[[pfas]], "linear")
    }
    
    # For "power", use log-transformed PFAS as a practical approximation.
    if (best_model_name == "linear") {
      focal_terms <- c("PFAS_scaled")
    } else if (best_model_name == "log") {
      focal_terms <- c("log_PFAS_scaled")
    } else if (best_model_name == "quadratic") {
      focal_terms <- c("PFAS_scaled", "PFAS_sq_scaled")
    } else if (best_model_name == "exp") {
      focal_terms <- c("exp_PFAS_scaled")
    } else if (best_model_name == "power") {
      focal_terms <- c("log_PFAS_scaled")
    } else {
      focal_terms <- c("PFAS_scaled")
    }
    
    df_model <- df_model %>% drop_na(all_of(focal_terms))
    
    formula_text <- paste0(
      immune,
      " ~ (", paste(focal_terms, collapse = " + "), ") * (DNA_sex + age) + (1|", nest_var, ")"
    )
    
    model <- try(lmer(as.formula(formula_text),
                      data=df_model,
                      REML=FALSE),
                 silent=TRUE)
    
    if(inherits(model,"try-error")) next
    
    tidy_res <- broom.mixed::tidy(model, effects="fixed") %>%
      mutate(
        best_shape_model = best_model_name,
        focal_terms = paste(focal_terms, collapse = " + "),
        pfas_var = pfas,
        immune_var = immune,
        n_obs = nrow(df_model)
      )
    
    results_list[[paste(immune,pfas)]] <- tidy_res
  }
}

final_results <- bind_rows(results_list)

# Separate corrections
final_results <- final_results %>%
  mutate(
    test_type = case_when(
      grepl(":", term) & grepl("PFAS_scaled|log_PFAS_scaled|PFAS_sq_scaled|exp_PFAS_scaled", term) ~ "interaction",
      term %in% c("PFAS_scaled", "log_PFAS_scaled", "PFAS_sq_scaled", "exp_PFAS_scaled") ~ "main_effect",
      TRUE ~ "covariate"
    )
  ) %>%
  group_by(test_type) %>%
  mutate(FDR_p = p.adjust(p.value, method="BH")) %>%
  ungroup()

write.csv(final_results,
          file=file.path(output_dir,"PFAS_interactions.csv"),
          row.names=FALSE)

print("Interaction modeling complete.")
