##########################################
# PFAS vs Clysis (cellular immunity)
# Which PFAS are most associated with cellular immunity?
##########################################

rm(list=ls())

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)


input_file <- "./input/2023_2024_2025_AllData_clean.csv"
output_dir <- "./output/4_PFAS_vs_Clysis"
best_shape_file <- "./output/2_clysis_pfas_relationship/clysis_best_shape_summary.csv"
key_file <- "./input/PFAS_subtypes.csv"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(log_clysis = log(clysis + 1e-6))

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

# Load best-fit shapes from Script 2 (log_clysis models)
best_shape_lookup <- setNames(character(0), character(0))
if (file.exists(best_shape_file)) {
  best_shape_tbl <- read.csv(best_shape_file, stringsAsFactors = FALSE) %>%
    filter(Response == "log_clysis") %>%
    distinct(PFAS, .keep_all = TRUE)
  best_shape_lookup <- setNames(best_shape_tbl$Best_Model, best_shape_tbl$PFAS)
}

# PFAS columns
pfas_cols <- unique(c(
  names(df %>% select(starts_with("PF"))),
  intersect(c("TOTAL_PFAS"), names(df)),
  sum_cols
))
pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal", "DNA_sex")]

covariates <- c("DNA_sex", "age")
nest_var <- "nest_no"

results_list <- list()

for(pfas in pfas_cols){
  
  if(!is.numeric(df[[pfas]]) || sum(!is.na(df[[pfas]])) < 6) next
  
  best_model_name <- ifelse(pfas %in% names(best_shape_lookup), best_shape_lookup[[pfas]], "linear")
  
  df_model <- df %>%
    filter(!is.na(.data[[pfas]]) & !is.na(log_clysis)) %>%
    mutate(
      PFAS_raw = .data[[pfas]],
      PFAS_scaled = as.numeric(scale(PFAS_raw)),
      log_PFAS_scaled = as.numeric(scale(log(PFAS_raw + 1e-6))),
      PFAS_sq_scaled = as.numeric(scale(PFAS_raw^2)),
      exp_PFAS_scaled = as.numeric(scale(exp(PFAS_scaled)))
    )
  
  # Map best-fit shape to mixed-model terms. For "power", use log-transformed PFAS
  # as a practical approximation in this mixed-model framework.
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
  
  df_model <- df_model %>% drop_na(all_of(c(covariates, focal_terms, nest_var)))
  if (nrow(df_model) < 8) next
  
  full_formula_text <- paste0(
    "log_clysis ~ ",
    paste(c(focal_terms, covariates), collapse = " + "),
    " + (1|", nest_var, ")"
  )
  null_formula_text <- paste0(
    "log_clysis ~ ",
    paste(covariates, collapse = " + "),
    " + (1|", nest_var, ")"
  )
  
  model_full <- try(lmer(as.formula(full_formula_text), data=df_model, REML=FALSE), silent=TRUE)
  model_null <- try(lmer(as.formula(null_formula_text), data=df_model, REML=FALSE), silent=TRUE)
  
  if(inherits(model_full, "try-error") || inherits(model_null, "try-error")) next
  
  lr_tbl <- anova(model_null, model_full)
  p_model <- lr_tbl$`Pr(>Chisq)`[2]
  
  full_tidy <- broom.mixed::tidy(model_full)
  focal_tidy <- full_tidy %>% filter(term %in% focal_terms)
  
  primary_row <- if (nrow(focal_tidy) > 0) focal_tidy[1, ] else tibble(
    effect = "fixed",
    group = NA_character_,
    term = paste(focal_terms, collapse = " + "),
    estimate = NA_real_,
    std.error = NA_real_,
    statistic = NA_real_,
    df = NA_real_,
    p.value = NA_real_
  )
  
  tidy_res <- primary_row %>%
    mutate(
      term = paste(focal_terms, collapse = " + "),
      p.value = p_model,
      p_term = ifelse(nrow(focal_tidy) > 0, focal_tidy$p.value[1], NA_real_),
      best_shape_model = best_model_name,
      pfas_var = pfas,
      immune_var = "Clysis",
      n_obs = nrow(df_model)
    )
  
  results_list[[pfas]] <- tidy_res
}

final_results <- bind_rows(results_list) %>% mutate(FDR_p = p.adjust(p.value, method="fdr"))

write.csv(final_results, file=file.path(output_dir,"PFAS_Clysis_mixed_model_results.csv"), row.names=FALSE)

cat("Done: PFAS vs Clysis\n")
