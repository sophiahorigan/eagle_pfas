##########################################
# PFAS vs White Blood Cells (mixed models)
# - Total_count: negative binomial mixed model
# - Cell proportions: beta mixed model
##########################################

rm(list = ls())

library(tidyverse)
library(glmmTMB)
library(broom.mixed)

# --- Settings ---
input_file <- "./input/2023_2024_2025_AllData_clean.csv"
output_dir <- "./output/9_WBC_PFAS_Mixed_Models"
key_file <- "./input/PFAS_subtypes.csv"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load data ---
df <- read.csv(input_file, stringsAsFactors = FALSE)

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

# Ensure TOTAL_PFAS exists
pfas_analytes <- names(df %>% select(starts_with("PF")))
if (!"TOTAL_PFAS" %in% names(df) && length(pfas_analytes) > 0) {
  pfas_matrix <- df[, pfas_analytes, drop = FALSE]
  df$TOTAL_PFAS <- rowSums(pfas_matrix, na.rm = TRUE)
  df$TOTAL_PFAS[rowSums(!is.na(pfas_matrix)) == 0] <- NA_real_
}

# --- PFAS variables ---
pfas_cols <- unique(c(
  names(df %>% select(starts_with("PF"))),
  intersect(c("TOTAL_PFAS"), names(df)),
  sum_cols
))
pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal", "DNA_sex")]

# --- WBC outcomes ---
wbc_props <- c("Lperc", "Hperc", "Eperc", "Mperc", "Bperc")
wbc_props <- wbc_props[wbc_props %in% names(df)]
wbc_total <- "Total_count"
has_total <- wbc_total %in% names(df)

covariates <- c("DNA_sex", "age")
nest_var <- "nest_no"

# Candidate PFAS terms for shape selection
shape_term_map <- list(
  linear = c("PFAS_scaled"),
  log = c("log_PFAS_scaled"),
  quadratic = c("PFAS_scaled", "PFAS_sq_scaled")
)

results_list <- list()

for (pfas in pfas_cols) {
  if (!is.numeric(df[[pfas]]) || sum(!is.na(df[[pfas]])) < 8) next
  
  # ---------- Proportion models (beta mixed) ----------
  for (wbc in wbc_props) {
    df_model <- df %>%
      select(all_of(c(pfas, wbc, covariates, nest_var))) %>%
      filter(
        !is.na(.data[[pfas]]),
        !is.na(.data[[wbc]]),
        !is.na(.data[[nest_var]]),
        !is.na(DNA_sex),
        !is.na(age)
      ) %>%
      mutate(
        PFAS_raw = .data[[pfas]],
        PFAS_scaled = as.numeric(scale(PFAS_raw)),
        log_PFAS_scaled = as.numeric(scale(log(PFAS_raw + 1e-6))),
        PFAS_sq_scaled = as.numeric(scale(PFAS_raw^2))
      )
    
    if (nrow(df_model) < 20) next
    
    # Squeeze to open interval (0,1) for beta family
    n_local <- nrow(df_model)
    df_model$y <- ((df_model[[wbc]] * (n_local - 1)) + 0.5) / n_local
    df_model <- df_model %>% filter(is.finite(y), y > 0, y < 1)
    if (nrow(df_model) < 20) next
    
    fitted_shapes <- list()
    
    for (shape_name in names(shape_term_map)) {
      focal_terms <- shape_term_map[[shape_name]]
      if (any(!focal_terms %in% names(df_model))) next
      
      full_formula <- as.formula(
        paste0("y ~ ",
               paste(c(focal_terms, covariates), collapse = " + "),
               " + (1|", nest_var, ")")
      )
      null_formula <- as.formula(
        paste0("y ~ ",
               paste(covariates, collapse = " + "),
               " + (1|", nest_var, ")")
      )
      
      full_fit <- try(
        glmmTMB(full_formula, data = df_model, family = beta_family(link = "logit")),
        silent = TRUE
      )
      null_fit <- try(
        glmmTMB(null_formula, data = df_model, family = beta_family(link = "logit")),
        silent = TRUE
      )
      
      if (inherits(full_fit, "try-error") || inherits(null_fit, "try-error")) next
      
      fitted_shapes[[shape_name]] <- list(
        full = full_fit,
        null = null_fit,
        terms = focal_terms,
        aic = AIC(full_fit)
      )
    }
    
    if (length(fitted_shapes) == 0) next
    
    best_shape <- names(which.min(sapply(fitted_shapes, function(z) z$aic)))
    best_obj <- fitted_shapes[[best_shape]]
    
    lr_tbl <- anova(best_obj$null, best_obj$full)
    p_model <- lr_tbl$`Pr(>Chisq)`[2]
    
    tidy_full <- broom.mixed::tidy(best_obj$full, effects = "fixed")
    focal_tidy <- tidy_full %>% filter(term %in% best_obj$terms)
    primary <- if (nrow(focal_tidy) > 0) focal_tidy[1, ] else tibble(
      effect = "fixed",
      component = "cond",
      group = NA_character_,
      term = paste(best_obj$terms, collapse = " + "),
      estimate = NA_real_,
      std.error = NA_real_,
      statistic = NA_real_,
      p.value = NA_real_
    )
    
    results_list[[paste("prop", wbc, pfas, sep = "_")]] <- primary %>%
      mutate(
        term = paste(best_obj$terms, collapse = " + "),
        p.value = p_model,
        p_term = ifelse(nrow(focal_tidy) > 0, focal_tidy$p.value[1], NA_real_),
        best_shape_model = best_shape,
        PFAS = pfas,
        WBC_var = wbc,
        outcome_type = "proportion_beta",
        n = nrow(df_model)
      )
  }
  
  # ---------- Total_count models (negative binomial mixed) ----------
  if (has_total) {
    df_model <- df %>%
      select(all_of(c(pfas, wbc_total, covariates, nest_var))) %>%
      filter(
        !is.na(.data[[pfas]]),
        !is.na(.data[[wbc_total]]),
        !is.na(.data[[nest_var]]),
        !is.na(DNA_sex),
        !is.na(age),
        .data[[wbc_total]] >= 0
      ) %>%
      mutate(
        PFAS_raw = .data[[pfas]],
        PFAS_scaled = as.numeric(scale(PFAS_raw)),
        log_PFAS_scaled = as.numeric(scale(log(PFAS_raw + 1e-6))),
        PFAS_sq_scaled = as.numeric(scale(PFAS_raw^2)),
        y = .data[[wbc_total]]
      )
    
    if (nrow(df_model) >= 20) {
      fitted_shapes <- list()
      
      for (shape_name in names(shape_term_map)) {
        focal_terms <- shape_term_map[[shape_name]]
        if (any(!focal_terms %in% names(df_model))) next
        
        full_formula <- as.formula(
          paste0("y ~ ",
                 paste(c(focal_terms, covariates), collapse = " + "),
                 " + (1|", nest_var, ")")
        )
        null_formula <- as.formula(
          paste0("y ~ ",
                 paste(covariates, collapse = " + "),
                 " + (1|", nest_var, ")")
        )
        
        full_fit <- try(
          glmmTMB(full_formula, data = df_model, family = nbinom2(link = "log")),
          silent = TRUE
        )
        null_fit <- try(
          glmmTMB(null_formula, data = df_model, family = nbinom2(link = "log")),
          silent = TRUE
        )
        
        if (inherits(full_fit, "try-error") || inherits(null_fit, "try-error")) next
        
        fitted_shapes[[shape_name]] <- list(
          full = full_fit,
          null = null_fit,
          terms = focal_terms,
          aic = AIC(full_fit)
        )
      }
      
      if (length(fitted_shapes) > 0) {
        best_shape <- names(which.min(sapply(fitted_shapes, function(z) z$aic)))
        best_obj <- fitted_shapes[[best_shape]]
        
        lr_tbl <- anova(best_obj$null, best_obj$full)
        p_model <- lr_tbl$`Pr(>Chisq)`[2]
        
        tidy_full <- broom.mixed::tidy(best_obj$full, effects = "fixed")
        focal_tidy <- tidy_full %>% filter(term %in% best_obj$terms)
        primary <- if (nrow(focal_tidy) > 0) focal_tidy[1, ] else tibble(
          effect = "fixed",
          component = "cond",
          group = NA_character_,
          term = paste(best_obj$terms, collapse = " + "),
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_
        )
        
        results_list[[paste("count", wbc_total, pfas, sep = "_")]] <- primary %>%
          mutate(
            term = paste(best_obj$terms, collapse = " + "),
            p.value = p_model,
            p_term = ifelse(nrow(focal_tidy) > 0, focal_tidy$p.value[1], NA_real_),
            best_shape_model = best_shape,
            PFAS = pfas,
            WBC_var = wbc_total,
            outcome_type = "total_count_nbinom2",
            n = nrow(df_model)
          )
      }
    }
  }
}

final_results <- bind_rows(results_list)

if (nrow(final_results) > 0) {
  final_results <- final_results %>%
    group_by(outcome_type, WBC_var) %>%
    mutate(FDR_p = p.adjust(p.value, method = "BH")) %>%
    ungroup()
}

write.csv(
  final_results,
  file = file.path(output_dir, "PFAS_WBC_mixed_model_results.csv"),
  row.names = FALSE
)

cat("Done: PFAS vs WBC mixed models\n")
