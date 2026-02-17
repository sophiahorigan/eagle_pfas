rm(list = ls())

library(nlme)
library(ggplot2)
library(dplyr)
library(tidyr)

# ========== SETTINGS ==========
add_CI <- FALSE
setwd("~/Desktop/EagleStats/")

if (!dir.exists("./output/PFAS_model_plots")) dir.create("./output/PFAS_model_plots")

# ========== LOAD DATA ==========
data <- read.csv("./input/2023_2024_2025_AllData.csv", stringsAsFactors = FALSE)
data_clean <- data %>% filter(!is.na(igy), !is.na(nest_no))

# ========== LOAD PFAS KEY ==========
key <- read.csv("./input/PFAS_subtypes.csv", stringsAsFactors = FALSE)

setdiff(key$variable, names(df))

# Ensure the key has a matching analyte name column called 'pfas_analyte'
# (adjust if your key uses a different column name)
if (!"pfas_analyte" %in% colnames(key)) stop("PFAS key must have column 'pfas_analyte'")

# ========== IDENTIFY PFAS COLUMNS IN MAIN DATA ==========
pfas_present <- intersect(key$pfas_analyte, colnames(data_clean))
if (length(pfas_present) == 0) stop("No PFAS analyte columns found in main data matching key$pfas_analyte")

# ========== WIDE -> LONG of PFAS measurements ==========
# We assume PFAS columns are numeric concentrations; pivot them long for summing by category.
long_pf <- data_clean %>%
  select(bird_id, all_of(pfas_present)) %>%
  pivot_longer(cols = all_of(pfas_present),
               names_to = "pfas_analyte",
               values_to = "conc") %>%
  # keep only measured values (NA or non-numeric removed)
  filter(!is.na(conc)) %>%
  left_join(key, by = "pfas_analyte")

# ========== BUILD CATEGORY SUMS (chain_cat, func_cat, pfas_age) ==========
chain_sums <- long_pf %>%
  filter(!is.na(chain_cat)) %>%
  group_by(bird_id, chain_cat) %>%
  summarize(chain_sum = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = chain_cat, values_from = chain_sum,
              values_fill = 0, names_prefix = "chain_")

func_sums <- long_pf %>%
  filter(!is.na(func_cat)) %>%
  group_by(bird_id, func_cat) %>%
  summarize(func_sum = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = func_cat, values_from = func_sum,
              values_fill = 0, names_prefix = "func_")

age_sums <- long_pf %>%
  filter(!is.na(pfas_age)) %>%
  group_by(bird_id, pfas_age) %>%
  summarize(age_sum = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pfas_age, values_from = age_sum,
              values_fill = 0, names_prefix = "age_")

# If any of these are empty (no categories present), create empty tibbles to avoid join failures
if (nrow(chain_sums) == 0) chain_sums <- tibble(bird_id = character(0))
if (nrow(func_sums) == 0)  func_sums  <- tibble(bird_id = character(0))
if (nrow(age_sums) == 0)   age_sums   <- tibble(bird_id = character(0))

# ========== MERGE CATEGORY SUMS BACK INTO MAIN DATA ==========
data_model <- data_clean %>%
  left_join(chain_sums, by = "bird_id") %>%
  left_join(func_sums,  by = "bird_id") %>%
  left_join(age_sums,   by = "bird_id")

# Replace any NA in the new sum-columns with 0 (no measured analytes in that category)
sum_cols <- grep("^(chain_|func_|age_)", colnames(data_model), value = TRUE)
if (length(sum_cols) > 0) data_model[sum_cols][is.na(data_model[sum_cols])] <- 0

# ========== BUILD LIST OF PREDICTORS ==========
pfas_individual <- pfas_present  # analytes present in the dataframe
pfas_chain_vars <- setdiff(colnames(chain_sums), "bird_id")
pfas_func_vars  <- setdiff(colnames(func_sums), "bird_id")
pfas_age_vars   <- setdiff(colnames(age_sums), "bird_id")

pfas_vars <- unique(c(pfas_individual, pfas_chain_vars, pfas_func_vars, pfas_age_vars))

cat("Modeling these predictors (count):", length(pfas_vars), "\n")
cat("Examples:", paste(head(pfas_vars, 10), collapse = ", "), "\n")

# ========== STORAGE ==========
all_aic_results <- list()
all_best_model_stats <- list()

# Helper: safe model fit checker
is_model_ok <- function(obj) inherits(obj, c("lme", "nlme"))

# ========== MODELING LOOP ==========
for (pfas_var in pfas_vars) {
  cat("\n==== IgY models for:", pfas_var, "====\n")
  
  # Build df for modeling
  if (!pfas_var %in% colnames(data_model)) {
    cat("  SKIP (variable not in data_model):", pfas_var, "\n")
    next
  }
  
  df <- data_model %>%
    select(igy, nest_no, all_of(pfas_var)) %>%
    filter(!is.na(.data[[pfas_var]]))
  
  if (nrow(df) < 6) { # small sample guard
    cat("  SKIP (too few observations):", nrow(df), "\n")
    next
  }
  
  df <- df %>% mutate(
    xvar_val = .data[[pfas_var]],
    log_x = log(xvar_val + 1e-6),
    x2 = xvar_val^2
  )
  
  # Fit candidate models (wrap in try/catch)
  models <- list(
    linear = try(lme(igy ~ xvar_val, random = ~1 | nest_no, data = df), silent = TRUE),
    log = try(lme(igy ~ log_x, random = ~1 | nest_no, data = df), silent = TRUE),
    quadratic = try(lme(igy ~ xvar_val + x2, random = ~1 | nest_no, data = df), silent = TRUE),
    exp = try(nlme(igy ~ B0 * exp(-B1 * xvar_val),
                   fixed = B0 + B1 ~ 1, random = B0 ~ 1 | nest_no,
                   start = c(B0 = 0.5, B1 = 0.01), data = df), silent = TRUE),
    power = try({
      df_power <- df %>% mutate(xvar_val_adj = xvar_val + 1e-6)
      nlme(igy ~ B0 * xvar_val_adj^(-B1),
           fixed = B0 + B1 ~ 1, random = B0 ~ 1 | nest_no,
           start = c(B0 = 0.5, B1 = 0.1), data = df_power)
    }, silent = TRUE)
  )
  
  valid_models <- models[sapply(models, is_model_ok)]
  if (length(valid_models) == 0) {
    cat("  No valid models for", pfas_var, "\n")
    next
  }
  
  model_aics <- sapply(valid_models, AIC)
  best_model_name <- names(which.min(model_aics))
  best_model <- valid_models[[best_model_name]]
  
  aic_df <- data.frame(PFAS = pfas_var, Model = names(model_aics), AIC = as.numeric(model_aics))
  all_aic_results[[pfas_var]] <- aic_df
  
  preds <- predict(best_model, level = 0)
  residuals <- df$igy - preds
  n <- nrow(df)
  p <- length(fixef(best_model))
  RSE <- sqrt(sum(residuals^2) / (n - p))
  
  p_values <- tryCatch({
    sm <- summary(best_model)
    if ("tTable" %in% names(sm)) sm$tTable[, "p-value"] else rep(NA, length(fixef(best_model)))
  }, error = function(e) rep(NA, length(fixef(best_model))))
  
  best_model_stats <- data.frame(
    PFAS = pfas_var,
    Best_Model = best_model_name,
    Residual_Std_Error = RSE,
    Fixed_Effect = names(p_values),
    P_value = as.numeric(p_values),
    stringsAsFactors = FALSE
  )
  all_best_model_stats[[pfas_var]] <- best_model_stats
  
  # Prediction grid
  new_x <- seq(min(df$xvar_val, na.rm = TRUE), max(df$xvar_val, na.rm = TRUE), length.out = 200)
  pred_df <- data.frame(
    xvar_val = new_x,
    log_x = log(new_x + 1e-6),
    x2 = new_x^2,
    xvar_val_adj = new_x + 1e-6,
    nest_no = df$nest_no[1]
  )
  
  pred_df$fit <- predict(best_model, newdata = pred_df, level = 0)
  
  # Overlay all valid models
  all_fit_df <- do.call(rbind, lapply(names(valid_models), function(mname) {
    mod <- valid_models[[mname]]
    tmp <- pred_df
    if (mname == "power") tmp$xvar_val <- tmp$xvar_val_adj
    yhat <- tryCatch(predict(mod, newdata = tmp, level = 0), error = function(e) rep(NA, nrow(tmp)))
    data.frame(x = new_x, y = yhat, model = mname)
  }))
  
  p_all <- ggplot(df, aes(x = xvar_val, y = igy)) +
    geom_point(alpha = 0.6) +
    geom_line(data = all_fit_df, aes(x = x, y = y, color = model), linewidth = 0.8, na.rm = TRUE) +
    labs(title = paste("IgY vs", pfas_var, "(All Models)"), x = pfas_var, y = "IgY") +
    theme_minimal()
  
  ggsave(filename = paste0("./output/PFAS_model_plots/", pfas_var, "_all_models.jpeg"),
         plot = p_all, width = 7, height = 5, dpi = 300)
  
  # Best model plot (with optional CI)
  p_best <- ggplot(df, aes(x = xvar_val, y = igy)) + geom_point(alpha = 0.6)
  if (add_CI) {
    # CI computation for nlme/lme is more involved; leave disabled unless you want this implemented
    # (we keep add_CI in code as a placeholder)
  }
  p_best <- p_best +
    geom_line(data = pred_df, aes(x = xvar_val, y = fit), inherit.aes = FALSE, color = "blue") +
    labs(title = paste("Best Model:", best_model_name, "-", pfas_var),
         subtitle = paste0("RSE = ", round(RSE, 3)),
         x = pfas_var, y = "IgY") +
    theme_minimal()
  
  ggsave(filename = paste0("./output/PFAS_model_plots/", pfas_var, "_best_model.jpeg"),
         plot = p_best, width = 7, height = 5, dpi = 300)
}

# ========== SAVE AIC + MODEL STATS ==========
if (length(all_aic_results) > 0) {
  aic_df <- do.call(rbind, all_aic_results)
  write.csv(aic_df, "./output/PFAS_model_plots/AIC_values_IgY.csv", row.names = FALSE)
}

if (length(all_best_model_stats) > 0) {
  best_model_df <- do.call(rbind, all_best_model_stats)
  write.csv(best_model_df, "./output/PFAS_model_plots/best_model_stats.csv", row.names = FALSE)
}

# ========== RESIDUAL PLOTS FOR BEST MODELS ==========
for (pfas_var in pfas_vars) {
  best_stats <- all_best_model_stats[[pfas_var]]
  if (is.null(best_stats)) next
  
  best_model_type <- unique(best_stats$Best_Model)
  
  if (!pfas_var %in% colnames(data_model)) next
  
  df <- data_model %>%
    select(igy, nest_no, all_of(pfas_var)) %>%
    filter(!is.na(.data[[pfas_var]])) %>%
    mutate(xvar_val = .data[[pfas_var]],
           log_x = log(xvar_val + 1e-6),
           x2 = xvar_val^2)
  
  if (nrow(df) < 6) next
  if (best_model_type == "power") df <- df %>% mutate(xvar_val_adj = xvar_val + 1e-6)
  
  # Re-fit the chosen model for residuals
  best_model <- switch(best_model_type,
                       linear = try(lme(igy ~ xvar_val, random = ~1 | nest_no, data = df), silent = TRUE),
                       log = try(lme(igy ~ log_x, random = ~1 | nest_no, data = df), silent = TRUE),
                       quadratic = try(lme(igy ~ xvar_val + x2, random = ~1 | nest_no, data = df), silent = TRUE),
                       exp = try(nlme(igy ~ B0 * exp(-B1 * xvar_val),
                                      fixed = B0 + B1 ~ 1, random = B0 ~ 1 | nest_no,
                                      start = c(B0 = 0.5, B1 = 0.01), data = df), silent = TRUE),
                       power = try(nlme(igy ~ B0 * xvar_val_adj^(-B1),
                                        fixed = B0 + B1 ~ 1, random = B0 ~ 1 | nest_no,
                                        start = c(B0 = 0.5, B1 = 0.1), data = df), silent = TRUE),
                       NULL
  )
  
  if (!is_model_ok(best_model)) next
  
  df$fitted <- predict(best_model, level = 0)
  df$residual <- df$igy - df$fitted
  
  p_resid <- ggplot(df, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Residuals vs Fitted:", pfas_var),
         x = "Fitted values", y = "Residuals") +
    theme_minimal()
  
  ggsave(filename = paste0("./output/PFAS_model_plots/", pfas_var, "_residuals.jpeg"),
         plot = p_resid, width = 6, height = 5, dpi = 300)
}

