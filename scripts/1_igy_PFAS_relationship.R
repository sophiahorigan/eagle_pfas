rm(list = ls())

library(nlme)
library(ggplot2)
library(dplyr)
library(tidyr)

# ========== SETTINGS ==========
add_CI <- FALSE
setwd("~/Desktop/EagleStats/")

if (!dir.exists("./output/PFAS_model_plots")) dir.create("./output/PFAS_model_plots", recursive = TRUE)

# ========== LOAD DATA ==========
data <- read.csv("./input/2023_2024_2025_AllData.csv", stringsAsFactors = FALSE)
data_clean <- data %>% filter(!is.na(igy), !is.na(nest_no))

# ========== LOAD PFAS KEY ==========
key <- read.csv("./input/PFAS_subtypes.csv", stringsAsFactors = FALSE)
if (!"pfas_analyte" %in% colnames(key)) stop("PFAS key must have column 'pfas_analyte'")

# ========== IDENTIFY PFAS COLUMNS ==========
pfas_present <- intersect(key$pfas_analyte, colnames(data_clean))
if (length(pfas_present) == 0) stop("No PFAS analytes found in data")

# ========== LONG FORMAT ==========
long_pf <- data_clean %>%
  select(bird_id, all_of(pfas_present)) %>%
  pivot_longer(cols = all_of(pfas_present),
               names_to = "pfas_analyte",
               values_to = "conc") %>%
  filter(!is.na(conc)) %>%
  left_join(key, by = "pfas_analyte")

# ========== CATEGORY SUMS ==========
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

if (nrow(chain_sums) == 0) chain_sums <- tibble(bird_id = character(0))
if (nrow(func_sums)  == 0) func_sums  <- tibble(bird_id = character(0))
if (nrow(age_sums)   == 0) age_sums   <- tibble(bird_id = character(0))

# ========== MERGE BACK ==========
data_model <- data_clean %>%
  left_join(chain_sums, by = "bird_id") %>%
  left_join(func_sums,  by = "bird_id") %>%
  left_join(age_sums,   by = "bird_id")

sum_cols <- grep("^(chain_|func_|age_)", colnames(data_model), value = TRUE)
data_model[sum_cols][is.na(data_model[sum_cols])] <- 0

# ========== PFAS VARIABLES ==========
pfas_individual <- pfas_present
pfas_chain_vars <- setdiff(colnames(chain_sums), "bird_id")
pfas_func_vars  <- setdiff(colnames(func_sums), "bird_id")
pfas_age_vars   <- setdiff(colnames(age_sums), "bird_id")

pfas_vars <- unique(c(pfas_individual, pfas_chain_vars, pfas_func_vars, pfas_age_vars))
cat("Modeling", length(pfas_vars), "PFAS predictors\n")

# ========== STORAGE ==========
all_aic_results <- list()
all_best_model_stats <- list()

is_model_ok <- function(obj) inherits(obj, c("lme", "nlme"))

# ========== MODEL LOOP (LOG PFAS) ==========
for (pfas_var in pfas_vars) {
  cat("\n==== Modeling", pfas_var, "====\n")
  
  if (!pfas_var %in% colnames(data_model)) next
  
  df <- data_model %>%
    select(igy, nest_no, all_of(pfas_var)) %>%
    filter(!is.na(.data[[pfas_var]]))
  
  if (nrow(df) < 6) {
    cat(" SKIP: n < 6\n")
    next
  }
  
  df <- df %>% mutate(
    x_raw = .data[[pfas_var]],
    log_pfas = log(x_raw + 1e-6),
    log_pfas2 = log_pfas^2
  )
  
  # Candidate models on log PFAS
  models <- list(
    linear = try(lme(igy ~ log_pfas, random = ~1 | nest_no, data = df), silent = TRUE),
    quadratic = try(lme(igy ~ log_pfas + log_pfas2, random = ~1 | nest_no, data = df), silent = TRUE),
    exp = try(nlme(igy ~ B0 * exp(-B1 * log_pfas),
                   fixed = B0 + B1 ~ 1, random = B0 ~ 1 | nest_no,
                   start = c(B0 = 0.5, B1 = 0.1), data = df), silent = TRUE)
  )
  
  valid_models <- models[sapply(models, is_model_ok)]
  if (length(valid_models) == 0) {
    cat(" No valid models\n")
    next
  }
  
  model_aics <- sapply(valid_models, AIC)
  best_model_name <- names(which.min(model_aics))
  best_model <- valid_models[[best_model_name]]
  
  # Save AICs
  all_aic_results[[pfas_var]] <- data.frame(PFAS = pfas_var, Model = names(model_aics), AIC = model_aics)
  
  # Residual SE
  preds <- predict(best_model, level = 0)
  residuals <- df$igy - preds
  n <- nrow(df)
  p <- length(fixef(best_model))
  RSE <- sqrt(sum(residuals^2) / (n - p))
  
  # p-values
  p_values <- tryCatch({
    sm <- summary(best_model)
    if ("tTable" %in% names(sm)) sm$tTable[, "p-value"] else rep(NA, length(fixef(best_model)))
  }, error = function(e) rep(NA, length(fixef(best_model))))
  
  all_best_model_stats[[pfas_var]] <- data.frame(
    PFAS = pfas_var,
    Best_Model = best_model_name,
    Residual_Std_Error = RSE,
    Fixed_Effect = names(p_values),
    P_value = as.numeric(p_values)
  )
  
  # Prediction grid in log-space
  new_log <- seq(min(df$log_pfas), max(df$log_pfas), length.out = 200)
  pred_df <- data.frame(
    log_pfas = new_log,
    log_pfas2 = new_log^2,
    nest_no = df$nest_no[1]
  )
  pred_df$fit <- predict(best_model, newdata = pred_df, level = 0)
  pred_df$pfas_raw <- exp(pred_df$log_pfas)
  
  # Overlay all models
  all_fit_df <- do.call(rbind, lapply(names(valid_models), function(m) {
    tmp <- pred_df
    yhat <- tryCatch(predict(valid_models[[m]], newdata = tmp, level = 0),
                     error = function(e) rep(NA, nrow(tmp)))
    data.frame(x = tmp$pfas_raw, y = yhat, model = m)
  }))
  
  # ========== PLOTS ==========
  p_all <- ggplot(df, aes(x = x_raw, y = igy)) +
    geom_point(alpha = 0.6) +
    geom_line(data = all_fit_df, aes(x = x, y = y, color = model)) +
    scale_x_log10() +
    theme_minimal() +
    labs(title = paste("IgY vs", pfas_var, "(log PFAS models)"),
         x = paste0(pfas_var, " (log scale)"), y = "IgY")
  
  ggsave(paste0("./output/PFAS_model_plots/", pfas_var, "_all_models_log.jpeg"),
         p_all, width = 7, height = 5, dpi = 300)
  
  p_best <- ggplot(df, aes(x = x_raw, y = igy)) +
    geom_point(alpha = 0.6) +
    geom_line(data = pred_df, aes(x = pfas_raw, y = fit), inherit.aes = FALSE, color = "blue") +
    scale_x_log10() +
    theme_minimal() +
    labs(title = paste("Best model:", best_model_name, "-", pfas_var),
         subtitle = paste0("RSE = ", round(RSE, 3)),
         x = paste0(pfas_var, " (log scale)"), y = "IgY")
  
  ggsave(paste0("./output/PFAS_model_plots/", pfas_var, "_best_model_log.jpeg"),
         p_best, width = 7, height = 5, dpi = 300)
}

# ========== SAVE TABLES ==========
if (length(all_aic_results) > 0) {
  write.csv(do.call(rbind, all_aic_results),
            "./output/PFAS_model_plots/AIC_values_IgY_log.csv", row.names = FALSE)
}

if (length(all_best_model_stats) > 0) {
  write.csv(do.call(rbind, all_best_model_stats),
            "./output/PFAS_model_plots/best_model_stats_log.csv", row.names = FALSE)
}

# ========== RESIDUAL DIAGNOSTICS ==========
for (pfas_var in pfas_vars) {
  best_stats <- all_best_model_stats[[pfas_var]]
  if (is.null(best_stats)) next
  best_model_type <- unique(best_stats$Best_Model)
  if (!pfas_var %in% colnames(data_model)) next
  
  df <- data_model %>%
    select(igy, nest_no, all_of(pfas_var)) %>%
    filter(!is.na(.data[[pfas_var]])) %>%
    mutate(x_raw = .data[[pfas_var]],
           log_pfas = log(x_raw + 1e-6),
           log_pfas2 = log_pfas^2)
  
  if (nrow(df) < 6) next
  
  best_model <- switch(best_model_type,
                       linear = lme(igy ~ log_pfas, random = ~1 | nest_no, data = df),
                       quadratic = lme(igy ~ log_pfas + log_pfas2, random = ~1 | nest_no, data = df),
                       exp = nlme(igy ~ B0 * exp(-B1 * log_pfas),
                                  fixed = B0 + B1 ~ 1, random = B0 ~ 1 | nest_no,
                                  start = c(B0 = 0.5, B1 = 0.1), data = df),
                       NULL)
  if (!is_model_ok(best_model)) next
  
  df$fitted <- predict(best_model, level = 0)
  df$residual <- df$igy - df$fitted
  
  p_resid <- ggplot(df, aes(x = fitted, y = residual)) +
    geom_point(alpha = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = paste("Residuals vs Fitted:", pfas_var))
  
  ggsave(paste0("./output/PFAS_model_plots/", pfas_var, "_residuals_log.jpeg"),
         p_resid, width = 6, height = 5, dpi = 300)
}
