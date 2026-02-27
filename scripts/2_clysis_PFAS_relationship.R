rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)

# ========= SETTINGS =========
out_dir <- "./output/2_clysis_pfas_relationship"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
raw_plot_dir <- file.path(out_dir, "raw")
log_plot_dir <- file.path(out_dir, "log")
if (!dir.exists(raw_plot_dir)) dir.create(raw_plot_dir, recursive = TRUE)
if (!dir.exists(log_plot_dir)) dir.create(log_plot_dir, recursive = TRUE)

# ========= LOAD DATA =========
data <- read.csv("./input/2023_2024_2025_AllData.csv", stringsAsFactors = FALSE)
data_clean <- data %>% filter(!is.na(clysis))

key <- read.csv("./input/PFAS_subtypes.csv", stringsAsFactors = FALSE)

# ========= IDENTIFY PFAS =========
pfas_present <- intersect(key$pfas_analyte, colnames(data_clean))

# ========= BUILD CATEGORY SUMS =========
long_pf <- data_clean %>%
  select(bird_id, all_of(pfas_present)) %>%
  pivot_longer(cols = all_of(pfas_present),
               names_to = "pfas_analyte",
               values_to = "conc") %>%
  filter(!is.na(conc)) %>%
  left_join(key, by = "pfas_analyte")

make_sum_table <- function(df, group_var, prefix) {
  df %>%
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

data_model <- data_clean %>%
  left_join(chain_sums, by = "bird_id") %>%
  left_join(func_sums,  by = "bird_id") %>%
  left_join(age_sums,   by = "bird_id")

sum_cols <- grep("^(chain_|func_|age_)", colnames(data_model), value = TRUE)
data_model[sum_cols][is.na(data_model[sum_cols])] <- 0

pfas_matrix <- data_model[, pfas_present, drop = FALSE]
data_model$TOTAL_PFAS <- rowSums(pfas_matrix, na.rm = TRUE)
data_model$TOTAL_PFAS[rowSums(!is.na(pfas_matrix)) == 0] <- NA_real_

pfas_vars <- unique(c(pfas_present, "TOTAL_PFAS", sum_cols))

# ========= MODELING FUNCTION =========

run_shape_models <- function(response_var, response_label) {
  
  results_list <- list()
  plot_list <- list()
  
  for (pfas_var in pfas_vars) {
    
    if (!pfas_var %in% colnames(data_model)) next
    
    df <- data_model %>%
      select(all_of(response_var), all_of(pfas_var)) %>%
      filter(!is.na(.data[[pfas_var]]))
    
    colnames(df)[1] <- "y"
    
    if (response_label == "log_clysis") {
      df <- df %>% filter(y > 0)
      df$y <- log(df$y)
    }
    
    if (nrow(df) < 8) next
    
    df <- df %>%
      mutate(
        x = .data[[pfas_var]],
        log_x = log(x + 1e-6),
        x2 = x^2,
        x_adj = x + 1e-6
      )
    
    # ---- Fit candidate models ----
    models <- list(
      linear = try(lm(y ~ x, data = df), silent = TRUE),
      log = try(lm(y ~ log_x, data = df), silent = TRUE),
      quadratic = try(lm(y ~ x + x2, data = df), silent = TRUE),
      exp = try(nls(y ~ a * exp(b * x),
                    start = list(a = mean(df$y), b = -0.01),
                    data = df), silent = TRUE),
      power = try(nls(y ~ a * x_adj^b,
                      start = list(a = mean(df$y), b = -0.5),
                      data = df), silent = TRUE)
    )
    
    valid_models <- models[!sapply(models, inherits, "try-error")]
    if (length(valid_models) == 0) next
    
    model_aics <- sapply(valid_models, AIC)
    best_model_name <- names(which.min(model_aics))
    best_model <- valid_models[[best_model_name]]
    
    # ---- Null model comparison ----
    null_model <- lm(y ~ 1, data = df)
    
    if (inherits(best_model, "lm")) {
      
      model_test <- anova(null_model, best_model)
      p_value <- model_test$`Pr(>F)`[2]
      r2 <- summary(best_model)$r.squared
      adj_r2 <- summary(best_model)$adj.r.squared
      
    } else if (inherits(best_model, "nls")) {
      
      delta_aic <- AIC(null_model) - AIC(best_model)
      lr_stat <- 2 * delta_aic
      df_diff <- length(coef(best_model))
      p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
      
      r2 <- NA
      adj_r2 <- NA
    }
    
    # ---- Prediction grid ----
    new_x <- seq(min(df$x), max(df$x), length.out = 200)
    
    pred_df <- data.frame(
      x = new_x,
      log_x = log(new_x + 1e-6),
      x2 = new_x^2,
      x_adj = new_x + 1e-6
    )
    
    pred_df$fit <- predict(best_model, newdata = pred_df)
    
    # ---- Determine Direction ----
    y_min_pred <- predict(best_model,
                          newdata = data.frame(
                            x = min(df$x),
                            log_x = log(min(df$x) + 1e-6),
                            x2 = min(df$x)^2,
                            x_adj = min(df$x) + 1e-6
                          ))
    
    y_max_pred <- predict(best_model,
                          newdata = data.frame(
                            x = max(df$x),
                            log_x = log(max(df$x) + 1e-6),
                            x2 = max(df$x)^2,
                            x_adj = max(df$x) + 1e-6
                          ))
    
    delta_y <- as.numeric(y_max_pred - y_min_pred)
    slope_global <- delta_y / (max(df$x) - min(df$x))
    
    direction <- ifelse(delta_y > 0, "Increasing",
                        ifelse(delta_y < 0, "Decreasing", "Flat"))
    
    row_key <- paste(pfas_var, response_label, sep = "_")
    
    results_list[[row_key]] <-
      data.frame(
        PFAS = pfas_var,
        Response = response_label,
        Best_Model = best_model_name,
        AIC = min(model_aics),
        R2 = r2,
        Adj_R2 = adj_r2,
        P_value = p_value,
        N = nrow(df),
        Slope_Global = slope_global,
        Delta_Y = delta_y,
        Direction = direction
      )
    
    plot_list[[row_key]] <- list(
      df = df,
      pred_df = pred_df,
      best_model_name = best_model_name,
      pfas_var = pfas_var,
      response_label = response_label,
      file_path = file.path(
        ifelse(grepl("^log_", response_label), log_plot_dir, raw_plot_dir),
        paste0(pfas_var, "_", response_label, "_shape.jpeg")
      )
    )
  }
  
  return(list(results = results_list, plots = plot_list))
}

# ========= RUN =========
results_raw <- run_shape_models("clysis", "clysis_raw")
results_log <- run_shape_models("clysis", "log_clysis")

all_results <- do.call(rbind, c(results_raw$results, results_log$results))

# ========= MULTIPLE TEST CORRECTION =========
all_results$FDR_P <- p.adjust(all_results$P_value, method = "BH")

# ========= PLOT WITH SIGNIFICANCE COLORS =========
all_plots <- c(results_raw$plots, results_log$plots)

for (row_key in names(all_plots)) {
  plot_info <- all_plots[[row_key]]
  result_row <- all_results[row_key, , drop = FALSE]
  
  is_significant <- isTRUE(result_row$P_value < 0.05) || isTRUE(result_row$FDR_P < 0.05)
  line_color <- ifelse(is_significant, "blue", "gray80")
  
  p <- ggplot(plot_info$df, aes(x = x, y = y)) +
    geom_point(alpha = 0.6) +
    geom_line(data = plot_info$pred_df,
              aes(x = x, y = fit),
              color = line_color,
              linewidth = 1) +
    labs(title = paste("Best Fit:", plot_info$best_model_name, "-", plot_info$pfas_var),
         subtitle = paste(plot_info$response_label,
                          "| p =", signif(result_row$P_value, 3),
                          "| FDR =", signif(result_row$FDR_P, 3)),
         x = plot_info$pfas_var,
         y = plot_info$response_label) +
    theme_minimal()
  
  ggsave(
    filename = plot_info$file_path,
    plot = p,
    width = 7, height = 5, dpi = 300
  )
}

# ========= SAVE =========
write.csv(all_results,
          file.path(out_dir, "clysis_best_shape_summary.csv"),
          row.names = FALSE)

print("Clysis modeling complete.")
