rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyr)

# ========= SETTINGS =========

if (!dir.exists("./output/2_clysis_pfas_relationship"))
  dir.create("./output/2_clysis_pfas_relationship")

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

chain_sums <- long_pf %>%
  filter(!is.na(chain_cat)) %>%
  group_by(bird_id, chain_cat) %>%
  summarize(sum_val = sum(conc), .groups = "drop") %>%
  pivot_wider(names_from = chain_cat,
              values_from = sum_val,
              values_fill = 0,
              names_prefix = "chain_")

func_sums <- long_pf %>%
  filter(!is.na(func_cat)) %>%
  group_by(bird_id, func_cat) %>%
  summarize(sum_val = sum(conc), .groups = "drop") %>%
  pivot_wider(names_from = func_cat,
              values_from = sum_val,
              values_fill = 0,
              names_prefix = "func_")

age_sums <- long_pf %>%
  filter(!is.na(pfas_age)) %>%
  group_by(bird_id, pfas_age) %>%
  summarize(sum_val = sum(conc), .groups = "drop") %>%
  pivot_wider(names_from = pfas_age,
              values_from = sum_val,
              values_fill = 0,
              names_prefix = "age_")

data_model <- data_clean %>%
  left_join(chain_sums, by = "bird_id") %>%
  left_join(func_sums,  by = "bird_id") %>%
  left_join(age_sums,   by = "bird_id")

# Replace NA sums with 0
sum_cols <- grep("^(chain_|func_|age_)", colnames(data_model), value = TRUE)
data_model[sum_cols][is.na(data_model[sum_cols])] <- 0

pfas_vars <- unique(c(pfas_present, sum_cols))

# ========= MODELING FUNCTION =========

run_shape_models <- function(response_var, response_label) {
  
  results_list <- list()
  
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
    
    if (nrow(df) < 6) next
    
    df <- df %>%
      mutate(
        x = .data[[pfas_var]],
        log_x = log(x + 1e-6),
        x2 = x^2,
        x_adj = x + 1e-6
      )
    
    # ---- Fit models ----
    
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
    
    # ---- Prediction grid ----
    
    new_x <- seq(min(df$x), max(df$x), length.out = 200)
    pred_df <- data.frame(
      x = new_x,
      log_x = log(new_x + 1e-6),
      x2 = new_x^2,
      x_adj = new_x + 1e-6
    )
    
    pred_df$fit <- predict(best_model, newdata = pred_df)
    
    # ---- Plot ----
    
    p <- ggplot(df, aes(x = x, y = y)) +
      geom_point(alpha = 0.6) +
      geom_line(data = pred_df,
                aes(x = x, y = fit),
                color = "blue",
                linewidth = 1) +
      labs(title = paste("Best Fit:", best_model_name, "-", pfas_var),
           subtitle = response_label,
           x = pfas_var,
           y = response_label) +
      theme_minimal()
    
    ggsave(
      filename = paste0("./output/2_clysis_pfas_relationship/",
                        pfas_var, "_", response_label, "_shape.jpeg"),
      plot = p,
      width = 7, height = 5, dpi = 300
    )
    
    results_list[[paste(pfas_var, response_label, sep = "_")]] <-
      data.frame(
        PFAS = pfas_var,
        Response = response_label,
        Best_Model = best_model_name,
        AIC = min(model_aics)
      )
  }
  
  return(results_list)
}

# ========= RUN BOTH RESPONSE TYPES =========

results_raw  <- run_shape_models("clysis", "clysis_raw")
results_log  <- run_shape_models("clysis", "log_clysis")

# ========= SAVE SUMMARY =========

all_results <- do.call(rbind, c(results_raw, results_log))

write.csv(all_results,
          "./output/2_clysis_pfas_relationship/clysis_best_shape_summary.csv",
          row.names = FALSE)
