##########################################
# Figures from Scripts 3-5
# - Volcano plots (Scripts 3 & 4)
# - Forest plot of significant PFAS hits (Scripts 3 & 4)
# - Interaction coefficient plot (Script 5)
# - Predicted interaction plot for top Script 5 signal
##########################################

rm(list = ls())

library(tidyverse)
library(lme4)
library(lmerTest)

# ---------- Paths ----------
input_data <- "./input/2023_2024_2025_AllData_clean.csv"
file3 <- "./output/3_PFAS_vs_IgY/PFAS_IgY_mixed_model_results.csv"
file4 <- "./output/4_PFAS_vs_Clysis/PFAS_Clysis_mixed_model_results.csv"
file5 <- "./output/5_PFAS_interactions/PFAS_interactions.csv"
out_dir <- "./output/3_4_5_summary_figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---------- Naming ----------
vars_renamed <- c(
  "FOSA" = "FOSA",
  "FTSA8_2" = "8:2 FTS",
  "N_ETFOSAA" = "NEtFOSAA",
  "N_MEFOSAA" = "NMeFOSAA",
  "PFBA" = "PFBA",
  "PFDA" = "PFDA",
  "PFDOA" = "PFDoA",
  "PFDS" = "PFDS",
  "PFHPA" = "PFHPA",
  "PFHPS" = "PFHPS",
  "PFHXA" = "PFHxA",
  "PFHXS" = "PFHxS",
  "PFNA" = "PFNA",
  "PFNS" = "PFNS",
  "PFOA" = "PFOA",
  "PFOS" = "PFOS",
  "PFTEDA" = "PFTeA",
  "PFTRDA" = "PFTriA",
  "PFUNDA" = "PFUnA",
  "TOTAL_PFAS" = "Total PFAS",
  "age_next_gen" = "Next-Gen",
  "age_Legacy" = "Legacy",
  "func_amine" = "Amine",
  "func_carboxyl" = "Carboxyl",
  "func_sulfonate" = "Sulfonic",
  "func_phosphoric" = "Phosphoric",
  "chain_short" = "Short-chain",
  "chain_long" = "Long-chain"
)

to_display <- function(x) {
  ifelse(x %in% names(vars_renamed), vars_renamed[x], gsub("_", " ", x))
}

# ---------- Load outputs ----------
m3 <- read.csv(file3, stringsAsFactors = FALSE) %>%
  mutate(endpoint = "IgY")
m4 <- read.csv(file4, stringsAsFactors = FALSE) %>%
  mutate(endpoint = "Clysis")
m34 <- bind_rows(m3, m4)

m5 <- read.csv(file5, stringsAsFactors = FALSE)

# ---------- Figure A: Volcano (Scripts 3 & 4) ----------
volcano_df <- m34 %>%
  filter(!is.na(FDR_p), !is.na(estimate)) %>%
  mutate(
    predictor = to_display(pfas_var),
    neglogFDR = -log10(pmax(FDR_p, 1e-300)),
    sig = FDR_p < 0.05
  )

p_volcano <- ggplot(volcano_df, aes(x = estimate, y = neglogFDR)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.35) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "gray45", linewidth = 0.35) +
  geom_point(aes(color = sig), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c("FALSE" = "#5c5c5c", "TRUE" = "#d73027"), guide = "none") +
  facet_wrap(~endpoint, scales = "free_x") +
  labs(
    x = "Effect estimate (beta)",
    y = expression(-log[10]("FDR"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "Fig_34_volcano.png"), p_volcano, width = 9, height = 4.5, dpi = 600)
ggsave(file.path(out_dir, "Fig_34_volcano.tiff"), p_volcano, width = 9, height = 4.5, dpi = 600, compression = "lzw")

# ---------- Figure B: Forest of significant hits (Scripts 3 & 4) ----------
forest_df <- m34 %>%
  filter(!is.na(FDR_p), FDR_p < 0.05, !is.na(estimate), !is.na(std.error)) %>%
  mutate(
    predictor = to_display(pfas_var),
    ci_low = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error
  )

if (nrow(forest_df) > 0) {
  forest_df <- forest_df %>%
    group_by(endpoint) %>%
    mutate(predictor = factor(predictor, levels = predictor[order(estimate)])) %>%
    ungroup()
  
  p_forest <- ggplot(forest_df, aes(x = estimate, y = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.14, color = "gray45", linewidth = 0.4) +
    geom_point(color = "#1f78b4", size = 2) +
    facet_wrap(~endpoint, scales = "free_x") +
    labs(
      x = "Effect estimate (95% CI)",
      y = "PFAS predictor"
    ) +
    theme_classic(base_size = 11) +
    theme(
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir, "Fig_34_forest_significant.png"), p_forest, width = 9, height = 5.5, dpi = 600)
  ggsave(file.path(out_dir, "Fig_34_forest_significant.tiff"), p_forest, width = 9, height = 5.5, dpi = 600, compression = "lzw")
}

# ---------- Figure C: Interaction coefficient plot (Script 5) ----------
int_df <- m5 %>%
  filter(test_type == "interaction") %>%
  filter(grepl("PFAS_scaled|log_PFAS_scaled|PFAS_sq_scaled|exp_PFAS_scaled", term)) %>%
  filter(!is.na(estimate), !is.na(std.error), !is.na(FDR_p)) %>%
  mutate(
    endpoint = ifelse(immune_var == "log_igy", "IgY", "Clysis"),
    label = paste0(to_display(pfas_var), " | ", gsub("^.*:", "", term)),
    ci_low = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error,
    sig = FDR_p < 0.05
  )

if (nrow(int_df) > 0) {
  int_df <- int_df %>%
    group_by(endpoint) %>%
    arrange(estimate, .by_group = TRUE) %>%
    mutate(label = factor(label, levels = unique(label))) %>%
    ungroup()
  
  p_int <- ggplot(int_df, aes(x = estimate, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.35) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.12, color = "gray45", linewidth = 0.35) +
    geom_point(aes(color = sig), size = 1.8) +
    scale_color_manual(values = c("FALSE" = "#5c5c5c", "TRUE" = "#d73027"), guide = "none") +
    facet_wrap(~endpoint, scales = "free_x") +
    labs(
      x = "Interaction estimate (95% CI)",
      y = "PFAS | modifier"
    ) +
    theme_classic(base_size = 10.5) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 7.2),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(file.path(out_dir, "Fig_5_interaction_coefficients.png"), p_int, width = 10, height = 8, dpi = 600)
  ggsave(file.path(out_dir, "Fig_5_interaction_coefficients.tiff"), p_int, width = 10, height = 8, dpi = 600, compression = "lzw")
}

# ---------- Figure D: Predicted interaction plot for top Script 5 signal ----------
top_int <- m5 %>%
  filter(test_type == "interaction", !is.na(FDR_p)) %>%
  arrange(FDR_p, p.value) %>%
  slice(1)

if (nrow(top_int) == 1) {
  data_raw <- read.csv(input_data, stringsAsFactors = FALSE) %>%
    mutate(
      log_igy = log(igy + 1e-6),
      log_clysis = log(clysis + 1e-6)
    )
  pfas_var <- top_int$pfas_var[1]
  immune <- top_int$immune_var[1]
  best_shape_model <- top_int$best_shape_model[1]
  
  if (pfas_var %in% names(data_raw) && immune %in% names(data_raw)) {
    df_model <- data_raw %>%
      filter(!is.na(.data[[pfas_var]]), !is.na(.data[[immune]]), !is.na(DNA_sex), !is.na(age), !is.na(nest_no)) %>%
      mutate(PFAS_raw = .data[[pfas_var]])
    
    if (nrow(df_model) >= 10) {
      mu_raw <- mean(df_model$PFAS_raw)
      sd_raw <- sd(df_model$PFAS_raw)
      log_raw <- log(df_model$PFAS_raw + 1e-6)
      mu_log <- mean(log_raw)
      sd_log <- sd(log_raw)
      raw_sq <- df_model$PFAS_raw^2
      mu_sq <- mean(raw_sq)
      sd_sq <- sd(raw_sq)
      
      PFAS_scaled_tmp <- (df_model$PFAS_raw - mu_raw) / sd_raw
      exp_tmp <- exp(PFAS_scaled_tmp)
      mu_exp <- mean(exp_tmp)
      sd_exp <- sd(exp_tmp)
      
      df_model <- df_model %>%
        mutate(
          PFAS_scaled = (PFAS_raw - mu_raw) / sd_raw,
          log_PFAS_scaled = (log(PFAS_raw + 1e-6) - mu_log) / sd_log,
          PFAS_sq_scaled = (PFAS_raw^2 - mu_sq) / sd_sq,
          exp_PFAS_scaled = (exp(PFAS_scaled) - mu_exp) / sd_exp
        )
      
      if (best_shape_model == "linear") {
        focal_terms <- c("PFAS_scaled")
      } else if (best_shape_model == "log") {
        focal_terms <- c("log_PFAS_scaled")
      } else if (best_shape_model == "quadratic") {
        focal_terms <- c("PFAS_scaled", "PFAS_sq_scaled")
      } else if (best_shape_model == "exp") {
        focal_terms <- c("exp_PFAS_scaled")
      } else if (best_shape_model == "power") {
        focal_terms <- c("log_PFAS_scaled")
      } else {
        focal_terms <- c("PFAS_scaled")
      }
      
      df_model <- df_model %>% drop_na(all_of(focal_terms))
      
      form <- as.formula(
        paste0(immune, " ~ (", paste(focal_terms, collapse = " + "), ") * (DNA_sex + age) + (1|nest_no)")
      )
      fit <- try(lmer(form, data = df_model, REML = FALSE), silent = TRUE)
      
      if (!inherits(fit, "try-error")) {
        x_grid <- seq(min(df_model$PFAS_raw), max(df_model$PFAS_raw), length.out = 120)
        age_q <- as.numeric(quantile(df_model$age, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
        age_labs <- c("Age Q1", "Age Median", "Age Q3")
        sex_mode <- names(sort(table(df_model$DNA_sex), decreasing = TRUE))[1]
        
        pred_age <- expand.grid(
          PFAS_raw = x_grid,
          age = age_q,
          DNA_sex = sex_mode,
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        ) %>%
          mutate(line = factor(age, levels = age_q, labels = age_labs), panel = "Age moderation")
        
        pred_sex <- expand.grid(
          PFAS_raw = x_grid,
          age = median(df_model$age, na.rm = TRUE),
          DNA_sex = sort(unique(df_model$DNA_sex)),
          KEEP.OUT.ATTRS = FALSE,
          stringsAsFactors = FALSE
        ) %>%
          mutate(line = paste0("Sex ", DNA_sex), panel = "Sex moderation")
        
        pred_all <- bind_rows(pred_age, pred_sex) %>%
          mutate(
            PFAS_scaled = (PFAS_raw - mu_raw) / sd_raw,
            log_PFAS_scaled = (log(PFAS_raw + 1e-6) - mu_log) / sd_log,
            PFAS_sq_scaled = (PFAS_raw^2 - mu_sq) / sd_sq,
            exp_PFAS_scaled = (exp(PFAS_scaled) - mu_exp) / sd_exp
          )
        
        pred_all$fit <- predict(fit, newdata = pred_all, re.form = NA, allow.new.levels = TRUE)
        
        p_pred <- ggplot(pred_all, aes(x = PFAS_raw, y = fit, color = line)) +
          geom_line(linewidth = 1.0) +
          facet_wrap(~panel, scales = "free_y") +
          labs(
            x = paste0(to_display(pfas_var), " (raw)"),
            y = immune,
            color = NULL
          ) +
          theme_classic(base_size = 11) +
          theme(
            axis.title = element_text(face = "bold"),
            strip.text = element_text(face = "bold")
          )
        
        ggsave(file.path(out_dir, "Fig_5_top_interaction_predictions.png"), p_pred, width = 9, height = 4.6, dpi = 600)
        ggsave(file.path(out_dir, "Fig_5_top_interaction_predictions.tiff"), p_pred, width = 9, height = 4.6, dpi = 600, compression = "lzw")
      }
    }
  }
}

cat("Done: figures for Scripts 3-5\n")
