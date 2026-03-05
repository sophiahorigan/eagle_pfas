##########################################
# Figure: Key take-homes from Script 1 (IgY)
# Publication-ready panel figure
##########################################

rm(list = ls())

library(tidyverse)

# ---------- Settings ----------
input_file <- "./input/2023_2024_2025_AllData_clean.csv"
key_file <- "./input/PFAS_subtypes.csv"
shape_file <- "./output/1_igy_pfas_relationship/best_shape_summary.csv"
out_dir <- "./output/1_igy_pfas_relationship/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# How many panels to show
n_panels <- 9

# If FALSE, prioritize individual PFAS analytes over grouped sums
allow_grouped_predictors <- TRUE

# Force-include these predictors in the main panel figure (if available)
forced_predictors <- c("TOTAL_PFAS", "func_sulfonate", "age_next_gen")

# Display-name convention
vars_renamed <- c(
  "FOSA"       = "FOSA",
  "FTSA8_2"    = "8:2 FTS",
  "N_ETFOSAA"  = "NEtFOSAA",
  "N_MEFOSAA"  = "NMeFOSAA",
  "PFBA"       = "PFBA",
  "PFDA"       = "PFDA",
  "PFDOA"      = "PFDoA",
  "PFDS"       = "PFDS",
  "PFHPA"      = "PFHPA",
  "PFHPS"      = "PFHPS",
  "PFHXA"      = "PFHxA",
  "PFHXS"      = "PFHxS",
  "PFNA"       = "PFNA",
  "PFNS"       = "PFNS",
  "PFOA"       = "PFOA",
  "PFOS"       = "PFOS",
  "PFTEDA"     = "PFTeA",
  "PFTRDA"     = "PFTriA",
  "PFUNDA"     = "PFUnA",
  "TOTAL_PFAS" = "Total PFAS",
  "age_next_gen" = "Next-Gen",
  "func_amine" = "Amine",
  "func_sulfonate" = "Sulfonic"
)

# ---------- Load data ----------
data <- read.csv(input_file, stringsAsFactors = FALSE)
key <- read.csv(key_file, stringsAsFactors = FALSE)
shape <- read.csv(shape_file, stringsAsFactors = FALSE)

data_clean <- data %>% filter(!is.na(igy), igy > 0)

# ---------- Build predictors exactly like Script 1 ----------
pfas_present <- intersect(key$pfas_analyte, colnames(data_clean))

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

# ---------- Select key relationships ----------
shape_log <- shape %>%
  filter(Response == "log_IgY") %>%
  mutate(
    predictor_type = case_when(
      PFAS %in% pfas_present ~ "analyte",
      PFAS == "TOTAL_PFAS" ~ "total",
      grepl("^(chain_|func_|age_)", PFAS) ~ "grouped",
      TRUE ~ "other"
    )
  ) %>%
  arrange(FDR_P, P_value)

sig_tbl <- shape_log %>%
  filter(!is.na(FDR_P), FDR_P < 0.05)

if (nrow(sig_tbl) == 0) {
  sig_tbl <- shape_log %>% filter(!is.na(P_value)) %>% slice_head(n = n_panels)
}

if (!allow_grouped_predictors) {
  analyte_tbl <- sig_tbl %>% filter(predictor_type %in% c("analyte", "total"))
  if (nrow(analyte_tbl) >= 1) sig_tbl <- analyte_tbl
}

forced_tbl <- shape_log %>%
  filter(PFAS %in% forced_predictors) %>%
  arrange(match(PFAS, forced_predictors))

remaining_tbl <- sig_tbl %>%
  filter(!PFAS %in% forced_tbl$PFAS) %>%
  arrange(FDR_P, P_value) %>%
  slice_head(n = max(0, n_panels - nrow(forced_tbl)))

selected <- bind_rows(forced_tbl, remaining_tbl) %>%
  distinct(PFAS, .keep_all = TRUE) %>%
  slice_head(n = n_panels)

if (nrow(selected) < n_panels) {
  fill_tbl <- shape_log %>%
    filter(!PFAS %in% selected$PFAS) %>%
    arrange(FDR_P, P_value) %>%
    slice_head(n = n_panels - nrow(selected))
  selected <- bind_rows(selected, fill_tbl) %>% distinct(PFAS, .keep_all = TRUE)
}

if (nrow(selected) == 0) stop("No predictors available for plotting.")

supp_selected <- shape_log %>%
  filter(!is.na(FDR_P), FDR_P < 0.05)

if (!allow_grouped_predictors) {
  supp_analyte <- supp_selected %>% filter(predictor_type %in% c("analyte", "total"))
  if (nrow(supp_analyte) >= 1) supp_selected <- supp_analyte
}

# ---------- Helpers ----------
fit_best_model <- function(df, best_model_name) {
  df <- df %>%
    mutate(
      log_x = log(x + 1e-6),
      x2 = x^2,
      x_adj = x + 1e-6
    )

  if (best_model_name == "linear") {
    fit <- lm(y ~ x, data = df)
  } else if (best_model_name == "log") {
    fit <- lm(y ~ log_x, data = df)
  } else if (best_model_name == "quadratic") {
    fit <- lm(y ~ x + x2, data = df)
  } else if (best_model_name == "exp") {
    fit <- nls(y ~ a * exp(b * x),
               start = list(a = mean(df$y), b = -0.01),
               data = df)
  } else if (best_model_name == "power") {
    fit <- nls(y ~ a * x_adj^b,
               start = list(a = mean(df$y), b = -0.5),
               data = df)
  } else {
    fit <- lm(y ~ x, data = df)
  }

  fit
}

get_display_name <- function(var_name) {
  if (var_name %in% names(vars_renamed)) return(unname(vars_renamed[[var_name]]))
  gsub("_", " ", var_name)
}

# ---------- Main panel data ----------
plot_dat_list <- list()
curve_dat_list <- list()

for (i in seq_len(nrow(selected))) {
  row <- selected[i, ]
  pfas_var <- row$PFAS

  if (!pfas_var %in% names(data_model)) next

  df_i <- data_model %>%
    transmute(
      x = .data[[pfas_var]],
      y = log(igy)
    ) %>%
    filter(!is.na(x), !is.na(y))

  if (nrow(df_i) < 8) next

  best_model_name <- row$Best_Model
  fit_i <- try(fit_best_model(df_i, best_model_name), silent = TRUE)
  if (inherits(fit_i, "try-error")) next

  new_x <- seq(min(df_i$x), max(df_i$x), length.out = 200)
  pred_df <- tibble(
    x = new_x,
    log_x = log(new_x + 1e-6),
    x2 = new_x^2,
    x_adj = new_x + 1e-6
  )
  pred_df$fit <- predict(fit_i, newdata = pred_df)

  panel_label <- paste0(
    get_display_name(pfas_var),
    "\n",
    best_model_name, " | FDR=", signif(row$FDR_P, 2),
    " | ", row$Direction
  )

  plot_dat_list[[pfas_var]] <- df_i %>% mutate(panel = panel_label)
  curve_dat_list[[pfas_var]] <- pred_df %>% transmute(x = x, fit = fit, panel = panel_label)
}

plot_dat <- bind_rows(plot_dat_list)
curve_dat <- bind_rows(curve_dat_list)

if (nrow(plot_dat) == 0 || nrow(curve_dat) == 0) {
  stop("No panel data were generated.")
}

# ---------- Main plot ----------
p <- ggplot(plot_dat, aes(x = x, y = y)) +
  geom_point(size = 1.4, alpha = 0.55, color = "gray35") +
  geom_line(data = curve_dat,
            aes(x = x, y = fit),
            linewidth = 1.0,
            color = "#1f78b4") +
  facet_wrap(~panel, scales = "free_x", ncol = 3) +
  labs(
    x = "PFAS predictor",
    y = "log(IgY)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.title = element_text(face = "bold"),
    panel.spacing = unit(0.9, "lines")
  )

ggsave(
  filename = file.path(out_dir, "Fig_1_key_takehomes_logIgY_panel.png"),
  plot = p, width = 10, height = 10, dpi = 600
)

ggsave(
  filename = file.path(out_dir, "Fig_1_key_takehomes_logIgY_panel.tiff"),
  plot = p, width = 10, height = 10, dpi = 600, compression = "lzw"
)

write.csv(
  selected,
  file.path(out_dir, "Fig_1_key_takehomes_selected_predictors.csv"),
  row.names = FALSE
)

# ---------- Supplemental plot: all significant ----------
if (nrow(supp_selected) > 0) {
  supp_plot_dat_list <- list()
  supp_curve_dat_list <- list()

  for (i in seq_len(nrow(supp_selected))) {
    row <- supp_selected[i, ]
    pfas_var <- row$PFAS

    if (!pfas_var %in% names(data_model)) next

    df_i <- data_model %>%
      transmute(
        x = .data[[pfas_var]],
        y = log(igy)
      ) %>%
      filter(!is.na(x), !is.na(y))

    if (nrow(df_i) < 8) next

    best_model_name <- row$Best_Model
    fit_i <- try(fit_best_model(df_i, best_model_name), silent = TRUE)
    if (inherits(fit_i, "try-error")) next

    new_x <- seq(min(df_i$x), max(df_i$x), length.out = 200)
    pred_df <- tibble(
      x = new_x,
      log_x = log(new_x + 1e-6),
      x2 = new_x^2,
      x_adj = new_x + 1e-6
    )
    pred_df$fit <- predict(fit_i, newdata = pred_df)

    panel_label <- paste0(
      get_display_name(pfas_var),
      "\n",
      best_model_name, " | FDR=", signif(row$FDR_P, 2),
      " | ", row$Direction
    )

    supp_plot_dat_list[[pfas_var]] <- df_i %>% mutate(panel = panel_label)
    supp_curve_dat_list[[pfas_var]] <- pred_df %>% transmute(x = x, fit = fit, panel = panel_label)
  }

  supp_plot_dat <- bind_rows(supp_plot_dat_list)
  supp_curve_dat <- bind_rows(supp_curve_dat_list)

  if (nrow(supp_plot_dat) > 0 && nrow(supp_curve_dat) > 0) {
    n_supp_panels <- length(unique(supp_plot_dat$panel))
    ncol_supp <- ifelse(n_supp_panels <= 9, 3, 4)
    supp_height <- ifelse(ncol_supp == 3,
                          3.2 * ceiling(n_supp_panels / 3),
                          2.8 * ceiling(n_supp_panels / 4))

    p_supp <- ggplot(supp_plot_dat, aes(x = x, y = y)) +
      geom_point(size = 1.1, alpha = 0.5, color = "gray35") +
      geom_line(data = supp_curve_dat,
                aes(x = x, y = fit),
                linewidth = 0.9,
                color = "#1f78b4") +
      facet_wrap(~panel, scales = "free_x", ncol = ncol_supp) +
      labs(
        x = "PFAS predictor",
        y = "log(IgY)"
      ) +
      theme_classic(base_size = 11) +
      theme(
        strip.text = element_text(face = "bold", size = 8.8),
        axis.title = element_text(face = "bold"),
        panel.spacing = unit(0.8, "lines")
      )

    ggsave(
      filename = file.path(out_dir, "Fig_S1_supp_all_significant_logIgY_panel.png"),
      plot = p_supp, width = 12, height = supp_height, dpi = 600
    )

    ggsave(
      filename = file.path(out_dir, "Fig_S1_supp_all_significant_logIgY_panel.tiff"),
      plot = p_supp, width = 12, height = supp_height, dpi = 600, compression = "lzw"
    )

    write.csv(
      supp_selected,
      file.path(out_dir, "Fig_S1_supp_all_significant_selected_predictors.csv"),
      row.names = FALSE
    )
  }
}

cat("Done: key take-home panel figure for Script 1\n")
