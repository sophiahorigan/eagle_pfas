##########################################
# Figure: Script 9 WBC mixed-model summaries
# - Heatmap of effect sizes
# - Forest-style dot plots with 95% CI
##########################################

rm(list = ls())

library(tidyverse)

input_file <- "./output/9_WBC_PFAS_Mixed_Models/PFAS_WBC_mixed_model_results.csv"
out_dir <- "./output/9_WBC_PFAS_Mixed_Models/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE)

if (nrow(df) == 0) stop("Script 9 output is empty.")

# Keep core plotting fields and remove rows without numeric estimates
plot_df <- df %>%
  filter(!is.na(estimate), !is.na(FDR_p)) %>%
  mutate(
    sig_fdr = FDR_p < 0.05,
    sig_nominal = p.value < 0.05,
    ci_low = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error
  )

if (nrow(plot_df) == 0) stop("No plottable rows (all estimates are NA).")

# Display-name convention for PFAS labels
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

plot_df <- plot_df %>%
  mutate(
    PFAS_display = to_display(PFAS),
    pfas_class = case_when(
      PFAS == "TOTAL_PFAS" ~ "total",
      grepl("^func_", PFAS) ~ "functional_group",
      TRUE ~ "other"
    ),
    WBC_display = recode(
      WBC_var,
      "Lperc" = "Lymphocyte",
      "Bperc" = "Basophil",
      "Eperc" = "Eosinophil",
      "Hperc" = "Heterophil",
      "Mperc" = "Monocyte",
      "Total_count" = "Total count"
    ),
    WBC_display = factor(
      WBC_display,
      levels = c("Lymphocyte", "Heterophil", "Eosinophil", "Monocyte", "Basophil", "Total count")
    )
  )

# Order PFAS so chain classes are together near the bottom, and Total PFAS is last
pfas_order_tbl <- plot_df %>%
  group_by(PFAS, PFAS_display, pfas_class) %>%
  summarize(mean_abs = mean(abs(estimate), na.rm = TRUE), .groups = "drop")

priority_bottom <- c("Total PFAS", "Legacy", "Next-Gen", "Short-chain", "Long-chain")

priority_tbl <- pfas_order_tbl %>%
  filter(PFAS_display %in% priority_bottom) %>%
  mutate(priority_rank = match(PFAS_display, priority_bottom)) %>%
  arrange(priority_rank)

other_bottom_tbl <- pfas_order_tbl %>%
  filter(pfas_class == "functional_group", !PFAS_display %in% priority_bottom) %>%
  arrange(mean_abs)

top_tbl <- pfas_order_tbl %>%
  filter(
    !PFAS_display %in% priority_bottom,
    pfas_class != "functional_group"
  ) %>%
  arrange(mean_abs)

# For y-axis factors in ggplot, early levels appear at the bottom.
pfas_order <- c(
  priority_tbl$PFAS_display,
  other_bottom_tbl$PFAS_display,
  top_tbl$PFAS_display
)

plot_df$PFAS_display <- factor(plot_df$PFAS_display, levels = pfas_order)

# ---------- 1) Heatmap ----------
heat <- ggplot(plot_df, aes(x = WBC_display, y = PFAS_display, fill = estimate)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_point(
    data = plot_df %>% filter(sig_fdr),
    shape = 21, size = 1.8, stroke = 0.35, fill = "black", color = "black"
  ) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    name = "Effect\n(beta)"
  ) +
  labs(
    x = "WBC endpoint",
    y = "PFAS predictor"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 7.7),
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.6)
  )

ggsave(
  file.path(out_dir, "Fig_WBC_heatmap_effects.png"),
  heat, width = 8.8, height = 10.5, dpi = 600
)
ggsave(
  file.path(out_dir, "Fig_WBC_heatmap_effects.tiff"),
  heat, width = 8.8, height = 10.5, dpi = 600, compression = "lzw"
)

# ---------- 2) Forest-style dot plot ----------
forest <- ggplot(plot_df, aes(x = estimate, y = PFAS_display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray55", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, color = "gray45", linewidth = 0.35) +
  geom_point(aes(color = sig_fdr), size = 1.5) +
  scale_color_manual(values = c("FALSE" = "#555555", "TRUE" = "#d73027"), guide = "none") +
  facet_wrap(~WBC_display, scales = "free_x", ncol = 3) +
  labs(
    x = "Model effect estimate (95% CI)",
    y = "PFAS predictor"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 6.9),
    strip.text = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "Fig_WBC_forest_effects.png"),
  forest, width = 12, height = 10, dpi = 600
)
ggsave(
  file.path(out_dir, "Fig_WBC_forest_effects.tiff"),
  forest, width = 12, height = 10, dpi = 600, compression = "lzw"
)

cat("Done: Script 9 visualization figures\n")
