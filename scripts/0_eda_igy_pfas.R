# ===============================================================
# Exploratory Data Analysis for Eagle IgY and PFAS dataset
# ===============================================================

# ---------- 0) Setup ----------
rm(list = ls())

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(GGally)
  library(viridis)
  library(tidyr)
  library(ggpubr)
  library(mgcv)
  library(cowplot)
})

# ---------- 1) Load data ----------
setwd("~/Desktop/EagleStats/")

data <- read.csv("data_clean_for_analysis.csv")

# Quick check
cat("Dataset has", nrow(data), "rows and", ncol(data), "columns\n")

# ---------- 2) Summary and structure ----------
summary(data)
str(data)

# ---------- 3) Basic distributions ----------
num_vars <- data %>%
  select(where(is.numeric)) %>%
  select(-lat, -long) # skip coordinates for now

# Histogram grid of numeric variables
pdf("plots/01_histograms.pdf", width = 10, height = 8)
num_vars %>%
  pivot_longer(everything()) %>%
  ggplot(aes(value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~name, scales = "free") +
  theme_bw() +
  labs(title = "Distributions of Numeric Variables")
dev.off()

# ---------- 4) Correlation heatmap ----------
pdf("plots/02_correlation_heatmap.pdf", width = 8, height = 6)
num_vars %>%
  ggcorr(label = TRUE, label_round = 2, low = "blue", high = "red", layout.exp = 2) +
  labs(title = "Correlation Matrix of Numeric Variables")
dev.off()

# ---------- 5) Key relationships ----------
# IgY vs TOTAL_PFAS
pdf("plots/03_IgY_vs_PFAS.pdf", width = 7, height = 5)
ggplot(data, aes(TOTAL_PFAS, igy)) +
  geom_point(aes(color = age_class), alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_color_viridis(discrete = TRUE) +
  theme_bw() +
  labs(title = "IgY vs Total PFAS", x = "Total PFAS (ng/g)", y = "IgY (units)")
dev.off()

# ---------- 6) Check linear vs nonlinear relationships ----------
pdf("plots/04_linear_vs_nonlinear.pdf", width = 7, height = 5)
p1 <- ggplot(data, aes(TOTAL_PFAS, igy)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_smooth(method = "gam", formula = y ~ s(x), color = "red", se = FALSE) +
  theme_bw() +
  labs(
    title = "Linear (blue) vs Nonlinear (red) fit: IgY ~ PFAS",
    subtitle = "Red line = GAM smoother; Blue line = linear model"
  )
print(p1)
dev.off()

# ---------- 7) IgY vs other predictors ----------
predictor_vars <- c("age", "body_condition", "PCA_score", "lat", "long")
pdf("plots/05_IgY_vs_predictors.pdf", width = 10, height = 8)
for (var in predictor_vars) {
  if (var %in% colnames(data)) {
    p <- ggplot(data, aes_string(x = var, y = "igy")) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", color = "black") +
      geom_smooth(method = "gam", formula = y ~ s(x), color = "red", se = FALSE) +
      theme_bw() +
      labs(
        title = paste("IgY vs", var),
        subtitle = "Linear (black) and Nonlinear (red) fits"
      )
    print(p)
  }
}
dev.off()

# ---------- 8) Pairwise scatterplots ----------
pdf("plots/06_pairwise_scatterplots.pdf", width = 10, height = 10)
vars_to_plot <- c("igy", "TOTAL_PFAS", "age", "body_condition", "PCA_score")
GGally::ggpairs(data[, vars_to_plot],
                upper = list(continuous = wrap("cor", size = 3)),
                lower = list(continuous = wrap("smooth_loess", alpha = 0.4))) +
  theme_bw()
dev.off()

# ---------- 9) Geographic visualization (optional) ----------
if(all(c("lat", "long") %in% colnames(data))) {
  pdf("plots/07_geographic_distribution.pdf", width = 7, height = 6)
  ggplot(data, aes(long, lat, color = TOTAL_PFAS)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_viridis(option = "C") +
    theme_bw() +
    labs(title = "Geographic Distribution of PFAS Levels",
         color = "Total PFAS")
  dev.off()
}

# ---------- 10) Model diagnostics ----------
# Example: IgY ~ PFAS + age + sex + body_condition
if(all(c("igy", "TOTAL_PFAS") %in% colnames(data))) {
  mod <- lm(igy ~ TOTAL_PFAS + age + body_condition, data = data)
  pdf("plots/08_model_diagnostics.pdf", width = 8, height = 8)
  par(mfrow = c(2, 2))
  plot(mod)
  dev.off()
}

cat("\nEDA complete — plots saved to /plots/\n")
