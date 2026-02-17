# ============================================================
# PFAS–IgY Analysis Script (Simplified, No LOD Handling)
# ============================================================

rm(list = ls())

# --- 1. Load packages ---
suppressPackageStartupMessages({
  library(tidyverse)
  library(car)
  library(ggpubr)
  library(ggcorrplot)
  library(lme4)       # For mixed-effects models
  library(lmerTest)   # Optional: gives p-values
})

# --- 2. Load data ---
setwd("~/Desktop/EagleStats/")
df <- read.csv("2023_2024_2025_AllData.csv")

# --- 3. Select relevant variables ---
pfas_cols <- grep("^PF", colnames(df), value = TRUE)

pfas_cols <- union(pfas_cols, "TOTAL_PFAS")

df_sub <- df %>%
  select(
    bird_id, igy, field_sex, age, year, lat, long, nest_name,
    all_of(pfas_cols)
  ) %>%
  drop_na(igy) %>%
  mutate(
    igy_log = log(igy + 1)
  )

# --- 4. Drop PFAS with zero or near-zero variance ---
pfas_good1 <- pfas_cols[
  sapply(pfas_cols, function(v) {
    x <- df_sub[[v]]
    xf <- x[is.finite(x)]
    length(unique(xf)) >= 3 && sd(xf) > 0
  })
]

dropped1 <- setdiff(pfas_cols, pfas_good1)
if (length(dropped1) > 0) {
  message("Dropping PFAS with too little variation: ",
          paste(dropped1, collapse = ", "))
}

# Keep only usable PFAS
df_sub <- df_sub %>% select(
  bird_id, igy, igy_log, field_sex, age, year, lat, long, nest_name,
  all_of(pfas_good1)
)

df_sub <- na.omit(df_sub)

# --- 5. Log-transform PFAS ---
df_sub <- df_sub %>%
  mutate(across(all_of(pfas_good1), ~log10(.x + 1), .names = "log10_{col}"))

pfas_log_cols <- grep("^log10_", colnames(df_sub), value = TRUE)

# --- 6. Shapiro-Wilk normality tests ---
normality <- df_sub %>%
  summarise(across(
    c(igy_log, all_of(pfas_log_cols)),
    ~{
      x <- .x[is.finite(.x)]
      if (length(unique(x)) >= 3) {
        tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      }
    },
    .names = "{.col}"
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "shapiro_p")

print("Shapiro-Wilk p-values:")
print(normality)

# Remove PFAS whose Shapiro test cannot run (NA)
pfas_log_cols2 <- normality %>%
  filter(!is.na(shapiro_p)) %>%
  pull(variable) %>%
  grep("^log10_", ., value = TRUE)

# --- 7. PFAS correlation matrix ---
cor_matrix <- cor(
  df_sub %>% select(all_of(pfas_log_cols2)),
  use = "pairwise.complete.obs",
  method = "spearman"
)

ggcorrplot(
  cor_matrix, hc.order = TRUE, type = "lower",
  lab = FALSE, title = "PFAS Correlation Matrix"
)

# --- 8. Correlation between IgY and PFAS ---
corrs <- map_dfr(
  pfas_log_cols2,
  function(var) {
    x <- df_sub[[var]]
    shp <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA)
    
    method <- if (!is.na(shp) && shp > 0.05) "pearson" else "spearman"
    
    ct <- suppressWarnings(cor.test(df_sub$igy_log, x, method = method))
    
    tibble(
      PFAS = var,
      method = method,
      cor = ct$estimate,
      p = ct$p.value
    )
  }
)

print(corrs)

# --- 9. Linear modeling ---
df_sub$TOTAL_PFAS_log10 <- log10(df_sub$TOTAL_PFAS + 1)

lm_mod <- lm(igy_log ~ TOTAL_PFAS_log10 + field_sex + age + year, data = df_sub)
summary(lm_mod)
car::Anova(lm_mod, type = "III")

# --- 10. Diagnostics ---
par(mfrow = c(2, 2))
plot(lm_mod)

# --- 11. Visualization ---
ggplot(df_sub, aes(x = TOTAL_PFAS_log10, y = igy_log, color = field_sex)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  labs(
    x = "log10(Total PFAS + 1)",
    y = "log(IgY + 1)",
    title = "Relationship Between IgY and PFAS Exposure"
  ) +
  theme_bw()

# --- 12. Save results ---
write_csv(corrs, "pfas_igy_correlations.csv")
saveRDS(lm_mod, "pfas_igy_lm_model.rds")


## --- MIXED EFFECT, NEST -----
# Ensure nest_name is a factor
df_sub$nest_name <- as.factor(df_sub$nest_name)

# Fit mixed-effects model
lmm_mod <- lmer(
  igy_log ~ TOTAL_PFAS_log10 + field_sex + age + year + (1 | nest_name),
  data = df_sub
)

# Model summary (includes fixed effects estimates)
summary(lmm_mod)

# Type III ANOVA for fixed effects
car::Anova(lmm_mod, type = "III")



## --- NEST -----
lmm_mod <- lmer(
  igy_log ~ TOTAL_PFAS_log10 + nest_name + (1 | field_sex),
  data = df_sub
)
summary(lmm_mod)
car::Anova(lmm_mod, type = "III")

