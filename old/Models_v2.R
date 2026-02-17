
################################################################################
# Eagle IgY ~ PFAS analysis pipeline
# - Sections:
#   1) libraries & helpers
#   2) data prep & nondetect handling
#   3) EDA & plots
#   4) derived variables (body condition, exposure metrics)
#   5) primary models: LMM + GAM
#   6) mixture models: WQS, BKMR
#   7) sensitivity analyses + diagnostics
################################################################################
rm(list = ls())

# ---------- 1) Libraries & helpers ----------
# Install missing packages if needed:
# pkg <- c("tidyverse","lubridate","lme4","lmerTest","mgcv","gWQS","bkmr","corrplot",
#          "car","cowplot","patchwork","ggpubr","broom.mixed")
# inst <- pkg[!(pkg %in% installed.packages()[,"Package"])]
# if(length(inst)) install.packages(inst, repos="https://cloud.r-project.org")

# load
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(mgcv)
library(gWQS)
library(bkmr)
library(corrplot)
library(car)
library(cowplot)
library(patchwork)
library(ggpubr)
library(broom.mixed)

# helper: safe log10 transform (adds tiny constant)
safe_log10 <- function(x, add = 1e-3) {
  # replace -Inf/NaN with NA after transform
  out <- log10(x + add)
  out[is.infinite(out) | is.nan(out)] <- NA
  out
}

# helper: replace zeros (possible nondetects) with half-minimum positive/sqrt(2)
replace_zeros_with_halfmin <- function(x) {
  if(!is.numeric(x)) return(x)
  pos <- x[x > 0 & !is.na(x)]
  if(length(pos) == 0) return(x) # nothing to replace against
  halfmin <- min(pos, na.rm=TRUE) / sqrt(2)
  x[x == 0 & !is.na(x)] <- halfmin
  return(x)
}

# ---------- 2) Load data ----------
setwd("~/Desktop/EagleStats/")

data <- read.csv("2023_2024_2025_AllData.csv")

# assuming object 'data' exists in environment
if(!exists("data")) stop("Please load your dataframe into the R object named `data` before running this script.")

# Quick check of required columns (based on your list)
required_cols <- c("bird_id","igy","age","DNA_sex","field_sex","weight_kg","culman_mm",
                   "nest_name","year","date","lat","long","TOTAL_PFAS")

missing_cols <- setdiff(required_cols, colnames(data))
if(length(missing_cols) > 0) {
  warning("The following 'required' cols are missing from your data: ", paste(missing_cols, collapse=", "))
  # continue — some models will be skipped if fields missing
}

# ---------- 3) Identify PFAS columns ----------
# Your PFAS columns appear to start with common prefixes: we'll pick columns that match PF or specific names
#pfas_candidates <- colnames(data)[grepl("^PF|^X[0-9]|^TOTAL_PFAS|PFOS|PFOA|HFPO.DA|ADONA|FOSA", colnames(data), ignore.case = FALSE)]
# more conservative: any column with "PF" substring or TOTAL_PFAS
pfas_candidates <- unique(colnames(data)[grepl("PF|TOTAL_PFAS|ADONA|FOSA|HFPO", colnames(data))])

# Filter to numeric columns only
pfas_cols <- pfas_candidates[sapply(data[pfas_candidates], is.numeric)]
pfas_cols <- sort(pfas_cols)
message("Detected PFAS columns: ", paste(head(pfas_cols,30), collapse=", "))
# ensure TOTAL_PFAS exists
if(!"TOTAL_PFAS" %in% colnames(data)) {
  warning("TOTAL_PFAS column not detected; script will compute a sum-of-PFAS automatically.")
}

# ---------- 4) Non-detect handling ----------
# Common field encoding: non-detects may be 0, NA, or flagged elsewhere.
# Strategy implemented here:
#  - Replace numeric zeros with half-minimum-positive / sqrt(2) per-analyte
#  - If a column is mostly zeros/NA (>50%), still replace zeros the same way but flag for sensitivity
data_clean <- data %>% mutate(across(all_of(pfas_cols), replace_zeros_with_halfmin))

# Flag analytes with many nondetects (>=50% zeros or NAs)
nd_flags <- sapply(data[pfas_cols], function(x) mean(is.na(x) | x == 0))
nd_flags_df <- tibble(analyte = names(nd_flags), prop_nd = as.numeric(nd_flags))
message("Analytes with high proportion non-detects (>=0.5):\n",
        paste(nd_flags_df %>% filter(prop_nd >= 0.5) %>% pull(analyte), collapse=", "))

# If TOTAL_PFAS missing, compute as sum of numeric PFAS columns
if(!"TOTAL_PFAS" %in% colnames(data_clean)) {
  data_clean <- data_clean %>%
    mutate(TOTAL_PFAS = rowSums(select(., all_of(pfas_cols)), na.rm=TRUE))
}

# ---------- 5) Derived variables ----------
# 5.1 Body condition index:
# Use log-weight ~ log(culmen_mm) (or other size metric) and residuals as condition
if(all(c("weight_kg","culman_mm") %in% colnames(data_clean))) {
  
  # Identify complete cases for the model
  complete_rows <- complete.cases(data_clean[, c("weight_kg", "culman_mm")])
  
  # Fit model only on complete rows
  bc_mod <- lm(log(weight_kg + 1e-6) ~ log(culman_mm + 1e-6),
               data = data_clean[complete_rows, ])
  
  # Initialize new column as NA
  data_clean$body_condition <- NA_real_
  
  # Fill residuals only for rows used in the model
  data_clean$body_condition[complete_rows] <- resid(bc_mod)
  
} else {
  data_clean <- data_clean %>% mutate(body_condition = NA_real_)
  warning("weight_kg or culman_mm missing -> body_condition set to NA")
}


# 5.2 Age / Sex harmonization
# Use DNA_sex if present, otherwise field_sex
if("DNA_sex" %in% colnames(data_clean)) {
  data_clean <- data_clean %>% mutate(sex = if_else(!is.na(DNA_sex) & DNA_sex != "", DNA_sex, field_sex))
} else {
  data_clean <- data_clean %>% mutate(sex = field_sex)
}

# 5.3 Log transforms for modeling
# We'll produce log10 of IgY and PFAS exposures
data_clean <- data_clean %>%
  mutate(log_igy = safe_log10(igy),
         log_TOTAL_PFAS = safe_log10(TOTAL_PFAS),
         log_TOTAL_PFAS_z = scale(log_TOTAL_PFAS))

# Also create log-transformed individual PFAS (and scaled) for mixture methods
for(col in pfas_cols) {
  newname <- paste0("log_", col)
  data_clean[[newname]] <- safe_log10(data_clean[[col]])
}
# scaled versions
log_pfas_cols <- paste0("log_", pfas_cols)
scaled_log_pfas_cols <- paste0(log_pfas_cols, "_z")
data_clean <- data_clean %>% mutate(across(all_of(log_pfas_cols), ~scale(.x), .names = "{.col}_z"))

# 5.4 PCA on PFAS (if many analytes)
# Drop analytes with near-zero variance
pfas_for_pca <- log_pfas_cols[sapply(data_clean[log_pfas_cols], function(x) sd(x, na.rm=TRUE) > 1e-6)]
if(length(pfas_for_pca) >= 2) {
  pca_mat <- data_clean %>% select(all_of(pfas_for_pca)) %>% drop_na()
  pca_obj <- prcomp(pca_mat, center=TRUE, scale.=TRUE)
  # project full dataset with predict (handles NA rows)
  pc_scores <- predict(pca_obj, newdata = data_clean %>% select(all_of(pfas_for_pca)))
  data_clean <- bind_cols(data_clean, tibble(PCA_score = pc_scores[,1])) # first PC
} else {
  data_clean$PCA_score <- NA_real_
}

# ---------- 6) Exploratory data analysis ----------
# 6.1 Distribution of IgY and TOTAL_PFAS
p1 <- ggplot(data_clean, aes(x=igy)) + geom_histogram(bins=30) + ggtitle("IgY distribution (raw)")
p2 <- ggplot(data_clean, aes(x=log_igy)) + geom_histogram(bins=30) + ggtitle("IgY distribution (log10)")
p3 <- ggplot(data_clean, aes(x=TOTAL_PFAS)) + geom_histogram(bins=30) + ggtitle("TOTAL_PFAS distribution")
p4 <- ggplot(data_clean, aes(x=log_TOTAL_PFAS)) + geom_histogram(bins=30) + ggtitle("log10 TOTAL_PFAS")
ggsave("igs_pfas_distributions.png", plot = (p1 + p2) / (p3 + p4), width = 12, height = 8)

# 6.2 Correlation among top PFAS
if(length(pfas_for_pca) >= 2) {
  cor_mat <- cor(data_clean %>% select(all_of(pfas_for_pca)), use="pairwise.complete.obs")
  png("pfas_correlation.png", width=1000, height=800)
  corrplot(cor_mat, method="color", type="upper", tl.cex=0.7, mar=c(1,1,1,1))
  #dev.off()
}

# 6.3 Scatter IgY vs TOTAL_PFAS with smoother
p5 <- ggplot(data_clean, aes(x=log_TOTAL_PFAS, y=log_igy)) +
  geom_point(alpha=0.6) + geom_smooth(method="loess") +
  ggtitle("log IgY vs log TOTAL_PFAS")
print(p5)
ggsave("igy_vs_totalpfas.png", p5, width=6, height=5)

# ---------- 7) Primary model: Linear Mixed Model (LMM) ----------
# Outcome: log_igy
# Primary exposure: log_TOTAL_PFAS_z (z-scored)
# Covariates: age_class, sex, body_condition, year, (1|nest_name)
# If nest_name misssing, fallback to field_id
rand_effect <- if("nest_name" %in% colnames(data_clean)) "nest_name" else if("field_id" %in% colnames(data_clean)) "field_id" else NULL
fix_covars <- c()
#if("age_class" %in% colnames(data_clean)) fix_covars <- c(fix_covars, "age_class")
if("sex" %in% colnames(data_clean)) fix_covars <- c(fix_covars, "sex")
#if("body_condition" %in% colnames(data_clean)) fix_covars <- c(fix_covars, "body_condition")
if("year" %in% colnames(data_clean)) fix_covars <- c(fix_covars, "year")
# construct formula
rhs <- paste(c("log_TOTAL_PFAS_z", fix_covars), collapse = " + ")
if(!is.null(rand_effect)) {
  form_lmm <- as.formula(paste0("log_igy ~ ", rhs, " + (1 | ", rand_effect, ")"))
} else {
  form_lmm <- as.formula(paste0("log_igy ~ ", rhs))
}

# Fit LMM
cat("Fitting LMM with formula:\n"); print(form_lmm)
lmm_fit <- tryCatch({
  lmer(form_lmm, data = data_clean, REML = FALSE)
}, error = function(e){
  warning("LMM failed: ", e$message)
  NULL
})

if(!is.null(lmm_fit)) {
  summary(lmm_fit)
  # extract tidy summary
  broom::tidy(lmm_fit, effects="fixed") %>% write_csv("lmm_fixed_effects.csv")
  saveRDS(lmm_fit, "lmm_fit.rds")
  # diagnostic plots
  png("lmm_resid_vs_fitted.png", width=800, height=600)
  plot(lmm_fit, which = 1)
  #dev.off()
  png("lmm_qq.png", width=800, height=600)
  qqnorm(resid(lmm_fit)); qqline(resid(lmm_fit))
  #dev.off()
}

# ---------- 8) Check multicollinearity if multiple congeners used ----------
# Example: simple LM with several top congeners (if present)
top_congeners <- c("PFOS","PFOA","PFNA","PFDA")[c("PFOS","PFOA","PFNA","PFDA") %in% pfas_cols]
if(length(top_congeners) >= 2) {
  lm_formula <- as.formula(paste0("log_igy ~ ", paste(paste0("log_", top_congeners), collapse = " + "), " + ", paste(fix_covars, collapse=" + ")))
  lm1 <- lm(lm_formula, data = data_clean)
  cat("VIFs for model with multiple PFAS congeners:\n")
  print(car::vif(lm1))
}

# ---------- 9) GAM to explore non-linearity ----------
# Add spline for log_TOTAL_PFAS
gam_formula <- as.formula(paste0("log_igy ~ s(log_TOTAL_PFAS) + ", paste(fix_covars, collapse = " + ")))
cat("Fitting GAM:\n"); print(gam_formula)
gam_fit <- tryCatch({
  gam(gam_formula, data = data_clean, method="REML")
}, error=function(e){ warning("GAM failed: ", e$message); NULL })

if(!is.null(gam_fit)) {
  png("gam_smooth_plot.png", width=800, height=600)
  plot(gam_fit, pages=1)
  dev.off()
  saveRDS(gam_fit, "gam_fit.rds")
  summary(gam_fit) %>% capture.output(file="gam_summary.txt")
}

# ---------- 11) Sensitivity analyses ----------
# 11.1 Alternate exposure: PCA_score
if(!all(is.na(data_clean$PCA_score))) {
  sens_form <- as.formula(paste0("log_igy ~ scale(PCA_score) + ", paste(fix_covars, collapse = " + "), 
                                 if(!is.null(rand_effect)) paste0(" + (1 | ", rand_effect, ")") else ""))
  try({
    sens_lmm <- lmer(sens_form, data=data_clean, REML = FALSE)
    saveRDS(sens_lmm, "sensitivity_pca_lmm.rds")
    broom::tidy(sens_lmm, effects="fixed") %>% write_csv("sensitivity_pca_fixed_effects.csv")
  }, silent = TRUE)
}

# 11.2 Exclude high non-detect analytes and recompute TOTAL_PFAS
high_nd_analytes <- nd_flags_df %>% filter(prop_nd >= 0.5) %>% pull(analyte)
if(length(high_nd_analytes) > 0) {
  df_no_highnd <- data_clean %>% select(-all_of(high_nd_analytes))
  # recompute TOTAL_PFAS_nohighnd
  pfas_small <- setdiff(pfas_cols, high_nd_analytes)
  if(length(pfas_small) >= 1) {
    df_no_highnd <- df_no_highnd %>% mutate(TOTAL_PFAS_no_high_nd = rowSums(select(., all_of(pfas_small)), na.rm=TRUE),
                                            log_TOTAL_PFAS_no_high_nd = safe_log10(TOTAL_PFAS_no_high_nd),
                                            log_TOTAL_PFAS_no_high_nd_z = scale(log_TOTAL_PFAS_no_high_nd))
    # rerun simple LMM
    try({
      form2 <- as.formula(paste0("log_igy ~ log_TOTAL_PFAS_no_high_nd_z + ", paste(fix_covars, collapse=" + "),
                                 if(!is.null(rand_effect)) paste0(" + (1 | ", rand_effect, ")") else ""))
      lmm2 <- lmer(form2, data = df_no_highnd, REML = FALSE)
      saveRDS(lmm2, "lmm_no_high_nd.rds")
      broom::tidy(lmm2, effects="fixed") %>% write_csv("lmm_no_high_nd_fixed_effects.csv")
    }, silent=TRUE)
  }
}

# 11.3 Handling lipid adjustment: if lipid measure exists, compare approaches
if("L" %in% colnames(data_clean)) {
  # Option A: include lipid as covariate
  form_lipid <- as.formula(paste0("log_igy ~ log_TOTAL_PFAS_z + L + ", paste(setdiff(fix_covars, "body_condition"), collapse=" + "),
                                  if(!is.null(rand_effect)) paste0(" + (1 | ", rand_effect, ")") else ""))
  try({
    lmm_lipid <- lmer(form_lipid, data = data_clean, REML = FALSE)
    saveRDS(lmm_lipid, "lmm_lipid.rds")
    broom::tidy(lmm_lipid, effects="fixed") %>% write_csv("lmm_lipid_fixed_effects.csv")
  }, silent=TRUE)
}

# ---------- 12) Effect modification tests ----------
# Example: PFAS * age_class, PFAS * sex
if("age_class" %in% colnames(data_clean)) {
  int_form <- as.formula(paste0("log_igy ~ log_TOTAL_PFAS_z * age_class + ", paste(setdiff(fix_covars, "age_class"), collapse=" + "),
                                if(!is.null(rand_effect)) paste0(" + (1 | ", rand_effect, ")") else ""))
  try({
    lmm_int <- lmer(int_form, data = data_clean, REML = FALSE)
    saveRDS(lmm_int, "lmm_pfasp_age_interaction.rds")
    broom::tidy(lmm_int, effects="fixed") %>% write_csv("lmm_pfasp_age_interaction_fixed_effects.csv")
  }, silent = TRUE)
}
if("sex" %in% colnames(data_clean)) {
  int_form2 <- as.formula(paste0("log_igy ~ log_TOTAL_PFAS_z * sex + ", paste(setdiff(fix_covars, "sex"), collapse=" + "),
                                 if(!is.null(rand_effect)) paste0(" + (1 | ", rand_effect, ")") else ""))
  try({
    lmm_int2 <- lmer(int_form2, data = data_clean, REML = FALSE)
    saveRDS(lmm_int2, "lmm_pfasp_sex_interaction.rds")
    broom::tidy(lmm_int2, effects="fixed") %>% write_csv("lmm_pfasp_sex_interaction_fixed_effects.csv")
  }, silent = TRUE)
}

# ---------- 13) Model summary output & plots ----------
# Produce marginal predicted plot for LMM (if fitted)
if(!is.null(lmm_fit)) {
  # create prediction grid across observed PFAS
  newdf <- data_clean %>% 
    summarize_at(vars(log_TOTAL_PFAS), list(min=min, max=max), na.rm=TRUE) %>%
    pivot_longer(everything(), names_to="var", values_to="value")
  # better: create seq across range
  rng <- range(data_clean$log_TOTAL_PFAS, na.rm=TRUE)
  pred_grid <- tibble(log_TOTAL_PFAS = seq(rng[1], rng[2], length.out = 50))
  # set covariates to modal or mean
  for(cv in fix_covars) {
    if(is.numeric(data_clean[[cv]])) {
      pred_grid[[cv]] <- mean(data_clean[[cv]], na.rm=TRUE)
    } else {
      pred_grid[[cv]] <- names(which.max(table(data_clean[[cv]])))
    }
  }
  # if random effect present, set to first level
  if(!is.null(rand_effect)) pred_grid[[rand_effect]] <- unique(na.omit(data_clean[[rand_effect]]))[1]
  # add z-scored exposure
  pred_grid$log_TOTAL_PFAS_z <- scale(data_clean$log_TOTAL_PFAS)[1] * 0 + (pred_grid$log_TOTAL_PFAS - mean(data_clean$log_TOTAL_PFAS, na.rm=TRUE)) / sd(data_clean$log_TOTAL_PFAS, na.rm=TRUE)
  # predict using lmer (predict will include random effect by default; include re.form=NA to get population-level)
  pred_grid$pred_log_igy <- predict(lmm_fit, newdata = pred_grid, re.form = NA)
  p_pred <- ggplot(pred_grid, aes(x=log_TOTAL_PFAS, y=pred_log_igy)) + geom_line() +
    labs(x="log10 TOTAL_PFAS", y="Predicted log10 IgY") + ggtitle("Predicted effect of TOTAL_PFAS on log IgY (LMM)")
  ggsave("predicted_lmm_totalpfas.png", p_pred, width=7, height=5)
  write_csv(pred_grid, "lmm_prediction_grid.csv")
}

# ---------- 14) Save cleaned data ----------
write_csv(data_clean, "data_clean_for_analysis.csv")

# ---------- 15) Quick report summary (console) ----------
cat("\n==== Analysis complete ====\n")
if(!is.null(lmm_fit)) {
  cat("LMM fixed effects:\n")
  print(broom::tidy(lmm_fit, effects="fixed"))
}
if(exists("wqs_mod")) cat("WQS model saved to wqs_model.rds and weights in wqs_weights.csv\n")
if(exists("fit_bkmr")) cat("BKMR demo run saved to bkmr_fit_demo.rds (short run - increase iterations for final analysis)\n")
cat("Cleaned data saved as data_clean_for_analysis.csv\n")
cat("Plots saved: igs_pfas_distributions.png, pfas_correlation.png (if PFAS >=2), igy_vs_totalpfas.png, predicted_lmm_totalpfas.png (if LMM fitted)\n")

################################################################################
# End of script
# Notes:
# - Adjust number of bootstrap iterations (b) in gwqs and iterations in bkmr for production runs.
# - Consider using brms for fully Bayesian hierarchical modeling (useful for small samples or complex censoring).
# - Document how non-detects were handled: zeros replaced with half-minimum-positive/sqrt(2) by analyte.
# - If you have LOD values available for each analyte, replace the zero-replacement strategy with LOD/√2 or use censored modeling approaches.
################################################################################
