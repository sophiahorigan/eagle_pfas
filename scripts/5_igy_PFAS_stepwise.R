rm(list = ls())

setwd("~/Desktop/EagleStats/")

# PFAS -> IgY mixed-model stepwise pipeline
# Requirements: dplyr, tidyr, lme4, lmerTest, MuMIn, broom.mixed, performance, readr, ggplot2
pkgs <- c("dplyr","tidyr","lme4","lmerTest","MuMIn","broom.mixed","performance","readr","ggplot2")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(missing_pkgs)) install.packages(missing_pkgs)

library(dplyr); library(tidyr); library(lme4); library(lmerTest)
library(MuMIn); library(broom.mixed); library(performance); library(readr)
library(ggplot2)

df <- read.csv("./input/2023_2024_2025_AllData.csv")
key <- read.csv("./input/PFAS_subtypes.csv")

# ---------------- User options ----------------
zero_cutoff <- 0.50          # drop analytes with >50% zeros
small_offset <- 1e-6         # offset for log transform
log_response <- TRUE         # log-transform IgY response?
random_effect_col <- "nest_name"
save_dir <- "./output/pfas_igy_stepwise_results"
dir.create(save_dir, showWarnings = FALSE)

# Fixed-effect covariates (always candidates)
fixed_covs <- c("clysis","PCR_avmal","DNA_sex","lat","long","weight_kg","age")

# ---------- Identify PFAS analyte columns ----------
analyte_cols <- intersect(names(df), key$pfas_analyte)
if(length(analyte_cols) == 0) stop("No analyte columns found; check names(df) vs key$pfas_analyte.")

# ---------- Drop analytes with too many zeros ----------
zero_stats <- data.frame(analyte = analyte_cols,
                         n_non_na = NA_integer_,
                         prop_zero = NA_real_,
                         stringsAsFactors = FALSE)

for(i in seq_along(analyte_cols)) {
  x <- df[[analyte_cols[i]]]
  zero_stats$n_non_na[i] <- sum(!is.na(x))
  zero_stats$prop_zero[i] <- mean(!is.na(x) & x == 0, na.rm = TRUE)
}

drop_list <- zero_stats %>% filter(prop_zero > zero_cutoff) %>% pull(analyte)
keep_analytes <- setdiff(analyte_cols, drop_list)

write_csv(zero_stats, file.path(save_dir, "analyte_zero_stats.csv"))
cat("Dropped analytes (> ", zero_cutoff*100, "% zeros):\n")
if(length(drop_list)==0) cat("(none)\n") else cat(paste(drop_list, collapse = "\n"), "\n")
cat("\nKept analytes:\n", paste(keep_analytes, collapse = ", "), "\n\n")

# ---------- Compute TOTAL_PFAS if missing ----------
if(!"TOTAL_PFAS" %in% names(df)) {
  message("TOTAL_PFAS not found; computing as sum of kept analytes.")
  df <- df %>% mutate(TOTAL_PFAS = rowSums(select(., all_of(keep_analytes)), na.rm = TRUE))
}

# ---------- Compute PFAS group sums ----------
long_pf <- df %>%
  select(bird_id = any_of("bird_id"), all_of(keep_analytes)) %>%
  pivot_longer(cols = -bird_id, names_to = "analyte", values_to = "conc") %>%
  left_join(key %>% select(pfas_analyte, chain_cat, func_cat, pfas_age),
            by = c("analyte" = "pfas_analyte"))

chain_sums <- long_pf %>% filter(!is.na(chain_cat)) %>%
  group_by(bird_id, chain_cat) %>% summarize(chain_sum = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = chain_cat, values_from = chain_sum, values_fill = 0, names_prefix = "chain_")
func_sums  <- long_pf %>% filter(!is.na(func_cat)) %>%
  group_by(bird_id, func_cat) %>% summarize(func_sum  = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = func_cat,  values_from = func_sum,  values_fill = 0, names_prefix = "func_")
age_sums   <- long_pf %>% filter(!is.na(pfas_age)) %>%
  group_by(bird_id, pfas_age) %>% summarize(age_sum   = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pfas_age, values_from = age_sum, values_fill = 0, names_prefix = "age_")

df2 <- df %>%
  left_join(chain_sums, by = "bird_id") %>%
  left_join(func_sums,  by = "bird_id") %>%
  left_join(age_sums,   by = "bird_id")

group_cols <- c(grep("^chain_", names(df2), value = TRUE),
                grep("^func_",  names(df2), value = TRUE),
                grep("^age_",   names(df2), value = TRUE))

df2 <- df2 %>% mutate_at(vars(all_of(group_cols)), ~replace_na(., 0))

# ---------- Candidate predictors ----------
pfas_predictors <- keep_analytes
group_predictors <- group_cols
total_predictor  <- "TOTAL_PFAS"
candidate_predictors <- unique(c(pfas_predictors, group_predictors, total_predictor, fixed_covs))
cat("Candidate fixed-effect predictors (count =", length(candidate_predictors), "):\n")
cat(paste(candidate_predictors, collapse = ", "), "\n\n")

# ---------- Prepare modeling dataset ----------
if(!"igy" %in% names(df2)) stop("df must contain 'igy'")
mod_df <- df2 %>% filter(!is.na(igy)) %>%
  filter_at(vars(all_of(fixed_covs)), all_vars(!is.na(.)))

cat("Rows after cleaning:", nrow(mod_df), "\n")
if(nrow(mod_df) < 20) warning("Fewer than 20 observations; model may not be reliable.")

# ---------- Clean dataset and remove empty / constant columns ----------
mod_df_clean <- mod_df %>%
  mutate(across(where(is.factor), as.character)) %>%
  mutate(across(where(is.character), ~na_if(.x, ""))) %>%
  mutate(across(where(is.character), ~na_if(.x, "#N/A"))) %>%
  select(where(~length(unique(na.omit(.x))) > 1)) %>%
  mutate(across(where(is.character), as.factor))

# ---------- Build global formula ----------
rhs <- paste(candidate_predictors, collapse = " + ")
response_expr <- if(log_response) paste0("log(igy + ", small_offset, ")") else "igy"

rhs_vars <- unlist(strsplit(rhs, " \\+ "))
existing_vars <- rhs_vars[rhs_vars %in% colnames(mod_df_clean)]
missing_vars <- setdiff(rhs_vars, existing_vars)
if(length(missing_vars) > 0) warning("Predictors removed (not in dataset): ", paste(missing_vars, collapse = ", "))

global_formula <- as.formula(
  paste0(response_expr, " ~ ", paste(existing_vars, collapse = " + "), " + (1 | ", random_effect_col, ")")
)
cat("Global model formula:\n"); print(global_formula)

# ---------- Fit global mixed model ----------
global_mod <- tryCatch({
  lmerTest::lmer(global_formula, data = mod_df_clean, REML = FALSE, control = lmerControl(check.nobs.vs.nlev = "ignore"))
}, error = function(e){
  stop("Global lmer failed: ", e$message)
})
cat("Global model fitted.\n")

# ---------- Stepwise selection ----------
step_mod <- tryCatch({
  lmerTest::step(global_mod, reduce.fixed = TRUE, reduce.random = FALSE, keep = NULL)
}, error = function(e){
  message("lmerTest::step failed: ", e$message, "\nFalling back to drop1-AIC backward selection.")
  current_mod <- global_mod
  improved <- TRUE
  while(improved){
    d1 <- drop1(current_mod, test = "Chisq")
    d1 <- d1[!rownames(d1) %in% "(Intercept)", , drop = FALSE]
    if(nrow(d1) == 0 || all(is.na(d1$AIC))){ improved <- FALSE; break }
    
    aic_vals <- d1$AIC
    term_names <- rownames(d1)
    aic_df <- data.frame(term = term_names, AIC = aic_vals, stringsAsFactors = FALSE)
    
    best_row <- which.min(aic_df$AIC)
    aic_current <- AIC(current_mod)
    aic_candidate <- aic_df$AIC[best_row]
    
    if(!is.na(aic_candidate) && (aic_candidate + 0) < aic_current - 1.999){
      term_to_drop <- aic_df$term[best_row]
      cat("Dropping term", term_to_drop, "because AIC decreases from", aic_current, "to", aic_candidate, "\n")
      current_mod <- update(current_mod, paste(". ~ . -", term_to_drop))
    } else improved <- FALSE
  }
  current_mod
})

cat("Stepwise selection finished.\n")
print(summary(step_mod))

# Check which fixed-effect predictors were actually retained
fixed_retained_df <- broom.mixed::tidy(step_mod, effects = "fixed", conf.int = TRUE)

cat("\n=== FIXED EFFECTS RETAINED IN STEPWISE MODEL ===\n")
if(nrow(fixed_retained_df) == 0){
  cat("No fixed effects retained (only intercept remains).\n")
} else {
  print(fixed_retained_df[, c("term", "estimate", "conf.low", "conf.high", "p.value")])
}

# Optional: list predictors that were dropped
all_candidate <- existing_vars  # variables used in global_formula
retained_terms <- fixed_retained_df$term
dropped_terms <- setdiff(all_candidate, retained_terms)
cat("\nPredictors dropped during stepwise selection:\n")
if(length(dropped_terms) == 0) cat("(none)\n") else cat(paste(dropped_terms, collapse = ", "), "\n")


# ---------- Save results ----------
saveRDS(global_mod, file = file.path(save_dir, "global_mixed_model.rds"))
saveRDS(step_mod,   file = file.path(save_dir, "stepwise_selected_mixed_model.rds"))

fixed_retained_df <- broom.mixed::tidy(step_mod, effects = "fixed", conf.int = TRUE)
write_csv(fixed_retained_df, file.path(save_dir, "selected_fixed_effects.csv"))
cat("Saved selected fixed effects to:", file.path(save_dir, "selected_fixed_effects.csv"), "\n")

# Diagnostics plots
resid_df <- data.frame(fitted = predict(step_mod), resid = residuals(step_mod))
p1 <- ggplot(resid_df, aes(x = fitted, y = resid)) + geom_point() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted", x = "Fitted values", y = "Residuals")
p1
ggsave(file.path(save_dir, "resid_vs_fitted.png"), plot = p1, width = 7, height = 5, dpi = 300)

p2 <- ggplot(resid_df, aes(sample = resid)) + stat_qq() + stat_qq_line() +
  labs(title = "QQ-plot of residuals")
p2
ggsave(file.path(save_dir, "qqplot_resid.png"), plot = p2, width = 6, height = 6, dpi = 300)

sink(file.path(save_dir, "stepwise_model_summary.txt"))
print(summary(step_mod))
sink()

cat("\nAll outputs saved to:", normalizePath(save_dir), "\n")
cat("Dropped analytes saved to:", file.path(save_dir, "analyte_zero_stats.csv"), "\n\n")

cat("==== QUICK REPORT ====\n")
cat("Dropped analytes (> ", zero_cutoff*100, "% zeros):\n")
if(length(drop_list)==0) cat("(none)\n") else cat(paste(drop_list, collapse = "\n"), "\n")
cat("\nSelected fixed effects:\n")
print(fixed_retained_df)
