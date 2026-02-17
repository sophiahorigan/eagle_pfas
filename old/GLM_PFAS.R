
rm(list = ls())

setwd("~/Desktop/EagleStats/")

# PFAS -> IgY mixed-model stepwise pipeline
# Requirements: dplyr, tidyr, lme4, lmerTest, glmmTMB (optional), MuMIn, broom.mixed, performance, readr
pkgs <- c("dplyr","tidyr","lme4","lmerTest","MuMIn","broom.mixed","performance","readr","ggplot2")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(missing_pkgs)) install.packages(missing_pkgs)

library(dplyr); library(tidyr); library(lme4); library(lmerTest)
library(MuMIn); library(broom.mixed); library(performance); library(readr)
library(ggplot2)

df <- read.csv("./input/2023_2024_2025_AllData.csv")
key <- read.csv("./input/PFAS_subtypes.csv")

# ---------------- User options ----------------
zero_cutoff <- 0.50          # drop analytes with >50% zeros (adjustable)
small_offset <- 1e-6         # offset used for logs to avoid -Inf
log_response <- TRUE         # TRUE: model log(igy + offset); FALSE: model igy directly
random_effect_col <- "nest_name"  # your nest grouping column (you used nest_name in df)
save_dir <- "pfas_igy_stepwise_results"
dir.create(save_dir, showWarnings = FALSE)

# Fixed-effect covariates (will always be included as candidates for stepwise)
fixed_covs <- c("clysis","PCR_avmal","DNA_sex","lat","long","weight_kg","age")

# ---------- Detect PFAS analyte columns ----------
# Assumes df and key are already loaded in environment (as in your previous message).
if(!exists("df") || !exists("key")) stop("Please load `df` and `key` into the environment before running this script.")

# Use key$pfas_analyte to identify analyte columns in df
analyte_cols <- intersect(names(df), key$pfas_analyte)
if(length(analyte_cols) == 0) stop("No analyte columns found; check names(df) vs key$pfas_analyte.")

# ---------- Compute zero prevalence and drop analytes ----------
zero_stats <- data.frame(analyte = analyte_cols,
                         n_non_na = NA_integer_,
                         prop_zero = NA_real_,
                         stringsAsFactors = FALSE)

for(i in seq_along(analyte_cols)) {
  col <- analyte_cols[i]
  x <- df[[col]]
  zero_stats$n_non_na[i] <- sum(!is.na(x))
  zero_stats$prop_zero[i] <- mean(!is.na(x) & x == 0, na.rm = TRUE)
}

# analytes to drop and keep
drop_list <- zero_stats %>% filter(prop_zero > zero_cutoff) %>% pull(analyte)
keep_analytes <- setdiff(analyte_cols, drop_list)

# Save zero stats and print explicit dropped list
write_csv(zero_stats, file.path(save_dir, "analyte_zero_stats.csv"))
cat("=== ANALYTES DROPPED (more than ", zero_cutoff*100, "% zeros) ===\n", sep = "")
if(length(drop_list)==0) cat("(none)\n") else cat(paste0(drop_list, collapse = "\n"), "\n")
cat("\nKept analytes (for modeling):\n", paste(keep_analytes, collapse = ", "), "\n\n")

# ---------- Build PFAS group sums and TOTAL_PFAS inclusion ----------
# If df already contains TOTAL_PFAS use it; otherwise compute from analyte columns present.
if(!"TOTAL_PFAS" %in% names(df)) {
  message("TOTAL_PFAS not found in df; computing TOTAL_PFAS as sum of all analytes present (kept analytes).")
  df <- df %>% mutate(TOTAL_PFAS = rowSums(select(., all_of(keep_analytes)), na.rm = TRUE))
}

# Build long-format join to compute group sums based on key
long_pf <- df %>%
  select(bird_id = any_of("bird_id"), all_of(keep_analytes)) %>%
  pivot_longer(cols = -bird_id, names_to = "analyte", values_to = "conc") %>%
  left_join(key %>% select(pfas_analyte, chain_cat, func_cat, pfas_age),
            by = c("analyte" = "pfas_analyte"))

# compute sums by bird & group (chain_cat, func_cat, pfas_age)
chain_sums <- long_pf %>% filter(!is.na(chain_cat)) %>%
  group_by(bird_id, chain_cat) %>% summarize(chain_sum = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = chain_cat, values_from = chain_sum, values_fill = 0, names_prefix = "chain_")
func_sums  <- long_pf %>% filter(!is.na(func_cat)) %>%
  group_by(bird_id, func_cat) %>% summarize(func_sum  = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = func_cat,  values_from = func_sum,  values_fill = 0, names_prefix = "func_")
age_sums   <- long_pf %>% filter(!is.na(pfas_age)) %>%
  group_by(bird_id, pfas_age) %>% summarize(age_sum   = sum(conc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pfas_age, values_from = age_sum, values_fill = 0, names_prefix = "age_")

# Merge group sums back into df by bird_id; assume df has bird_id column
if(!"bird_id" %in% names(df)) stop("df must contain a bird_id column for joining group sums.")
df2 <- df %>%
  left_join(chain_sums, by = c("bird_id" = "bird_id")) %>%
  left_join(func_sums,  by = c("bird_id" = "bird_id")) %>%
  left_join(age_sums,   by = c("bird_id" = "bird_id"))

# Replace any NA group sums with 0 (meaning no analytes of that group measured for that bird)
group_cols <- c(grep("^chain_", names(df2), value = TRUE),
                grep("^func_",  names(df2), value = TRUE),
                grep("^age_",   names(df2), value = TRUE))
df2 <- df2 %>% mutate_at(vars(all_of(group_cols)), ~replace_na(., 0))

# ---------- Build full candidate predictor list ----------
pfas_predictors <- keep_analytes           # individual analytes (kept ones)
group_predictors <- group_cols             # e.g., chain_*, func_*, age_*
total_predictor  <- "TOTAL_PFAS"           # already ensured present

# Combine predictors into candidate list (avoid duplicates)
candidate_predictors <- unique(c(pfas_predictors, group_predictors, total_predictor, fixed_covs))

cat("Candidate fixed-effect predictors (count =", length(candidate_predictors), "):\n")
cat(paste(candidate_predictors, collapse = ", "), "\n\n")

# ---------- Prepare data for modeling ----------
# Select only rows with igy (the response) non-missing
if(!"igy" %in% names(df2)) stop("df must contain column 'igy' as the IgY response.")
mod_df <- df2 %>% filter(!is.na(igy))

# Remove rows with NA for any of the fixed covariates (we'll let stepwise pick among PFAS but covariates should be complete)
mod_df <- mod_df %>% filter_at(vars(all_of(fixed_covs)), all_vars(!is.na(.)))

cat("Modeling dataset rows after cleaning:", nrow(mod_df), "\n")

# If too few rows, stop
if(nrow(mod_df) < 20) {
  warning("Fewer than 20 observations remain after cleaning. Stepwise mixed model may not be reliable.")
}

# ---------- Construct global formula ----------
# Build RHS with all candidate fixed-effect predictors joined by +
rhs <- paste(candidate_predictors, collapse = " + ")

if(log_response) {
  response_expr <- paste0("log(igy + ", small_offset, ")")
} else {
  response_expr <- "igy"
}

global_formula <- as.formula(paste0(response_expr, " ~ ", rhs, " + (1 | ", random_effect_col, ")"))
cat("Global model formula:\n")
print(global_formula)

# ---------- Fit global mixed model (lmerTest::lmer for stepwise) ----------
# Note: lmerTest::lmer returns lmertest object enabling lmerTest::step
# Use tryCatch to handle potential convergence or rank-deficiency problems
cat("Fitting global mixed model (this may take some time)...\n")
global_mod <- tryCatch({
  lmerTest::lmer(global_formula, data = mod_df, REML = FALSE, control = lmerControl(check.nobs.vs.nlev = "ignore"))
}, error = function(e){
  stop("Global lmer failed: ", e$message)
})

cat("Global model fitted. Summary:\n")
print(summary(global_mod))

# ---------- Stepwise selection on mixed model ----------
# lmerTest::step performs stepwise selection of fixed effects while keeping random effects structure
# By default it performs backward selection; we run both directions using scope argument derived from global model.
cat("Running stepwise selection (lmerTest::step)...\n")
step_mod <- tryCatch({
  # lmerTest::step will use the model's fixed effects as the starting set and attempt backward selection.
  # For a forward/backward you can provide scope, but with many predictors the default backward is typical.
  lmerTest::step(global_mod, reduce.fixed = TRUE, reduce.random = FALSE, keep = NULL) 
}, error = function(e){
  # If step fails, fall back to manual backward via drop1 with AIC
  message("lmerTest::step failed: ", e$message, "\nFalling back to drop1-based backward AIC on the fitted model.")
  current_mod <- global_mod
  improved <- TRUE
  while(improved) {
    d1 <- drop1(current_mod, test = "Chisq")
    # drop1 returns AIC; find term whose removal reduces AIC the most (lowest AIC)
    # Choose term to drop with highest p-value OR AIC increase? We'll pick based on AIC.
    aic_vals <- d1$AIC; term_names <- rownames(d1)
    # Remove the "(Intercept)" row and any NA entries
    aic_df <- data.frame(term = term_names, AIC = aic_vals, stringsAsFactors = FALSE)
    aic_df <- aic_df[!aic_df$term %in% "(Intercept)", ]
    # if all NA or fewer than 1 candidate, break
    if(all(is.na(aic_df$AIC)) || nrow(aic_df) == 0) break
    # choose term whose removal leads to best (lowest) AIC
    best_row <- which.min(aic_df$AIC)
    # if AIC decreases by at least 2, update model
    aic_current <- AIC(current_mod)
    aic_candidate <- aic_df$AIC[best_row]
    if(!is.na(aic_candidate) && (aic_candidate + 0) < aic_current - 1.999) {
      term_to_drop <- aic_df$term[best_row]
      cat("Dropping term", term_to_drop, "because AIC decreases from", aic_current, "to", aic_candidate, "\n")
      new_formula <- update(formula(current_mod), paste(". ~ . -", term_to_drop))
      current_mod <- update(current_mod, new_formula)
    } else {
      improved <- FALSE
    }
  }
  current_mod
})

cat("Stepwise selection finished. Final model summary:\n")
print(summary(step_mod))

# ---------- Save final model and summary artifacts ----------
saveRDS(global_mod, file = file.path(save_dir, "global_mixed_model.rds"))
saveRDS(step_mod,   file = file.path(save_dir, "stepwise_selected_mixed_model.rds"))

# Extract fixed effects retained
fixed_retained <- fixef(step_mod)
fixed_retained_df <- broom.mixed::tidy(step_mod, effects = "fixed", conf.int = TRUE)
write_csv(fixed_retained_df, file.path(save_dir, "selected_fixed_effects.csv"))
cat("Saved selected fixed effects to:", file.path(save_dir, "selected_fixed_effects.csv"), "\n")

# ---------- Diagnostics: residuals vs fitted and qqplot ----------
# Residuals vs fitted
resid_df <- data.frame(fitted = predict(step_mod), resid = residuals(step_mod))
p1 <- ggplot(resid_df, aes(x = fitted, y = resid)) + geom_point() + geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuals vs Fitted (selected mixed model)", x = "Fitted values", y = "Residuals")
ggsave(file.path(save_dir, "resid_vs_fitted.png"), plot = p1, width = 7, height = 5, dpi = 300)

# QQ plot of residuals
p2 <- ggplot(resid_df, aes(sample = resid)) + stat_qq() + stat_qq_line() +
  labs(title = "QQ-plot of residuals")
ggsave(file.path(save_dir, "qqplot_resid.png"), plot = p2, width = 6, height = 6, dpi = 300)

# Save model summary text
sink(file.path(save_dir, "stepwise_model_summary.txt"))
print(summary(step_mod))
sink()

cat("\nAll outputs saved to directory:", normalizePath(save_dir), "\n")
cat("Dropped analytes list saved earlier to:", file.path(save_dir, "analyte_zero_stats.csv"), "\n\n")

# Final quick-report printing
cat("==== QUICK REPORT ====\n")
cat("Dropped analytes (> ", zero_cutoff*100, "% zeros):\n", sep = "")
if(length(drop_list)==0) cat("(none)\n") else cat(paste(drop_list, collapse = "\n"), "\n")
cat("\nSelected fixed effects (from CSV):\n")
print(fixed_retained_df)

# End script
