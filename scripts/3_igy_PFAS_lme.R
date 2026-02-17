rm(list = ls())

# ------------------- Setup -------------------
setwd("~/Desktop/EagleStats/")

# Required packages
pkgs <- c("dplyr","tidyr","lme4","lmerTest","broom.mixed","readr","ggplot2")
missing_pkgs <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if(length(missing_pkgs)) install.packages(missing_pkgs)
lapply(pkgs, library, character.only = TRUE)

# Input/output
df <- read.csv("./input/2023_2024_2025_AllData.csv")
key <- read.csv("./input/PFAS_subtypes.csv")
save_dir <- "./output/pfas_igy_significance_results"
dir.create(save_dir, showWarnings = FALSE)

# ------------------- User options -------------------
zero_cutoff <- 0.50        # drop analytes with >50% zeros
small_offset <- 1e-6       # offset for log
log_response <- TRUE        # log-transform igy
random_effect_col <- "nest_name"
fixed_covs <- c("lat","long","weight_kg","age")  # always included

# ------------------- Identify analytes -------------------
analyte_cols <- intersect(names(df), key$pfas_analyte)
if(length(analyte_cols) == 0) stop("No analyte columns found!")

# Drop analytes with too many zeros
zero_stats <- sapply(df[analyte_cols], function(x) mean(x==0, na.rm = TRUE))
drop_list <- names(zero_stats[zero_stats > zero_cutoff])
keep_analytes <- setdiff(analyte_cols, drop_list)

# ------------------- Compute group sums -------------------
long_pf <- df %>%
  select(bird_id = any_of("bird_id"), all_of(keep_analytes)) %>%
  pivot_longer(-bird_id, names_to = "analyte", values_to = "conc") %>%
  left_join(key %>% select(pfas_analyte, chain_cat, func_cat, pfas_age),
            by = c("analyte" = "pfas_analyte"))

# Summarize by group
chain_sums <- long_pf %>% filter(!is.na(chain_cat)) %>%
  group_by(bird_id, chain_cat) %>% summarize(chain_sum = sum(conc, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=chain_cat, values_from=chain_sum, values_fill=0, names_prefix="chain_")

func_sums <- long_pf %>% filter(!is.na(func_cat)) %>%
  group_by(bird_id, func_cat) %>% summarize(func_sum = sum(conc, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=func_cat, values_from=func_sum, values_fill=0, names_prefix="func_")

age_sums <- long_pf %>% filter(!is.na(pfas_age)) %>%
  group_by(bird_id, pfas_age) %>% summarize(age_sum = sum(conc, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from=pfas_age, values_from=age_sum, values_fill=0, names_prefix="age_")

# Merge sums back into df
df2 <- df %>%
  left_join(chain_sums, by="bird_id") %>%
  left_join(func_sums,  by="bird_id") %>%
  left_join(age_sums,   by="bird_id")

# Replace NA in sums with 0
group_cols <- c(grep("^chain_", names(df2), value=TRUE),
                grep("^func_",  names(df2), value=TRUE),
                grep("^age_",   names(df2), value=TRUE))
df2[group_cols][is.na(df2[group_cols])] <- 0

# ------------------- Clean modeling dataframe -------------------
if(!"igy" %in% names(df2)) stop("Column 'igy' missing.")
mod_df <- df2 %>%
  filter(!is.na(igy)) %>%
  filter_at(vars(all_of(fixed_covs)), all_vars(!is.na(.)))  # remove NA in covariates

# Convert character to factor
mod_df_clean <- mod_df %>%
  mutate(across(where(is.character), as.factor))

cat("Modeling dataset rows:", nrow(mod_df_clean), "\n")

# ------------------- Function to fit single predictor model -------------------
fit_single_predictor <- function(pred) {
  response <- if(log_response) paste0("log(igy + ", small_offset, ")") else "igy"
  fmla <- as.formula(paste0(response, " ~ ", pred, " + ", paste(fixed_covs, collapse=" + "),
                            " + (1 | ", random_effect_col, ")"))
  mod <- tryCatch(lmerTest::lmer(fmla, data=mod_df_clean, REML=FALSE), error=function(e) NULL)
  if(is.null(mod)) return(NULL)
  tidy_mod <- broom.mixed::tidy(mod, effects="fixed", conf.int=TRUE)
  tidy_mod$predictor_tested <- pred
  return(tidy_mod)
}

# ------------------- Loop over all analytes and group sums -------------------
all_preds <- c(keep_analytes, group_cols)
results_list <- lapply(all_preds, fit_single_predictor)
results_list <- results_list[!sapply(results_list, is.null)]
all_results <- bind_rows(results_list) %>%
  select(predictor_tested, term, estimate, std.error, statistic, p.value, conf.low, conf.high)

# ------------------- Save results -------------------
write_csv(all_results, file.path(save_dir, "single_predictor_results.csv"))
cat("Saved results to:", file.path(save_dir, "single_predictor_results.csv"), "\n")

# ------------------- Quick report of significant predictors -------------------
sig_results <- all_results %>% filter(p.value < 0.05)
cat("Significant predictors (p < 0.05):\n")
if(nrow(sig_results)==0) cat("(none)\n") else print(sig_results[,c("predictor_tested","term","estimate","p.value")])

# ------------------- Plot significant predictors -------------------
if(nrow(sig_results) > 0) {
  p_sig <- ggplot(sig_results, aes(x = reorder(predictor_tested, estimate), y = estimate, fill = estimate > 0)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
    coord_flip() +
    scale_fill_manual(values = c("red", "blue"), labels = c("Negative", "Positive")) +
    labs(title = "Significant predictors of IgY",
         x = "Predictor",
         y = "Estimate (with 95% CI)",
         fill = "Effect direction") +
    theme_minimal(base_size = 12)
  
  print(p_sig)
  
  ggsave(file.path(save_dir, "significant_predictors_plot.png"), plot = p_sig, width = 7, height = 6, dpi = 300)
  cat("Saved plot of significant predictors to:", file.path(save_dir, "significant_predictors_plot.png"), "\n")
} else {
  cat("No significant predictors to plot.\n")
}
