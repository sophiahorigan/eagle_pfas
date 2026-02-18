##########################################
# PFAS vs IgY (humoral immunity)
# Which pfas are most strongly associated with humoral immunity?
##########################################

rm(list=ls())

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)


input_file <- "./input/2023_2024_2025_AllData.csv"
output_dir <- "./output/3_PFAS_vs_IgY"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(log_igy = log(igy + 1e-6))

# Identify PFAS columns
pfas_cols <- df %>% select(starts_with("PF"), TOTAL_PFAS, PCA_score) %>% colnames()
pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal", "DNA_sex")]

# Covariates
covariates <- c("DNA_sex", "age")
nest_var <- "nest_no"

results_list <- list()

for(pfas in pfas_cols){
  
  if(!is.numeric(df[[pfas]]) || sum(!is.na(df[[pfas]])) < 6) next
  
  df_model <- df %>% filter(!is.na(.data[[pfas]]) & !is.na(log_igy)) %>% mutate(PFAS_scaled = scale(.data[[pfas]]))
  
  formula_text <- paste0("log_igy ~ PFAS_scaled + ", paste(covariates, collapse=" + "), " + (1|", nest_var,")")
  
  model <- try(lmer(as.formula(formula_text), data=df_model, REML=FALSE), silent=TRUE)
  
  if(inherits(model, "try-error")) next
  
  tidy_res <- broom.mixed::tidy(model) %>% filter(term=="PFAS_scaled") %>%
    mutate(pfas_var = pfas, immune_var = "IgY", n_obs = nrow(df_model))
  
  results_list[[pfas]] <- tidy_res
}

final_results <- bind_rows(results_list) %>% mutate(FDR_p = p.adjust(p.value, method="fdr"))

write.csv(final_results, file=file.path(output_dir,"PFAS_IgY_mixed_model_results.csv"), row.names=FALSE)

cat("Done: PFAS vs IgY\n")
