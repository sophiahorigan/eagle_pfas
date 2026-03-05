##########################################
# Nest-level / geographic effects
##########################################

rm(list=ls())

library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)

input_file <- "./input/2023_2024_2025_AllData_clean.csv"
output_dir <- "./output/6_Nest_Geo_Effects"
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

df <- read.csv(input_file, stringsAsFactors=FALSE) %>%
  mutate(log_igy=log(igy+1e-6),
         log_clysis=log(clysis+1e-6))

pfas_cols <- df %>% select(starts_with("PF"),TOTAL_PFAS,PCA_score) %>% colnames()
pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal","DNA_sex")]

# Nest-level/geography predictors
geo_vars <- c("nest_no","county","lat","long")

results_list <- list()

# PFAS response
for(pfas in pfas_cols){
  for(pred in geo_vars){
    if(!is.numeric(df[[pfas]]) & !is.factor(df[[pred]])) next
    df_model <- df %>% filter(!is.na(.data[[pfas]]) & !is.na(.data[[pred]]))
    formula_text <- paste0(pfas," ~ ",pred," + (1|nest_no)")
    model <- try(lmer(as.formula(formula_text), data=df_model, REML=FALSE), silent=TRUE)
    if(inherits(model,"try-error")) next
    tidy_res <- broom.mixed::tidy(model) %>% filter(term==pred) %>% mutate(response_var=pfas, predictor=pred, n_obs=nrow(df_model))
    results_list[[paste(pfas,pred)]] <- tidy_res
  }
}

# Immune response
immune_vars <- c("log_igy","log_clysis")
for(immune in immune_vars){
  for(pred in geo_vars){
    df_model <- df %>% filter(!is.na(.data[[immune]]) & !is.na(.data[[pred]]))
    formula_text <- paste0(immune," ~ ",pred," + (1|nest_no)")
    model <- try(lmer(as.formula(formula_text), data=df_model, REML=FALSE), silent=TRUE)
    if(inherits(model,"try-error")) next
    tidy_res <- broom.mixed::tidy(model) %>% filter(term==pred) %>% mutate(response_var=immune, predictor=pred, n_obs=nrow(df_model))
    results_list[[paste(immune,pred)]] <- tidy_res
  }
}

final_results <- bind_rows(results_list) %>% mutate(FDR_p=p.adjust(p.value, method="fdr"))
write.csv(final_results, file=file.path(output_dir,"Nest_Geo_vs_PFAS_Immune.csv"), row.names=FALSE)
cat("Done: Nest/geography vs PFAS and immune function\n")
