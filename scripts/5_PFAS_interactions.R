##########################################
# Interaction: PFAS × sex/age
##########################################

rm(list=ls())
library(tidyverse)
library(lme4)
library(lmerTest)
library(broom.mixed)

input_file <- "./input/2023_2024_2025_AllData.csv"
output_dir <- "./output/5_PFAS_interactions"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

df <- read.csv(input_file, stringsAsFactors = FALSE) %>%
  mutate(
    log_igy = log(igy + 1e-6),
    log_clysis = log(clysis + 1e-6)
  )

pfas_cols <- df %>%
  select(starts_with("PF"), TOTAL_PFAS, PCA_score) %>%
  colnames()

pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal", "DNA_sex")]

nest_var <- "nest_no"

immune_vars <- c("log_igy","log_clysis")

results_list <- list()

for(immune in immune_vars){
  for(pfas in pfas_cols){
    
    if(!is.numeric(df[[pfas]]) || sum(!is.na(df[[pfas]])) < 8) next
    
    df_model <- df %>%
      filter(!is.na(.data[[pfas]]),
             !is.na(.data[[immune]]),
             !is.na(DNA_sex),
             !is.na(age),
             !is.na(.data[[nest_var]])) %>%
      mutate(PFAS_scaled = as.numeric(scale(.data[[pfas]])))
    
    formula_text <- paste0(
      immune,
      " ~ PFAS_scaled * (DNA_sex + age) + (1|", nest_var, ")"
    )
    
    model <- try(lmer(as.formula(formula_text),
                      data=df_model,
                      REML=FALSE),
                 silent=TRUE)
    
    if(inherits(model,"try-error")) next
    
    tidy_res <- broom.mixed::tidy(model, effects="fixed") %>%
      mutate(
        pfas_var = pfas,
        immune_var = immune,
        n_obs = nrow(df_model)
      )
    
    results_list[[paste(immune,pfas)]] <- tidy_res
  }
}

final_results <- bind_rows(results_list)

# Separate corrections
final_results <- final_results %>%
  mutate(
    test_type = case_when(
      grepl(":", term) ~ "interaction",
      term == "PFAS_scaled" ~ "main_effect",
      TRUE ~ "covariate"
    )
  ) %>%
  group_by(test_type) %>%
  mutate(FDR_p = p.adjust(p.value, method="BH")) %>%
  ungroup()

write.csv(final_results,
          file=file.path(output_dir,"PFAS_interactions.csv"),
          row.names=FALSE)

print("Interaction modeling complete.")