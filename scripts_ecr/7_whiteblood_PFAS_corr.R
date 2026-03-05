##########################################
# Correlations: PFAS vs White Blood Cells
##########################################

rm(list=ls())
library(tidyverse)

# --- Settings ---
input_file <- "./input/2023_2024_2025_AllData_clean.csv"
output_dir <- "./output/7_PFAS_WBC_Correlations"
dir.create(output_dir, showWarnings=FALSE, recursive=TRUE)

# --- Load data ---
df <- read.csv(input_file, stringsAsFactors=FALSE)

# --- PFAS variables ---
pfas_cols <- df %>% select(starts_with("PF"), TOTAL_PFAS) %>% colnames()
pfas_cols <- pfas_cols[!pfas_cols %in% c("PCR_avmal","DNA_sex")]

# --- WBC variables ---
wbc_props <- c("Lperc","Hperc","Eperc","Mperc","Bperc")
wbc_totals <- c("Total_count")

wbc_vars <- c(wbc_props, wbc_totals)

# --- Prepare results storage ---
results_list <- list()

# --- Loop through all PFAS x WBC combinations ---
for(pfas in pfas_cols){
  for(wbc in wbc_vars){
    df_pair <- df %>% select(all_of(c(pfas,wbc))) %>% filter(!is.na(.data[[pfas]]) & !is.na(.data[[wbc]]))
    if(nrow(df_pair) < 6) next  # skip tiny samples
    cor_test <- try(cor.test(df_pair[[pfas]], df_pair[[wbc]], method="spearman"), silent=TRUE)
    if(inherits(cor_test,"try-error")) next
    results_list[[paste(pfas,wbc)]] <- data.frame(
      PFAS = pfas,
      WBC_var = wbc,
      rho = cor_test$estimate,
      p.value = cor_test$p.value,
      n = nrow(df_pair)
    )
  }
}

# --- Combine results ---
results_df <- bind_rows(results_list)

# --- FDR adjustment ---
results_df <- results_df %>%
  mutate(FDR_p = p.adjust(p.value, method="fdr"))

# --- Save output ---
write.csv(results_df, file=file.path(output_dir,"PFAS_WBC_correlations.csv"), row.names=FALSE)

cat("Done: PFAS vs WBC correlations\n")
