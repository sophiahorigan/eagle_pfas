rm(list = ls())

# ------------------- Setup -------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(nlme)
library(car)

setwd("~/Desktop/EagleStats/")

# ------------------- Load Data -------------------
data <- read.csv("./input/2023_2024_2025_AllData.csv")
key  <- read.csv("./input/PFAS_subtypes.csv")

# Only keep birds with IgY
data <- data %>% filter(!is.na(igy))



# Check which PFAS analytes are in the dataset
pfas_cols_in_data <- intersect(key$pfas_analyte, names(data))

# Optional: show which are missing
missing <- setdiff(key$pfas_analyte, names(data))
cat("PFAS not in data:\n")
print(missing)

# Pivot only existing columns
long_pf <- data %>%
  pivot_longer(cols = all_of(pfas_cols_in_data),
               names_to = "pfas_analyte", values_to = "conc") %>%
  left_join(key, by = "pfas_analyte")

# Chain sums
chain_sums <- long_pf %>% 
  filter(!is.na(chain_cat)) %>%
  group_by(bird_id, chain_cat) %>%
  summarize(chain_sum = sum(conc, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = chain_cat, values_from = chain_sum, 
              values_fill = 0, names_prefix="chain_")

# Function sums
func_sums <- long_pf %>% 
  filter(!is.na(func_cat)) %>%
  group_by(bird_id, func_cat) %>%
  summarize(func_sum = sum(conc, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = func_cat, values_from = func_sum,
              values_fill = 0, names_prefix="func_")

# Age sums
age_sums <- long_pf %>% 
  filter(!is.na(pfas_age)) %>%
  group_by(bird_id, pfas_age) %>%
  summarize(age_sum = sum(conc, na.rm=TRUE), .groups="drop") %>%
  pivot_wider(names_from = pfas_age, values_from = age_sum,
              values_fill = 0, names_prefix="age_")

# Merge sums into main dataset
data2 <- data %>%
  left_join(chain_sums, by="bird_id") %>%
  left_join(func_sums, by="bird_id") %>%
  left_join(age_sums, by="bird_id")

# ------------------- Transformations -------------------
# Identify PFAS columns (raw + summed)
pfas_cols <- c(key$pfas_analyte, 
               grep("^chain_|^func_|^age_", names(data2), value = TRUE))

# Proportion columns
prop_cols <- c("Lperc","Hperc","Eperc","Mperc","Bperc")

# Apply transformations and track
transforms <- list()
for(col in c(pfas_cols, prop_cols)){
  if(!col %in% names(data2)) next
  x <- data2[[col]]
  if(all(is.na(x) | x == x[1])){  # skip if no variation
    cat("Skipping", col, "- insufficient variation\n")
    next
  }
  # Proportions: arcsine sqrt
  if(col %in% prop_cols){
    data2[[paste0(col,"_T")]] <- asin(sqrt(x))
    transforms[[col]] <- "arcsine"
  } else {  # PFAS: log10
    # Replace non-detects with zero (already done) and add small constant for log
    data2[[paste0(col,"_T")]] <- log10(x + 1e-9)
    transforms[[col]] <- "log10"
  }
}

# ------------------- Modeling -------------------
response <- "igy"
predictors <- c(pfas_cols, prop_cols)
results_list <- list()

for(pred in predictors){
  pred_T <- paste0(pred,"_T")
  if(!pred_T %in% names(data2)) next
  df <- data2 %>% 
    select(all_of(c(response,"nest_no",pred_T))) %>% 
    filter(!is.na(.data[[response]]) & !is.na(.data[[pred_T]]))
  
  if(length(unique(df[[pred_T]])) < 2){
    cat("Skipping predictor", pred_T, "- insufficient variation\n")
    next
  }
  
  # Linear mixed-effects model
  form <- as.formula(paste(response,"~",pred_T))
  mod <- try(lme(form, random = ~1|nest_no, data=df), silent=TRUE)
  
  if(inherits(mod,"try-error")){
    cat("Model failed for", pred_T,"\n")
    next
  }
  
  # Residuals normality
  sw <- shapiro.test(residuals(mod, level=0))
  normal_resid <- sw$p.value > 0.05
  
  # Correlation
  cor_test <- if(normal_resid){
    cor.test(df[[pred_T]], df[[response]], method="pearson")
  } else {
    cor.test(df[[pred_T]], df[[response]], method="spearman")
  }
  
  # Store results
  results_list[[pred]] <- list(
    model = mod,
    residual_normal = normal_resid,
    cor_method = if(normal_resid) "pearson" else "spearman",
    cor_estimate = cor_test$estimate,
    cor_p = cor_test$p.value
  )
  
  # Plot
  p <- ggplot(df, aes_string(x=pred_T, y=response)) +
    geom_point(alpha=0.6) +
    geom_smooth(method="lm", se=TRUE, color="blue") +
    labs(title=paste(response,"vs",pred_T),
         subtitle=paste("Residuals normal?", normal_resid, "| Cor:", results_list[[pred]]$cor_method)) +
    theme_minimal()
  ggsave(paste0("PFAS_model_plots/",pred_T,"_plot.jpeg"), plot=p, width=6, height=5, dpi=300)
}

# ------------------- Save Summary -------------------
summary_df <- do.call(rbind, lapply(names(results_list), function(pred){
  r <- results_list[[pred]]
  data.frame(
    Predictor = pred,
    Residuals_Normal = r$residual_normal,
    Cor_Method = r$cor_method,
    Cor_Estimate = r$cor_estimate,
    Cor_p = r$cor_p
  )
}))

write.csv(summary_df, "PFAS_model_plots/IgY_PFAS_summary.csv", row.names=FALSE)
