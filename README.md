README: Eagle PFAS & Immune Function Scripts

This repository contains R scripts for analyzing relationships between PFAS exposure and immune function, morphometrics, and nest/geographic effects in your bird dataset (2023–2025). Each script is standalone and outputs CSVs and/or plots to corresponding folders in the ./output/ directory.

Script Overview
0_eda_igy_pfas.R

Purpose: Exploratory data analysis (EDA) for IgY and PFAS.
Question: What are the distributions of IgY and PFAS variables, and are there obvious patterns, outliers, or missing data?
Output: Histograms, scatterplots, summary statistics.

0_Misc_Figs_PFAS.R

Purpose: Generate miscellaneous PFAS figures used in reports.
Question: Provide visual summaries of PFAS distributions or correlations that don’t fit into other scripts.
Output: PNG/JPEG plots in ./output/misc_figs/.

1_igy_PFAS_relationship.R

Purpose: Examine the shape of relationships between IgY (humoral immunity) and individual or grouped PFAS.
Question: What is the functional form (linear, log, quadratic, exponential, power) of IgY-PFAS relationships?
Output:

Plots of IgY vs each PFAS with all candidate model fits.

Plots of the “best-fit” model per PFAS.

Stored in ./output/1_igy_pfas_relationship/.

2_clysis_PFAS_relationship.R

Purpose: Examine the shape of relationships between clysis (cellular immunity) and individual or grouped PFAS.
Question: What is the functional form (linear, log, quadratic, exponential, power) of clysis-PFAS relationships?
Output:

Plots of clysis vs each PFAS with all candidate model fits.

Plots of the “best-fit” model per PFAS.

Stored in ./output/2_clysis_pfas_relationship/.

3_pfas_vs_igy.R

Purpose: Test whether PFAS predict IgY after accounting for covariates (age, sex) and random nest effects.
Question: Which PFAS are significantly associated with humoral immunity?
Output:

CSV with model estimates, p-values, FDR-adjusted p-values.

Plots of predicted vs observed IgY.

Stored in ./output/3_pfas_vs_igy/.

4_clysis_vs_PFAS.R

Purpose: Test whether PFAS predict clysis after accounting for covariates (age, sex) and random nest effects.
Question: Which PFAS are significantly associated with cellular immunity?
Output:

CSV with model estimates, p-values, FDR-adjusted p-values.

Plots of predicted vs observed clysis.

Stored in ./output/4_clysis_vs_PFAS/.

5_PFAS_interactions.R

Question: Does sex or age modify PFAS-immune relationships?
Purpose: Fit mixed models for IgY and clysis, testing interactions between each PFAS and sex (DNA_sex) and age while accounting for random nest effects (nest_no).
Outputs:

CSV with model estimates, p-values, and FDR-adjusted p-values for all PFAS × sex and PFAS × age interaction terms.

Saved to ./output/5_PFAS_interactions/ as PFAS_interactions.csv.

6_geo_interactions.R

Purpose: Test nest-level and geographic effects on PFAS exposure and immune function.
Question: Do nest identity, county, or latitude/longitude predict PFAS levels or immune function?
Output:

CSV with fixed effect estimates and p-values (FDR-corrected).

Stored in ./output/6_Nest_Geo_Effects/.

7_whiteblood_PFAS_corr.R

Purpose: Examine correlations between PFAS and white blood cell measures.
Question: Are proportions of lymphocytes (Lperc), heterophils (Hperc), eosinophils (Eperc), monocytes (Mperc), basophils (Bperc), or total WBC count associated with PFAS exposure?
Output:

CSV with Spearman correlation coefficients, p-values, and FDR-adjusted p-values.

Stored in ./output/7_PFAS_WBC_Correlations/.