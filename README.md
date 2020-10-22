# multiomics_ELGAN
Code for the working paper: ["Evidence for the Placenta-Brain Axis: Multi-Omic Kernel Aggregation Predicts Intellectual and Social Impairment in Children Born Extremely Preterm" (Santos and Bhattacharya et al, 2020)](https://www.biorxiv.org/content/10.1101/2020.07.19.211029v2)

- kernel_regressions.R provides functions that were used for linear kernel, Gaussian kernel, and elastic net regression for outcomes
- EWAS_example.R provides functions that were used to conduct EWAS to find outcome-associated CpG site
- deg_example.R provides functions for differential expression analysis for either mRNA or miRNA expression.
- power_calculations.R provides functions and code for sensitivity analyses for power
- power_linear_regression.tsv provides power results from simulations for the multiple regression of SRS/IQ versus demographics
- KernelModels.xlsx provides point estimates from KRLS models on full data
- Unmix_Deconvolution.R provides general code for reference-based deconvolution using unmix
