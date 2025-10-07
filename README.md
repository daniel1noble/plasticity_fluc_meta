# Temperature variability has limited effects on phenotypic plasticity in ectotherms – a meta-analysis

[![DOI](https://zenodo.org/badge/979073901.svg)](https://doi.org/10.5281/zenodo.17290563)

This repository contains data and R code for the associated manuscript.  

Clayton W. Stocker, Stephanie M. Bamford, Miki Jahn, Geoffrey P. F. Mazué, Amanda K. Pettersen, Daniel Ritchie, Alexander Rubin, Daniel W.A. Noble*, Frank Seebacher*. (2025) Temperature variability has limited effects on phenotypic plasticity in ectotherms – a meta-analysis. Ecology Letters.

Please email frank.seebacher@sydney.edu.au or daniel.noble@anu.edu.au for any inquiries. A release is also lodged on [Zenodo](https://doi.org/10.5281/zenodo.17290563)

## How to use this repository?

1. Clone or download the repository to your local machine.
2. Open the `Plasticity_Fluctuation_Meta.Rproj` file which will open RStudio and set your working directory.
3. Ensure you have all the required R packages installed (see Software Dependencies section below).
4. Run through the `Analysis_Final.R` script to reproduce the analyses and figures from the manuscript. Note that you do not need to re-run the entire analysis if you use the saved model objects in the `output/models/` directory. But, you need to set `rerun` = FALSE in the `Analysis_Final.R` script to use the saved models.

There are additional scripts and resources in this repository that are important and are detailed below:

- `Complexity_Final_data.csv`: This is the final processed dataset used for the meta-analysis.
- `complexity_tree.nwk`: This file contains the phylogenetic tree used in the analysis.
- `func.R`: This script contains custom functions used in the analysis.
- `output/figs/`: This directory contains all figures generated from the analysis.
- `output/tables/`: This directory contains all tables generated from the analysis.
- `output/models/`: This directory contains saved model objects from the analysis. Use these so you don't have to re-run the entire analysis.

## Analyses

All code for analyses reproducing the results and figures is found in `Analysis_Final.R` script. The code also generates supplementary tables and results, which are saved to the `output/` directory.

## Data
The data is the `Complex_Final_Data.csv` file. This file contains the data used in the meta-analysis. The meta-data that describes the columns are as follows. Note that some columns are created during the analysis (e.g., transformed means and SDs, effect sizes).

| Column Name                        | Description |
|-----------------------------------|-------------|
| Study_ID                           | Unique identifiers for each paper. |
| Species_ID                         | Identifier for each species that a study investigates. |
| Treatment_ID_T1                    | Identifier for the studies’ first constant and fluctuating temperature treatment pair (fluctuating mean matching constant temperature). |
| Treatment_ID_T2                    | Identifier for the studies’ second constant and fluctuating temperature treatment pair (fluctuating mean matching constant temperature). |
| Trait_ID                           | Identifier for each phenotypic trait a study measures. |
| First_Author                       | Initials and surname for the first author of the study. |
| Title                              | Title of the study. |
| Year                               | Year of publication. |
| Journal                            | Name of the journal the study was published in. |
| Journal_Impact_Factor-2021         | Journal impact factor as of 2021 (most recent records at time of data collection). |
| Kingdom                            | Kingdom for each species that a study investigates. |
| Phylum                             | Phylum for each species that a study investigates. |
| Class                              | Class for each species that a study investigates. |
| Order                              | Order for each species that a study investigates. |
| Family                             | Family for each species that a study investigates. |
| Scientific_Name                    | Binomial nomenclature for species. |
| Ecosystem                          | Ecosystem the species is naturally observed (Aquatic and Terrestrial). Amphibians considered aquatic. |
| Plasticity_Mechanism               | Mechanism of plasticity the treatments were exposed to (Acclimation or Developmental Plasticity). |
| Developmental_Exposure_Time_Category | Categorisation of Developmental Exposure Time. |
| Developmental_Exposure_Time        | The period of exposure for treatments imposed during development. |
| Acclimation_Exposure_Time          | The duration of exposure for acclimation treatments. |
| Exposure_Units                     | Units of Acclimation Exposure Time (Days). |
| T1_constant                        | Temperature of the first constant temperature treatment. |
| T1_fluctuation                     | Mean temperature of the first fluctuating temperature treatment. |
| T2_constant                        | Temperature of the second constant temperature treatment. |
| T2_fluctuation                     | Mean temperature of the second fluctuating temperature treatment. |
| Fluctuation_Magnitude              | Amplitude of the two fluctuating temperature treatments. |
| Fluctuation_Category               | Type of fluctuations imposed (Sinusoidal, Alternating, Stepwise, Stochastic). |
| Fluctuation_Period                 | Period of one fluctuation oscillation. |
| Fluctuation_Unit                   | Units of Fluctuation Period (Days). |
| Number_Of_Fluctuations             | Acclimation Exposure Time/Fluctuation Period for acclimation treatments. |
| Acclimation_Life-History_Stage     | Life-history stage of organisms for acclimation treatments. |
| Acclimation_Life-History_Stage_Category | Categorisation of Acclimation Life-history Stages. |
| Trait_Category                     | Categorisation of Trait Measurements. |
| Measurement                        | Phenotypic traits measured following treatment exposure. |
| Trait_Unit                         | Units for measurements. |
| per_transform                      | Indicator whether means and SDs need to be transformed from proportion scale (“Yes” or “No”). |
| ln_transform                       | Indicator whether means and SDs need to be transformed from log scale (“Yes” or “No”). |
| Sex                                | Sex of the organisms investigated (Both, Female or Male). NA = not specified. |
| Performance_Curve                  | Whether a performance curve was recorded in the study (Yes, No). |
| Complex_Design                     | Whether a comparison between constant and fluctuating treatments was made at multiple temperatures (Yes). |
| Shared_Control_Number              | Unique identifier for shared control codes across effect sizes. |
| Shared_Animal_Number               | Unique identifier for shared animal codes across effect sizes. |
| n_T1_C                             | Sample size of the constant treatment (first pair). |
| Mean_T1_C                          | Mean response of the constant treatment (first pair). |
| SD_Final_T1_C                      | Standard deviation of the constant treatment (first pair). |
| n_T1_F                             | Sample size of the fluctuating treatment (first pair). |
| Mean_T1_F                          | Mean response of the fluctuating treatment (first pair). |
| SD_Final_T1_F                      | Standard deviation of the fluctuating treatment (first pair). |
| n_T2_C                             | Sample size of the constant treatment (second pair). |
| Mean_T2_C                          | Mean response of the constant treatment (second pair). |
| SD_Final_T2_C                      | Standard deviation of the constant treatment (second pair). |
| n_T2_F                             | Sample size of the fluctuating treatment (second pair). |
| Mean_T2_F                          | Mean response of the fluctuating treatment (second pair). |
| SD_Final_T2_F                      | Standard deviation of the fluctuating treatment (second pair). |
| obs                                | Unique observation identifier (for random effect). |
| phylo                              | Open Tree of Life identifier used for phylogenetic random effect. |
| vert_invert                        | Whether species is “Invertebrate” or “Vertebrate” (for MLMR models). |
| Mean_T1_C_trans                    | Corrected mean of constant T1 after transformation. |
| Mean_T1_F_trans                    | Corrected mean of fluctuating T1 after transformation. |
| Mean_T2_C_trans                    | Corrected mean of constant T2 after transformation. |
| Mean_T2_F_trans                    | Corrected mean of fluctuating T2 after transformation. |
| SD_Final_T1_C_trans                | Corrected SD of constant T1 after transformation. |
| SD_Final_T1_F_trans                | Corrected SD of fluctuating T1 after transformation. |
| SD_Final_T2_C_trans                | Corrected SD of constant T2 after transformation. |
| SD_Final_T2_F_trans                | Corrected SD of fluctuating T2 after transformation. |
| PRRD                               | Uncorrected plasticity response ratio difference (PRRD). |
| v_PRRD                             | Sampling variance for PRRD. |
| lnRR1                              | Log response ratio of high vs low temperature in constant treatment. |
| lnRR2                              | Log response ratio of high vs low temperature in fluctuating treatment. |
| combined                           | Reaction norm direction (Positive, Negative, Opposite). |
| meaning                            | Whether slope is steeper in constant (C) or fluctuating (F) treatment. |
| PRRD_cor                           | Corrected PRRD with consistent sign across studies. |
| Year_Z                             | Z-transformed year of publication. |
| inv_n_eff                          | Inverse of the effective sample size (publication bias metric). |
| sqrt_inv_n_eff                     | Square root of inverse effective sample size (publication bias metric). |

## Software Dependencies

R version 4.5.1 (2025-06-13)
Platform: aarch64-apple-darwin20
Running under: macOS Tahoe 26.0.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

locale:
[1] C.UTF-8/C.UTF-8/C.UTF-8/C/C.UTF-8/C.UTF-8

time zone: Australia/Sydney
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] rphylopic_1.5.0     janitor_2.2.1       gridExtra_2.3      
 [4] ggbeeswarm_0.7.2    ggridges_0.5.6      ggpmisc_0.6.2      
 [7] ggpp_0.5.9          robumeta_2.1        orchaRd_2.0        
[10] metaAidR_0.0.0.9000 MCMCglmm_2.36       coda_0.19-4.1      
[13] phytools_2.4-4      maps_3.4.3          flextable_0.9.9    
[16] brms_2.22.0         Rcpp_1.0.14         metafor_4.8-0      
[19] numDeriv_2016.8-1.1 metadat_1.4-0       Matrix_1.7-3       
[22] latex2exp_0.9.6     patchwork_1.3.1     emmeans_1.11.1     
[25] ape_5.8-1           DescTools_0.99.60   rotl_3.1.0         
[28] gtsummary_2.3.0     readxl_1.4.5        lubridate_1.9.4    
[31] forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4        
[34] purrr_1.0.4         readr_2.1.5         tidyr_1.3.1        
[37] tibble_3.3.0        ggplot2_3.5.2       tidyverse_2.0.0    

loaded via a namespace (and not attached):
  [1] cubature_2.1.4          splines_4.5.1           cellranger_1.1.0       
  [4] XML_3.99-0.18           lifecycle_1.0.4         doParallel_1.0.17      
  [7] lattice_0.22-7          MASS_7.3-65             backports_1.5.0        
 [10] magrittr_2.0.3          rmarkdown_2.29          grImport2_0.3-3        
 [13] zip_2.3.3               askpass_1.2.1           gld_2.6.7              
 [16] pbapply_1.7-2           RColorBrewer_1.1-3      multcomp_1.4-28        
 [19] abind_1.4-8             expm_1.0-0              quadprog_1.5-8         
 [22] TH.data_1.1-3           tensorA_0.36.2.1        sandwich_3.1-1         
 [25] gdtools_0.4.2           rentrez_1.2.4           MatrixModels_0.5-4     
 [28] bridgesampling_1.1-2    codetools_0.2-20        xml2_1.3.8             
 [31] tidyselect_1.2.1        bayesplot_1.13.0        farver_2.1.2           
 [34] matrixStats_1.5.0       base64enc_0.1-3         mathjaxr_1.8-0         
 [37] jsonlite_2.0.0          e1071_1.7-16            survival_3.8-3         
 [40] iterators_1.0.14        systemfonts_1.2.3       foreach_1.5.2          
 [43] tools_4.5.1             progress_1.2.3          ragg_1.4.0             
 [46] glue_1.8.0              mnormt_2.1.1            xfun_0.52              
 [49] distributional_0.5.0    loo_2.8.0               withr_3.0.2            
 [52] combinat_0.0-8          fastmap_1.2.0           boot_1.3-31            
 [55] SparseM_1.84-2          openssl_2.3.3           digest_0.6.37          
 [58] timechange_0.3.0        R6_2.6.1                estimability_1.5.1     
 [61] colorspace_2.1-1        textshaping_1.0.1       rsvg_2.6.2             
 [64] jpeg_0.1-11             generics_0.1.4          fontLiberation_0.1.0   
 [67] data.table_1.17.6       corpcor_1.6.10          class_7.3-23           
 [70] clusterGeneration_1.3.8 prettyunits_1.2.0       httr_1.4.7             
 [73] scatterplot3d_0.3-44    pkgconfig_2.0.3         gtable_0.3.6           
 [76] Exact_3.3               htmltools_0.5.8.1       fontBitstreamVera_0.1.1
 [79] scales_1.4.0            lmom_3.2                png_0.1-8              
 [82] posterior_1.6.1         snakecase_0.11.1        knitr_1.50             
 [85] rstudioapi_0.17.1       tzdb_0.5.0              rncl_0.8.7             
 [88] uuid_1.2-1              checkmate_2.3.2         nlme_3.1-168           
 [91] curl_6.4.0              proxy_0.4-27            zoo_1.8-14             
 [94] DEoptim_2.2-8           rootSolve_1.8.2.4       parallel_4.5.1         
 [97] vipor_0.4.7             pillar_1.10.2           grid_4.5.1             
[100] vctrs_0.6.5             xtable_1.8-4            beeswarm_0.4.0         
[103] evaluate_1.0.4          mvtnorm_1.3-3           cli_3.6.5              
[106] compiler_4.5.1          rlang_1.1.6             crayon_1.5.3           
[109] rstantools_2.4.0        fs_1.6.6                stringi_1.8.7          
[112] optimParallel_1.0-2     Brobdingnag_1.2-9       pacman_0.5.1           
[115] quantreg_6.1            fontquiver_0.2.1        hms_1.1.3              
[118] haven_2.5.5             igraph_2.1.4            RcppParallel_5.1.10    
[121] phangorn_2.12.1         fastmatch_1.1-6         officer_0.6.10         
[124] polynom_1.4-1     
