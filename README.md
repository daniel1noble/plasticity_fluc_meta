# Temperature variability does not impact the capacity for phenotypic plasticity in ectotherms – a meta-analysis

This repository contains data and R code for the associated manuscript.  

Clayton W. Stocker, Stephanie M. Bamford, Miki Jahn, Geoffrey P. F. Mazué, Amanda K. Pettersen, Daniel Ritchie, Alexander Rubin, Daniel W.A. Noble, Frank Seebacher. (2025) Temperature variability does not impact the capacity for phenotypic plasticity in ectotherms – a meta-analysis.

Please email clayton.stocker@sydney.edu.au for any inquiries. 

## Analyses

All code for analyses reproducing the results and figures is found in `Analysis_Final.R` script. The code also generates supplementary tables and results.

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