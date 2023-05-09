# UNICEF-WHO-World Bank Joint Child Malnutrition Estimates

## Acknowledgments 

## Procedure for Generating Stunting and Overweight Modelled Estimates


## Contents
1. Sample Input Data *All the JME input files*
  - List of high-income-countries (**“HIC_list_Nov_2020.csv”**)
  - MCI and SDI covariate information from IHME (**“GBD 2022 MCI and SDI.csv”**)
  - Various economic and health related covariate data used for imputation
  (**“API_NY.GDP.MKTP.CD_DS2_en_csv_v2_4701247.csv**” and “**WPP2022_Demographic_Indicators_Medium.csv**”)

2. Preparing Covariates
  -  R code to perform single imputation of the covariate data
 **(“Imputation of IHME covariate data_single_impute.R”)**

3.  Preparing Primary Data
  - R code to calculate (if possible) and impute missing SE information (“**1 - SE_clean_and_impute.R”**)
  - R code to cross-walk surveys with partial age ranges **(“2 - Age_range_analysis.R”)**.
  - R code to merge covariate and survey data **(“3 - Merging of covariate and survey data_single_impute.R”)**
  - R code to cross-walk surveys with partial sex coverage, and remove redundant sex combinations **(i.e., remove “both” when “male” and “female” are included) (“4 - Sex_cross_walk_single_impute.R”)**

4. Model
  - R functions needed for various modeling stages **(“Programs_Feb_2020.R” and “Programs_Cleaning_SE.R”)**
  - Programs for analyzing the data **(“Overweight_analy.R” and “Stunting_analy.R”)**
  - Programs for plotting the estimates **(“Plotting_estimates.R”)**
