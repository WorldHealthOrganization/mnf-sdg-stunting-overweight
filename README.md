# UNICEF-WHO-World Bank Joint Child Malnutrition Estimates Stunting and Overweight Global Health Estimates
This GitHub repository contains code and sample data for models used by UNICEF, WHO and the World Bank to generate global health estimates for stunting and overweight. 
These estimates are published biennially as the UNICEF-WHO-World Bank Joint Malnutrition Estimates. 

https://data.unicef.org/resources/jme

www.who.int/teams/nutrition-and-food-safety/monitoring-nutritional-status-and-food-safety-and-events/joint-child-malnutrition-estimates

https://datatopics.worldbank.org/child-malnutrition/

## Editions
Details about specific rounds can be found here:
- 2025 Edition: https://github.com/WorldHealthOrganization/mnf-sdg-stunting-overweight-2025
- 2023 Edition: https://github.com/WorldHealthOrganization/mnf-sdg-stunting-overweight-2023

## Input Data Validation and Processing
Details about the input data validation, processing and cleaning can be found here: 
- https://data.unicef.org/resources/jme-standard-methodology/
- https://www.who.int/publications/i/item/9789240100190

## Model
The model used in these analyses is similar to that proposed in McLain et al. (2019) [^1]   
The general statistical model is a penalized longitudinal mixed model with a heterogeneous error term. 
The non-linear longitudinal patterns in the outcomes were captured using penalized cubic B-splines, with country-specific intercepts and random cubic B-splines. 

### Updates to the Methodology since 2023 Round
In prior rounds, a dummy variable indicating whether the survey was based on the SMART methodology was used. This is no longer used in the model

## Procedure for Generating Stunting and Overweight Modelled Estimates

### Contents
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


### Set-up
The programs should be run in the following order:
1. Run the imputation code in folder 2
2. Run programs 1-4 in folder 3
   _Programs 3-4 in folder 3 require the covariate data. _
3. Run analysis programs in folder 4 (Overweight_analy.R and Stunting_analy.R) c
All of the outputs are written to the ```“1- Sample Input Data/Analysis files/”``` folder

## Acknowledgments

### Code
We thank Alexander McLain (@alexmclain) for the development and implementation of the model code.

### Conceptualization of the Model
We thank Alexander McLain, Edward Frongillo, Monika Blössner, Juan Feng, Elaine Borghi, Chika Hayashi, Julia Krasevec, Gretchen Stevens, Mariel Finucane, Leontine Alkema and Simon Cousens, Nicholas Kassebaum for their helpful comments and contributions to this project. 

### Input Data
#### Joint Malnutrition Estimates 2025 Edition
UNICEF: Joel Conkle, Chika Hayashi, Julia Krasevec, Robert Johnston and Vrinda Mehra

WHO: Elaine Borghi, Caroline Dos Santos Costa, Elisa Dominguez, Monica Flores-Urrutia and Giovanna Gatica-Domínguez

World Bank Group: Umar Serajuddin and Emi Suzuki

#### Covariates
Global Burden of Disease Collaborative Network Institute for Health Metrics and Evaluation (IHME)

[^1]: •	McLain A.C., E.A. Frongillo, E. Borghi, and J. Feng (2019). Prediction intervals for heterogeneous penalized longitudinal models with multi-source summary measures: an application to estimating child malnutrition rates. Statistics in Medicine 38:1 1002–1012. GitHub Repo: https://github.com/alexmclain/PHMM
