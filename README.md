# UNICEF-WHO-World Bank Joint Child Malnutrition Estimates- Stunting and Overweight Global Health Estimates
This GitHub repository contains code and sample data for models used by [UNICEF](https://data.unicef.org/resources/jme), [WHO](https://www.who.int/teams/nutrition-and-food-safety/monitoring-nutritional-status-and-food-safety-and-events/joint-child-malnutrition-estimates) and the [World Bank](https://datatopics.worldbank.org/child-malnutrition/) to generate global health estimates for stunting and overweight. 

These estimates are published biennially as the UNICEF-WHO-World Bank Joint Malnutrition Estimates. 

## Editions
Details about specific rounds can be found here:
- 2025 Edition: https://github.com/WorldHealthOrganization/mnf-sdg-stunting-overweight-2025
- 2023 Edition: https://github.com/WorldHealthOrganization/mnf-sdg-stunting-overweight-2023

## Input Data Validation and Processing
Details about the input data validation, processing and cleaning can be found here: 
- https://data.unicef.org/resources/jme-standard-methodology/
- https://www.who.int/publications/i/item/9789240100190

## Model
The model used in these analyses is based on McLain et al. (2019) [^1]  and and application of methods published in Saraswati et al. (2022) [^2], and a corresponding editorial by Finaret (2022) [^3].

The general statistical model is a penalized longitudinal mixed model with a heterogeneous error term. 
The non-linear longitudinal patterns in the outcomes were captured using penalized cubic B-splines, with country-specific intercepts and random cubic B-splines. 

### Updates to the Methodology since 2023 Round
In prior rounds, a dummy variable indicating whether the survey was based on the SMART methodology was used. This is no longer used in the model

## Procedure for Generating Stunting and Overweight Modelled Estimates

### Contents

+ Data 
  - Economic Covariate Data from the World Bank `\WB\API_NY.GDP.MKTP.CD_DS2_en_csv_v2_4701247.csv`
  - Demographic Covariate Data from the United Nations Population Division (UNPD) `\WPP\WPP2022_Demographic_Indicators_Medium.csv`
  - MCI and SDI covariate information from IHME `\IHME_covs\GBD 2022 MCI and SDI.csv`
  - Regional Grouping `\Country\Crosswalk_Jan_2025.xlsx`
  - Sample Input Data `\JME\2024\Raw`

+ Utils
  - Contains various programs and functions required to run the analyses.

+ 0 - Setup
  - R code to install required packages `1- Install Required Packages.R`

+ 1 - Preparing Covariates
  - R code to perform multiple imputation of the covariate data   `Imputation of IHME covariate data_multiple_impute.R`

+ 2 - Preparing Primary Data
  - R code to calculate (if possible) and impute missing SE information `1 - SE_clean_and_impute.R`
  - R code to integrate data sources that were originally collected using different (non-standard) age ranges `2 - Age_range_analysis.R`
  - R code to merge covariate and survey data `3 - Merging Survey_Covariate_Data.R`
  - R code to generates approximate sex-specific prevalence rates for specific years and countries (solely for use in visual figures) `4 - Sex_cross_walk.R`

+ 3 - Model
  - R code for performing a multiple imputed analysis of the data for stunting or overweight `1 - Analysis_MI.R`

+ 4 - Results
  - R code for pooling imputed estimates `1 - Summarizing_imputed_results.R`

### Set-up
The programs should be run in the following order:
1. Set-Up: Run `1- Install Required Packages.R`
2. Preparing Covariates: Run `1- Imputation of IHME covariate data_multiple_impute.R`
3. Preparing Primary Data: Run `1 - SE_clean_and_impute.R`
4. Preparing Primary Data: Run `2 - Age_range_analysis.R`
5. Preparing Primary Data: Run `3 - Merging Survey_Covariate_Data.R`
6. Preparing Primary Data: Run `4 - Sex_cross_walk.R`
7. Model: Run `1 - Analysis_MI.R`
8. Results: Run `1 - Summarizing_imputed_results.R`

All the programs in the folders `2- Preparing Primary Data`, `3- Model`, and `4- Results` need to be run separately for both Stunting and Overweight.

Ensure that you read each file carefully to verify that the files are named and stored correctly.

You can do this by setting marker to "Stunting" or "Overweight" in each script
```
marker<-"Stunting" 
marker<-"Overweight"
```

### Scaling Outputs
All outputs are further rescaled to align with 2024 UNPD World Population Prospects.

## Acknowledgments

### Code
We thank Alexander McLain (@alexmclain) for the development and implementation of the model code.

### Conceptualization of the Model
We thank Alexander McLain, Edward Frongillo, Monika Blössner, Juan Feng, Elaine Borghi, Chika Hayashi, Julia Krasevec, Gretchen Stevens, Mariel Finucane, Leontine Alkema and Simon Cousens, Nicholas Kassebaum for their helpful comments and contributions to this project. 

### Update to include Sex Specific Estimates
We thank Alexander McLain, Edward Frongillo, Elaine Borghi, Elisa Dominguez, Caroline Dos Santos Costa, Monica Flores-Urrutia, Giovanna Gatica-Domínguez, Chika Hayashi, Robert Johnston, Julia Krasevec, Richard Kumapley and Emi Suzuki for their efforts and contributions to update the methodology to include sex-specific estimates.

### Input Data
#### Joint Malnutrition Estimates 2025 Edition
UNICEF: Joel Conkle, Chika Hayashi, Julia Krasevec, Robert Johnston and Vrinda Mehra

WHO: Elaine Borghi, Caroline Dos Santos Costa, Elisa Dominguez, Monica Flores-Urrutia and Giovanna Gatica-Domínguez

World Bank Group: Umar Serajuddin and Emi Suzuki

#### Covariates
Global Burden of Disease Collaborative Network Institute for Health Metrics and Evaluation (IHME)

United Nations Population Division (UNPD)

World Bank Group (WB)

[^1]: •	McLain A.C., E.A. Frongillo, E. Borghi, and J. Feng (2019). Prediction intervals for heterogeneous penalized longitudinal models with multi-source summary measures: an application to estimating child malnutrition rates. Statistics in Medicine 38:1 1002–1012. <a href="https://doi.org/10.1093/jn/nxac072">doi.org/10.1093/jn/nxac072</a> GitHub Repo: https://github.com/alexmclain/PHMM
[^2]: •	Saraswati, Chitra M, Elaine Borghi, João JR da Silva Breda, Monica C Flores-Urrutia, Julianne Williams, Chika Hayashi, Edward A Frongillo, and Alexander C McLain. 2022. “Estimating Childhood Stunting and Overweight Trends in the European Region from Sparse Longitudinal Data.” *The Journal of Nutrition* 152 (7): 1773–82 <a href="https://doi.org/10.1002/sim.8024">doi.org/10.1002/sim.8024</a>
[^3]: •	Finaret, Amelia B. 2022. “Advancing Nutritional Epidemiology by Linking Datasets and Addressing Data Quality.” *The Journal of Nutrition* 152 (7): 1595–96. <a  ref="https://doi.org/10.1093/jn/nxac092">doi.org/10.1093/jn/nxac092</a>


