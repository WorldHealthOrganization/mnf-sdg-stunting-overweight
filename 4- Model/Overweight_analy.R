remove(list=ls())

wd <- "~path to root directory"
setwd(wd)

source("4- Model/Programs_Feb_2020.R")
library(tidyverse)

marker <- "Overweight"
## Date needs to match what was used 
## in "5 - Sex_cross_walk_single_impute.R"
date <- paste0("NS_Mar23")
data_one <- readRDS("1- Sample Input Data/Over_data_w_cov_Mar23.rds") 

## Path to output folder
path = "1- Sample Input Data/Analysis files/"


all_data <- data_one %>% 
  mutate(
    Sex = case_when(
      Sex == "Both" ~ 0,
      Sex == "Female" ~ 1,
      Sex == "Male" ~ -1
    ),
    Region = factor(Region), 
    "SMART" = case_when(
      ShortSource=="SMART" ~ 1,
      TRUE ~ 0
    ),
    "Surveillance" = case_when(
      ShortSource=="Surveillance" ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  select(c("country","year", "Point.Estimate.NS", "Point.Estimate.Imp","SE_val", "ShortSource",
           "Region","SEV", "Sex", "SMART", "Surveillance", "MCI_5_yr", "SDI")) %>% 
  rename("Y" = "Point.Estimate.NS", 
         "Y_all" = "Point.Estimate.Imp",
         "SE_var" = "SE_val") %>% 
  arrange(country, year) 

View(all_data)
length(unique(all_data$country))
length(c(all_data$Y[!is.na(all_data$Y)]))

all_data$P_Region <- all_data$Region
levels(all_data$P_Region)
referent_lev <- 2
all_data$P_Region <- relevel(all_data$P_Region,ref = referent_lev)
levels(all_data$P_Region)


## Transforming outcome ##
data_w_out <- all_data %>% 
  select(c("country", "year", "Y", "SE_var")) %>% 
  mutate(
    SE_pred = 0, 
    SE_var = SE_var*((1/Y) + 1/(1-Y)), 
    Y = log(Y/(1-Y))
  )


Pcov_data <- model.matrix(~P_Region,data=all_data,contrasts.arg = list(P_Region=diag(nlevels(all_data$P_Region))))[,-1]

TRANS=FALSE
##Number of penalized and random splines.  Equally spaced throughout the B.knots.
DF_P <- seq(1993,2022,5)
DF_R <- NULL 
B.knots <- range(data_w_out$year[!is.na(data_w_out$Y)])
B.knots[1] <- B.knots[1] - 1
B.knots[2] <- B.knots[2] + 10
q.order = 2
plots=FALSE
slope = TRUE

#Specifying the covariance matrix of the random effects.  
#   - cov_mat="CS" (default) -> compound symmetric covariance matrix
#   - cov_mat="UN"  -> unstructured symmetric covariance matrix
#   - cov_mat="VC"  -> Diagonal only covariance matrix
cov_mat <- "VC"


cov_data <- as.matrix(data.frame(model.matrix(
  ~ SMART + Sex + P_Region + MCI_5_yr, 
  data = all_data)[,-1]
))
zero_covs <- colnames(cov_data)[1]


##################### Covariate analysis with multiple penalized functions #################################
Estimation <- cmnpe(data_w_out, DF_P, DF_R, B.knots,q.order , cov_data=cov_data, Pcov_data = Pcov_data, cov_mat = cov_mat,plots=plots,TRANS=TRANS,zero_covs=zero_covs,slope = slope)

summary(Estimation$result$model)

### Outputting the data to "path" folder.
output_function(Estimation, all_data, data_one, marker, date, path)

