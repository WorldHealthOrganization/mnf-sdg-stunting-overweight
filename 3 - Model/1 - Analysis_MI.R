remove(list=ls())

library(this.path)
wd <- dirname(this.path::here())
print(wd)
setwd(wd)
source("Utils/Programs_Feb_2020.R")
library(tidyverse)

marker <- as.character(commandArgs(trailingOnly = TRUE))
year <- "2024"
month <- "Dec" #for finding the data
marker_f <- paste0(marker," December Model New Regions")

#### 1. Setting Model Parameters #### 
if(grepl("St",marker)){
  measure = "Stunting"
  model_formula <- ~Sex + Region + MCI_5_yr+ I(MCI_5_yr^2) + SDI + 
    MCI_5_yr:Sex+ I(MCI_5_yr^2):Sex + Region:Sex
  
  ### Pcov_data: penalized covariate data (optional).  Usually a regional groupings.
  ##    The columns should consist of 1's if the row is in the group and zero otherwise.
  ##    For example, if there are three regions Pcov_data should have three rows, where
  ##    row i column j has a 1 is that observations is in region j.  It's important that 
  ##    every group has a row.
  Pcov_data <- NULL
  
  ##Penalized and random splines. 
  DF_P <- seq(1993,2023,1)
  
  #Specifying the covariance matrix of the random effects.  
  #   - cov_mat="CS" (default) -> compound symmetric covariance matrix
  #   - cov_mat="UN"  -> unstructured symmetric covariance matrix
  #   - cov_mat="VC"  -> Diagonal only covariance matrix
  cov_mat <- "VC"
  slope <- TRUE
  
}else{
  measure = "Overweight"
  model_formula <- ~ Sex + P_Region + MCI_5_yr
  
  ##Penalized and random splines. 
  DF_P <- 2005
  
  #Specifying the covariance matrix of the random effects.  
  #   - cov_mat="CS" (default) -> compound symmetric covariance matrix
  #   - cov_mat="UN"  -> unstructured symmetric covariance matrix
  #   - cov_mat="VC"  -> Diagonal only covariance matrix
  cov_mat <- "VC"
  slope <- TRUE
}
## Path to output folder
path = paste0("Data/Analysis files/",measure,"/")

TRANS <- FALSE
q.order <- 2
plots <- FALSE
B <- 20
boot_vals <- 1:B
#### 2. Running and Outputting the Models #### 
for(j in 1:B){
  
  all_data <- readRDS(paste0("Data/Merged/",year,"/",marker,"_",
                             month,"_final_multiple_impute.rds"))  %>% 
    filter(.imp == j | .imp == 0)  %>% 
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
    dplyr::select(c(".imp","country","year", "Point.Estimate.NS", "Point.Estimate.Imp","SE_val", "ShortSource",
                    "Region","SEV", "Sex", "SMART", "Surveillance", "MCI_5_yr", "SDI")) %>% 
    rename("Y" = "Point.Estimate.NS", 
           "Y_all" = "Point.Estimate.Imp",
           "SE_var" = "SE_val") %>% 
    arrange(country, year) #%>% filter(!is.na(Y))
  
  cat(j,length(unique(all_data$country)), 
      length(c(all_data$Y)),
      length(c(all_data$Y[!is.na(all_data$Y)])), 
      table(all_data$Region),"\n") 
  
  if(measure == "Overweight"){
    all_data$P_Region <- all_data$Region
    levels(all_data$P_Region)
    referent_lev <- 2
    all_data$P_Region <- relevel(all_data$P_Region,ref = referent_lev)
    levels(all_data$P_Region)
    
    Pcov_data <- model.matrix(~P_Region,data=all_data,contrasts.arg = list(P_Region=diag(nlevels(all_data$P_Region))))[,-1]
  }
  
  ## Transforming outcome ##
  data_w_out <- all_data %>% 
    dplyr::select(c("country", "year", "Y", "SE_var")) %>% 
    mutate(
      SE_pred = 0, 
      SE_var = SE_var*((1/Y) + 1/(1-Y)), 
      Y = log(Y/(1-Y))
    )
  
  if(measure == "Stunting"){
    DF_R <- quantile(data_w_out$year,probs = c(0.5))
    B.knots <- c(min(data_w_out$year[!is.na(data_w_out$Y)])-1, 
                 max(data_w_out$year[!is.na(data_w_out$Y)])+2)
  }else{
    DF_R <- NULL  
    B.knots <- range(data_w_out$year[!is.na(data_w_out$Y)])
    B.knots[1] <- B.knots[1] - 1
    B.knots[2] <- B.knots[2] + 7
  }
  
  cov_data <- as.matrix(data.frame(model.matrix(
    model_formula,
    data=all_data)[,-1]
  ))
  zero_covs <- NULL
  
  remove(all_data)
  
  ##################### Covariate analysis with multiple penalized functions #################################
  t1 <- try(Estimation <- cmnpe(data_w_out, DF_P, DF_R, B.knots, q.order, 
                                cov_data = cov_data, Pcov_data = Pcov_data, 
                                cov_mat = cov_mat, plots = plots, TRANS=TRANS,
                                zero_covs = zero_covs, slope = slope))
  
  if(is.null(attr(t1,"class"))){
    cat(j,2*Estimation$df - 2*c(summary(Estimation$result$model)$logLik) + 
          summary(Estimation$result$model)$AIC + 
          2*c(summary(Estimation$result$model)$logLik),"\n")
    
    t_plot_data <- Estimation$plot_data
    saveRDS(t_plot_data, 
            file = paste0(path,"MI_files/Plot data for ",marker_f," imputation ",j,".rds"))
    
    outfile <- list(
      gamma = c(Estimation$result$model$modelStruct$varStruct[1]),
      sigma2 = Estimation$result$model$sigma^2,
      smart_cov = Estimation$result$model$coefficients$fixed[which(colnames(cov_data) %in% zero_covs) +2] ,
      zero_covs = zero_covs, 
      boot_vals = boot_vals[boot_vals<=j])
    
    saveRDS(outfile, file = paste0(path,"MI_files/Estimation ",marker_f,".rds"))
    remove(Estimation)
    
  }else{
    boot_vals <- boot_vals[boot_vals != j]
    cat(j,"Model Error.\n")
  }
}



