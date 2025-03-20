remove(list=ls())

library(this.path)
wd <- dirname(this.path::here())
print(wd)
setwd(wd)

source("Utils/Programs_Feb_2020.R")
library(tidyverse)

marker <- as.character(commandArgs(trailingOnly = TRUE))

month <- "May" #Month Key
year <- "2024" #Year Key
date <- "MultiImp_noSMARTcorr_May24" #for appending filename
marker_f <- paste0(marker,"_noSMARTcorr")

## Path to output folder
if(grepl("St",marker)){
  measure = "Stunting"
}else{
  measure = "Overweight"
}

path = paste0("Data/Analysis files/",measure,"/")


Estimation <- readRDS(paste0(path,"/Mi_files/","Estimation ",marker_f,".rds"))

### Summarizing imputed samples:
boot_vals <- Estimation$boot_vals
Point.Est_mat <- NULL
Point.Est_fixed_mat <- NULL
Point.Est_fixpen_mat <- NULL
Var_mean_mat <- NULL
Var_pred_mat <- NULL
B <- length(boot_vals)
for(j in boot_vals){
  plot_data <- read_rds(paste0(path,"/Mi_files/","/Plot data for ",marker_f," imputation ",j,".rds"))
  cat(dim(plot_data),"\n")
  Point.Est_mat <- cbind(Point.Est_mat,plot_data$pred)
  Point.Est_fixed_mat <- cbind(Point.Est_fixed_mat,plot_data$pred_fixed)
  Point.Est_fixpen_mat <- cbind(Point.Est_fixpen_mat,plot_data$pred_fixpen)
  Var_mean_mat  <- cbind(Var_mean_mat ,plot_data$sigma_T_est^2)
  Var_pred_mat  <- cbind(Var_pred_mat ,plot_data$sigma_Y_est^2)
}

plot_data$pred <- apply(Point.Est_mat,1,mean)
plot_data$pred_fixed <- apply(Point.Est_fixed_mat,1,mean)
plot_data$pred_fixpen <- apply(Point.Est_fixpen_mat,1,mean)
B_Var <- apply(Point.Est_mat,1,var)
W_mean <- apply(Var_mean_mat,1,mean)
W_pred <- apply(Var_pred_mat,1,mean)
plot_data$SE_mean_pred <- sqrt(W_mean + (B+1)/B*B_Var)
plot_data$SE_pred <- sqrt(W_pred + (B+1)/B*B_Var)
plot_data$lower_CI <- plot_data$pred - 1.96*plot_data$SE_mean_pred
plot_data$upper_CI <- plot_data$pred + 1.96*plot_data$SE_mean_pred
plot_data$lower_CI2 <- plot_data$pred - 1.96*plot_data$SE_pred
plot_data$upper_CI2 <- plot_data$pred + 1.96*plot_data$SE_pred

remove(list = c("Point.Est_mat","Point.Est_fixed_mat", 
                "Point.Est_fixpen_mat", "Var_mean_mat",
                "Var_pred_mat"))


data_one <- readRDS(paste0("Data/Merged/",year,"/",marker,"_",month,"_final_multiple_impute.rds")) %>% 
  filter(.imp == j | .imp == 0) %>% 
  dplyr::select(-".imp")

### Outputting the data to "path" folder.
output_function(plot_data, Estimation, data_one, marker, date, path)

