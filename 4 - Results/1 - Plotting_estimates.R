remove(list=ls())

library(tidyverse)
library(this.path)
wd <- dirname(this.path::here())
print(wd)
setwd(wd)


marker <- as.character(commandArgs(trailingOnly = TRUE))

month <- "May" #Month Key
year <- "2024" #Year Key
date <- "MultiImp_noSMARTcorr_May24" #for appending filename


if(grepl("St",marker)){
  measure = "Stunting"
}else{
  measure = "Overweight"
}

path_data = paste0("Data/Analysis files/",measure,"/")
path_fig = paste0("Figures/",measure,"/",month,year,"/")

all_data <- readRDS(paste0("Data/Merged/",year,"/",marker,"_",
                           month,"_final_multiple_impute.rds")) %>% 
  filter(.imp == 1 | .imp == 0) %>% 
  mutate(
    Sex = case_when(
      Sex == "Overall" ~ 0,
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


filename <- paste(path_data,marker," results ",date,".csv",sep="")

PP_plot_data <- read.csv(filename)

PP_plot_data <- PP_plot_data %>% 
  rename(
    "Y"="Point.Estimate.Imp","SE_var"="SSE_imp","pred"="Prediction","lower_CI2"="lower_PI","upper_CI2"="upper_PI"
  ) 

Sex_m <- c("Overall", "Male", "Female")

for(j in Sex_m){
  
  ### Limiting to appropriate sex and creating "Type" variable
  P_plot_data <- PP_plot_data %>% 
    filter(Sex == j) %>% 
      mutate(Type = case_when(
        (resid>3) ~ "Outlier",
        TRUE ~ "Survey"
      )) %>% 
      mutate(
      Type = factor(Type, levels = c("Survey","Outlier"))
      )
  
  date_j <- paste0(j,"_",date)
  
  ################## ONLY INCLUDING COUNTRIES WITH DATA #####################
  
  Countries_w_data <- unique(all_data$country[all_data$country %in% unique(all_data$country[!is.na(all_data$Y)])])
  
  reg_vals <- sort(as.character(unique(all_data$Region)))
  Lim_U <- max(P_plot_data$upper_CI2)
  
  ############################### FIXED SCALES ###############################
  
  pdf(paste(path_fig,marker," estimates ",date_j,"unadj.pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(P_plot_data$country[P_plot_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    plot_data <- plot_data[plot_data$country %in% Countries_w_data,]
    data_two <- plot_data
    
    ymn <- plot_data$Y-2*(plot_data$SE_var)
    ymn[ymn<0] <- 0
    limits <- aes(y = Y, x = year, ymax = Y+2*(SE_var), ymin=Y - SE_var,color=Type) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + 
      geom_line() + 
      facet_wrap(~country,scales="fixed") +  
      ylim(0,Lim_U) + 
      labs(x="Year", y=paste(measure,"Proportion"), 
           title = paste(marker,"estimates for",reg_nice),size=60) + 
      theme_bw() + 
      theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), 
            axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p + geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI), 
                         fill="blue", linetype=2, alpha=0.2) + 
      geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + 
      geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2)+ 
      geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2)+ 
      geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  
      geom_pointrange(limits,fatten=0.5,size = 0.7) 
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  if(!is.null(P_plot_data$adj_pred)){
  
  pdf(paste(path_fig,marker," estimates ",date_j,".pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(P_plot_data$country[P_plot_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    plot_data <- plot_data[plot_data$country %in% Countries_w_data,]
    data_two <- plot_data
    
    ymn <- plot_data$adj_pred-2*(plot_data$SE_var)
    ymn[ymn<0] <- 0
    limits <- aes(y = adj_pred, x = year, ymax = adj_pred+2*(SE_var), ymin=ymn,color=Type) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + 
      ylim(0,Lim_U) + geom_line() + facet_wrap(~country,scales="fixed") 
    
    p <- p + labs(x="Year", y=paste(measure,"Proportion"),title = paste(marker,"estimates for",reg_nice),size=60) + 
      theme_bw() + 
      theme(axis.text = element_text(family = "Helvetica", color="#666666",size=10), 
            axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p + geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI),fill="blue", linetype=2, alpha=0.2) + 
      geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + 
      geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  
      geom_pointrange(limits,fatten=0.5,size = 0.7)
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  }
  
  
  
  
  
  ######################### INCLUDING ALL COUNTRIES #########################
  ############################### FREE SCALES ###############################
  
  pdf(paste(path_fig,marker," estimates all countries ",date_j,".pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(P_plot_data$country[P_plot_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    data_two <- plot_data
    
    ymn <- plot_data$Y-2*(plot_data$SE_var)
    ymn[ymn<0] <- 0
    limits <- aes(y = Y, x = year, ymax = Y+2*(SE_var), ymin=ymn,color=Type) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + 
      ylim(0,Lim_U) + geom_line() + facet_wrap(~country,scales="fixed")
    
    p <- p +  labs(x="Year", y=paste(measure,"Proportion"), 
                   title = paste(marker,"estimates for",reg_nice),size=60) + 
      theme_bw() + 
      theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), 
            axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p + geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI),fill="blue", linetype=2, alpha=0.2) + 
      geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + 
      geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2) + 
      geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  
      geom_pointrange(limits,fatten=0.5,size = 0.7)
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  ############################## Plotting fixed, pen-fixed and full #############################
  
  Lim_U <- max(P_plot_data$pred)
  pdf(paste(path_fig,marker," fixed penfixed and full estimates ",date_j,".pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(all_data$country[all_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    data_two <- plot_data
    
    ### Stacking the data to get different types
    
    full_est <-   plot_data %>% dplyr::select(c(country,year,pred)) %>% mutate(Pred_Type = "Full")
    fixed_est <-   plot_data %>% dplyr::select(c(country,year,pred_fixed)) %>% mutate(pred=pred_fixed,Pred_Type = "Fixed Only") %>%  dplyr::select(-pred_fixed)
    fixedpen_est <-   plot_data %>% dplyr::select(c(country,year,pred_fixpen)) %>% mutate(pred=pred_fixpen,Pred_Type = "Fixed w/ Penalized") %>%  dplyr::select(-pred_fixpen)
    
    all_est <- rbind(full_est,fixed_est,fixedpen_est)
    
    p <- ggplot(data=all_est,aes(x=year,y=pred,col=Pred_Type)) + 
      ylim(0,Lim_U) + geom_line() + facet_wrap(~country,scales="fixed")
    
    p <- p +  labs(x="Year", y=paste(measure,"Proportion"),
                   title = paste(marker,"estimates for",reg_nice),size=60) + 
      theme_bw() + 
      theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), 
            axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p + geom_point(data=plot_data, aes(x=year,y=Y),color="black") 
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  
  
}








####### Plots comparing Female, Male, and Overall estimates

  P_plot_data <- as_tibble(PP_plot_data) 
  Countries_w_data <- unique(P_plot_data$country[P_plot_data$country %in% unique(P_plot_data$country[!is.na(P_plot_data$Y)])])
  Lim_U <- max(P_plot_data$upper_CI2)
  
  pdf(paste0(path_fig,marker," comparison by Sex NS_",date,".pdf"),width = 15, height = 8)
  
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- as.character(unique(P_plot_data$country[P_plot_data$Region==k]))
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    plot_data <- plot_data[plot_data$country %in% Countries_w_data,]
    
    
    ymn <- plot_data$Y
    ymn[ymn<0] <- 0
    limits <- aes(y = plot_data$Y, x = year, ymax = plot_data$Y, ymin=plot_data$Y, color=Sex) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred), color=Sex)) + 
      geom_line() + 
      ylim(0,Lim_U) + 
      facet_wrap(~country,scales="fixed")
    p <- p +  labs(x="Year", y=paste(measure,"Proportion"), 
                   title = paste(marker,"estimates for",reg_nice),size=60) + 
      theme_bw() + 
      theme(axis.text = element_text(family = "Helvetica", color="#666666",size=10), 
            axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p +  geom_pointrange(limits,fatten=0.5,size = 0.7) 
    
    print(p) # This will give a Warning message "Removed XXX rows contiaining...." ignore this.
  }
  
  dev.off()














