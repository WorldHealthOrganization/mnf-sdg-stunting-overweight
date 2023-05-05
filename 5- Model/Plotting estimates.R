
marker <- "Stunting"
Sex_val <- "NS"
date_val <- "Mar23"
date <- paste0(Sex_val,"_",date_val)
imp <- FALSE

### Path to output plots
path <- "1- Sample Input Data/Analysis files/"


library(tidyverse)

if(Sex_val == "NS"){
  
  if(marker == "Overweight"){
    data_one <- readRDS(paste0("1- Sample Input Data/Over_data_w_cov_",date_val,".rds"))
  }
  if(marker == "Stunting"){
    data_one <- readRDS(paste0("1- Sample Input Data/Stunt_data_w_cov_",date_val,".rds"))
  }
}

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
  arrange(country, year) #%>% filter(!is.na(Y))


filename <- paste(path,marker," results ",date,".csv",sep="")
if(imp){
  filename <- paste(path,marker," results ",paste0(Sex_val,"_imp_",date_val),".csv",sep="")
}
PP_plot_data <- read.csv(filename)

PP_plot_data <- PP_plot_data %>% 
  rename(
    "Y"="Point.Estimate.Imp","SE_var"="SSE_imp","pred"="Prediction","lower_CI2"="lower_PI","upper_CI2"="upper_PI"
  ) 


Sex_m <- Sex_val
if(Sex_val=="NS"){
  Sex_m <- c("Both", "Male", "Female")
  PP_plot_data$resid[PP_plot_data$Sex=="Both"] <- 0
}

for(j in Sex_m){
  
  P_plot_data <- PP_plot_data %>% filter(Sex == j) 
  
  date_j <- paste0(j,"_",date)
  if(j == Sex_val){date_j <- paste0(Sex_val,"_noNS_",date)}
  
  ################## ONLY INCLUDING COUNTRIES WITH DATA #####################
  
  Countries_w_data <- unique(all_data$country[all_data$country %in% unique(all_data$country[!is.na(all_data$Y)])])
  
  reg_vals <- sort(as.character(unique(all_data$Region)))
  Lim_U <- max(P_plot_data$upper_CI2)
  
  ############################### FIXED SCALES ###############################
  
  pdf(paste(path,marker," estimates ",date_j,"unadj.pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(P_plot_data$country[P_plot_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    plot_data <- plot_data[plot_data$country %in% Countries_w_data,]
    data_two <- plot_data
    
    
    plot_data <- plot_data %>% 
      mutate(Type = case_when(
        (!is.na(adj_pred) & is.na(adj_pred_ns)) ~ "Not Used",
        (ShortSource=="SMART" & resid>3) ~ "SMART Outlier",
        (ShortSource=="SMART") ~ "SMART",
        (resid>3) ~ "Outlier",
        TRUE ~ "Other"
      )) %>% 
      mutate(SE_var = case_when(
        (!is.na(adj_pred) & is.na(adj_pred_ns))  ~ 0,
        TRUE ~ SE_var
      ),
      Type = factor(Type, levels = c("Other","Not Used","SMART", "Outlier","SMART Outlier"))
      )
    
    
    
    ymn <- plot_data$Y-2*(plot_data$SE_var)
    ymn[ymn<0] <- 0
    limits <- aes(y = Y, x = year, ymax = Y+2*(SE_var), ymin=Y - SE_var,color=Type) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + 
      geom_line() + 
      facet_wrap(~country,scales="fixed") +  
      ylim(0,Lim_U) + 
      labs(x="Year", y=paste(marker,"Proportion"),title = paste(marker,"estimates for",reg_nice),size=60) + 
      theme_bw() + 
      theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p+geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI),fill="blue", linetype=2, alpha=0.2) + geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  geom_pointrange(limits,fatten=0.5,size = 0.7) 
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  pdf(paste(path,marker," estimates ",date_j,".pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(P_plot_data$country[P_plot_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    plot_data <- plot_data[plot_data$country %in% Countries_w_data,]
    data_two <- plot_data
    
    plot_data <- plot_data %>% 
      mutate(Type = case_when(
        (!is.na(adj_pred) & is.na(adj_pred_ns)) ~ "Not Used",
        (ShortSource=="SMART" & resid>3) ~ "SMART Outlier",
        (ShortSource=="SMART") ~ "SMART",
        (resid>3) ~ "Outlier",
        TRUE ~ "Other"
      )) %>% 
      mutate(SE_var = case_when(
        (!is.na(adj_pred) & is.na(adj_pred_ns))  ~ 0,
        TRUE ~ SE_var
      ),
      Type = factor(Type, levels = c("Other","Not Used","SMART", "Outlier","SMART Outlier"))
      )
    
    
    ymn <- plot_data$adj_pred-2*(plot_data$SE_var)
    ymn[ymn<0] <- 0
    limits <- aes(y = adj_pred, x = year, ymax = adj_pred+2*(SE_var), ymin=ymn,color=Type) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + 
      ylim(0,Lim_U) + geom_line() + facet_wrap(~country,scales="fixed") 
    
    p <- p +  labs(x="Year", y=paste(marker,"Proportion"),title = paste(marker,"estimates for",reg_nice),size=60) + theme_bw() + theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p+geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI),fill="blue", linetype=2, alpha=0.2) + geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  geom_pointrange(limits,fatten=0.5,size = 0.7)
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  
  
  
  
  ######################### INCLUDING ALL COUNTRIES #########################
  ############################### FREE SCALES ###############################
  
  pdf(paste(path,marker," estimates all countries ",date_j,".pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(P_plot_data$country[P_plot_data$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    data_two <- plot_data
    
    plot_data <- plot_data %>% 
      mutate(Type = case_when(
        (!is.na(adj_pred) & is.na(adj_pred_ns)) ~ "Not Used",
        (ShortSource=="SMART" & resid>3) ~ "SMART Outlier",
        (ShortSource=="SMART") ~ "SMART",
        (resid>3) ~ "Outlier",
        TRUE ~ "Other"
      )) %>% 
      mutate(SE_var = case_when(
        (!is.na(adj_pred) & is.na(adj_pred_ns))  ~ 0,
        TRUE ~ SE_var
      ),
      Type = factor(Type, levels = c("Other","Not Used","SMART", "Outlier","SMART Outlier"))
      )
    
    
    ymn <- plot_data$Y-2*(plot_data$SE_var)
    ymn[ymn<0] <- 0
    limits <- aes(y = Y, x = year, ymax = Y+2*(SE_var), ymin=ymn,color=Type) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred))) + 
      ylim(0,Lim_U) + geom_line() + facet_wrap(~country,scales="fixed")
    
    p <- p +  labs(x="Year", y=paste(marker,"Proportion"),title = paste(marker,"estimates for",reg_nice),size=60) + theme_bw() + theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p+geom_ribbon(data = plot_data, aes(ymin=lower_CI, ymax=upper_CI),fill="blue", linetype=2, alpha=0.2) + geom_line(data=plot_data,aes(x=year,y=(pred)),color="blue")  + geom_line(data=plot_data,aes(x=year,y=lower_CI),color="blue",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI),color="blue",linetype = 2) + geom_line(data=plot_data,aes(x=year,y=lower_CI2),color="green",linetype = 2)+ geom_line(data=plot_data,aes(x=year,y=upper_CI2),color="green",linetype = 2) +  geom_pointrange(limits,fatten=0.5,size = 0.7)
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  ############################## Plotting fixed, pen-fixed and full #############################
  
  Lim_U <- max(P_plot_data$pred)
  pdf(paste(path,marker," fixed penfixed and full estimates ",date_j,".pdf",sep = ""),width = 15, height = 8)
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- unique(data_one$country[data_one$Region==k])
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    data_two <- plot_data
    
    ### Stacking the data to get different types
    
    full_est <-   plot_data %>% dplyr::select(c(country,year,pred)) %>% mutate(Pred_Type = "Full")
    fixed_est <-   plot_data %>% dplyr::select(c(country,year,pred_fixed)) %>% mutate(pred=pred_fixed,Pred_Type = "Fixed Only") %>%  dplyr::select(-pred_fixed)
    fixedpen_est <-   plot_data %>% dplyr::select(c(country,year,pred_fixpen)) %>% mutate(pred=pred_fixpen,Pred_Type = "Fixed w/ Penalized") %>%  dplyr::select(-pred_fixpen)
    
    all_est <- rbind(full_est,fixed_est,fixedpen_est)
    
    p <- ggplot(data=all_est,aes(x=year,y=pred,col=Pred_Type)) + 
      ylim(0,Lim_U) + geom_line() + facet_wrap(~country,scales="fixed")
    
    p <- p +  labs(x="Year", y=paste(marker,"Proportion"),title = paste(marker,"estimates for",reg_nice),size=60) + theme_bw() + theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p + geom_point(data=plot_data, aes(x=year,y=adj_pred),color="black") 
    
    print(p) # This will give a Warning message "Removed XXX rows containing...." ignore this.
  }
  dev.off()
  
  
  
  
}








####### Plots comparing Female Male and Both estimates

  P_plot_data <- as_tibble(PP_plot_data)
  Countries_w_data <- unique(P_plot_data$country[P_plot_data$country %in% unique(P_plot_data$country[!is.na(P_plot_data$Y)])])
  Lim_U <- max(P_plot_data$upper_CI2)
  
  pdf(paste(path,marker," comparison by Sex NS_",date_val,"_2022.pdf",sep = ""),width = 15, height = 8)
  
  for(k in reg_vals){
    reg_nice <- reg_vals[k==reg_vals]
    t_country <- as.character(unique(P_plot_data$country[P_plot_data$Region==k]))
    plot_data <- P_plot_data[P_plot_data$country %in% t_country,]
    plot_data <- plot_data[plot_data$country %in% Countries_w_data,]
    
    
    ymn <- plot_data$adj_pred_ns
    ymn[ymn<0] <- 0
    limits <- aes(y = adj_pred_ns, x = year, ymax = adj_pred_ns, ymin=adj_pred_ns, color=Sex) 
    
    p <- ggplot(data=plot_data, aes(x=year,y=(pred), color=Sex)) + 
      geom_line() + 
      ylim(0,Lim_U) + 
      facet_wrap(~country,scales="fixed")
    p <- p +  labs(x="Year", y=paste(marker,"Proportion"),title = paste(marker,"estimates for",reg_nice),size=60) + theme_bw() + theme(axis.text=element_text(family = "Helvetica", color="#666666",size=10), axis.title = element_text(family = "Helvetica", color="#666666", face="bold", size=22)) 
    p <- p +  geom_pointrange(limits,fatten=0.5,size = 0.7) 
    
    print(p) # This will give a Warning message "Removed XXX rows contiaining...." ignore this.
  }
  
  dev.off()














