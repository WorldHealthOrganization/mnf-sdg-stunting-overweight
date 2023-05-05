
wd <- "~path to root directory"
setwd(wd)

source("5- Model/Programs_Cleaning_SE.R")
library(readxl)
library(nlme)
library(tidyverse)
library(dplyr)


marker <- "Stunting"
if(marker == "Overweight"){
  marker = "Over"
}else{
  if(marker == "Stunting"){
    marker = "Stunt"
  }else{
    stop("Marker doesn't match options.")
  }
}
#
month = "Mar23" #for appending filename
surv_data <- read_rds(paste0("1- Sample Input Data/Stunt_data_w_cov_",month,"_reg_all.rds"))

##### Creating a dataset where both is missing when male and female are observed.

surv_data <- surv_data %>% 
  mutate(UNICEFSurveyID = case_when(
    UNICEFSurveyID == 0 ~ NA_real_,
    TRUE ~ UNICEFSurveyID)
  ) %>% 
  group_by(ISO.code, year) %>% 
  fill(UNICEFSurveyID,
       .direction = "updown") 
  

surv_data <- surv_data %>% 
  group_by(UNICEFSurveyID) %>% 
  mutate(check = all(!is.na(Point.Estimate)) & n()==3) %>% 
  ungroup() %>% 
  mutate(Point.Estimate.NS = case_when(
    (check==TRUE & Sex == "Both") ~ NA_real_,
    TRUE ~ Point.Estimate
  ))

#surv_data %>%  filter(country == "Bangladesh") %>%View()
#surv_data %>%  filter(country == "Uganda") %>%View()
########################### Now to the Sex cross walk

surv_wide_one <- surv_data %>% 
  dplyr::select(
    c(ISO.code, country,  year, Sex, UNICEFSurveyID,
      Point.Estimate, SE_val, SEV, All_africa_HIW)
  ) %>% 
  mutate(
    year_ID = case_when(
      is.na(UNICEFSurveyID) ~ year,
      !is.na(UNICEFSurveyID) ~ year+UNICEFSurveyID
    ), 
    UNICEFSurveyID = case_when(
      is.na(UNICEFSurveyID) ~ 0,
      !is.na(UNICEFSurveyID) ~ UNICEFSurveyID
    )) %>% 
  group_by(ISO.code, country, year, UNICEFSurveyID, All_africa_HIW, year_ID, Sex) %>%
  filter(row_number() == 1) %>%  
  pivot_wider(
    names_from = Sex, 
    values_from = c(Point.Estimate, SE_val, SEV)
  ) %>% 
  mutate(
    SE_val_Female = case_when(
      is.na(SE_val_Female) ~ SE_val_Both, 
      TRUE ~ SE_val_Female
    ), 
    SE_val_Male = case_when(
      is.na(SE_val_Male) ~ SE_val_Both, 
      TRUE ~ SE_val_Male
    )
  ) %>% 
  filter(!is.na(Point.Estimate_Both))

surv_wide_B <- surv_data %>% 
  filter(Sex == "Both") %>% 
  dplyr::select(
    c(.imp, ISO.code, country,  year, UNICEFSurveyID, ShortSource)
  ) %>% 
  right_join(surv_wide_one) 


################# Look at variable relationships #################
if(FALSE){
  ggplot(data=surv_wide_B, aes(x=Point.Estimate_Both,y=Point.Estimate_Female)) + 
    geom_point() + 
    theme(text=element_text(family="Garamond", size=14))
  
  ggplot(data=surv_wide_B, aes(x=log(Point.Estimate_Both),y=log(Point.Estimate_Female))) + 
    geom_point() + 
    theme(text=element_text(family="Garamond", size=14))
  
  ggplot(data=surv_wide_B, aes(x=sqrt(Point.Estimate_Both),y=sqrt(Point.Estimate_Female))) + 
    geom_point() + 
    theme(text=element_text(family="Garamond", size=14))
  
  
  ggplot(data=surv_wide_B, aes(x=Point.Estimate_Both,y=Point.Estimate_Male)) + 
    geom_point() + 
    theme(text=element_text(family="Garamond", size=14)) 
  
  ggplot(data=surv_wide_B, aes(x=log(Point.Estimate_Both),y=log(Point.Estimate_Male))) + 
    geom_point() + 
    theme(text=element_text(family="Garamond", size=14)) 
  ggplot(data=surv_wide_B, aes(x=sqrt(Point.Estimate_Both),y=sqrt(Point.Estimate_Male))) + 
    geom_point() + 
    theme(text=element_text(family="Garamond", size=14)) 
  
  ggplot(data=surv_wide_B, aes(y = sqrt(Point.Estimate_Female) - sqrt(Point.Estimate_Both),
                               x = sqrt(Point.Estimate_Both))) + 
    geom_point() +
    geom_smooth(method = "lm") +
    theme(text=element_text(family="Garamond", size=14))
  
  
  ggplot(data=surv_wide_B, aes(y = sqrt(Point.Estimate_Male) - sqrt(Point.Estimate_Both),
                               x = sqrt(Point.Estimate_Both))) + 
    geom_point() +
    geom_smooth(method = "lm") +
    theme(text=element_text(family="Garamond", size=14))
}


### Predicting for Females
surv_wide <-  surv_wide_B %>%  
  filter(.imp == 1) %>% 
  mutate(Source = relevel(factor(All_africa_HIW), ref = "South Asia")) %>% 
  dplyr::select(-c(".imp"))

surv_wide$SEV_Both[is.na(surv_wide$SEV_Both)] <- mean(surv_wide$SEV_Both[!is.na(surv_wide$SEV_Both)])
surv_wide$SEV_Female[is.na(surv_wide$SEV_Female)] <- mean(surv_wide$SEV_Female[!is.na(surv_wide$SEV_Female)])
surv_wide$SEV_Male[is.na(surv_wide$SEV_Male)] <- mean(surv_wide$SEV_Male[!is.na(surv_wide$SEV_Male)])

if(marker=="Over"){
  formula <- sqrt(Point.Estimate_Female) ~ 
    sqrt(Point.Estimate_Both)
  formula_x <- ~ 
    sqrt(Point.Estimate_Both) 
  cnt_list = list(maxIter = 5000, opt = "optim")
}
if(marker=="Stunt"){
  formula <- sqrt(Point.Estimate_Female) ~ 
    sqrt(Point.Estimate_Both) 
  formula_x <- ~ 
    sqrt(Point.Estimate_Both) 
  cnt_list = list(maxIter = 5000)
}



#### Fitting the model

fitted_model <- lme(formula,
                    random=~1|country,
                    weights = varIdent(form=~1|Source),
                    data = surv_wide, na.action = na.omit, 
                    control = cnt_list)

all_data <- surv_wide
data_frame <- surv_wide %>% 
  mutate(Y = sqrt(Point.Estimate_Female)) %>% 
  filter(!is.na(Y))

Xi.matrix <- cbind(model.matrix(formula_x ,data=data_frame))
SE.matrix <- cbind(model.matrix(~Source,data=data_frame)[,-1])
N <- length(unique(data_frame$country))
X_orig <- Xi.matrix
Z_orig <- make_Z_noPcov(data_frame,data_frame)

# Estimate EBLUPs, predicted means, and variances of the predicted means.
pred_mat <- pred_fun_noPcov(X=X_orig,Z=Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=SE.matrix,full_b = 1)
S <- pred_mat$S
b <- pred_mat$b


# Extracting only those countries that were used in the original analysis
cov_data_used_only <- all_data[all_data$country %in% data_frame$country,]
# Ordering the data so the random effects align.
div=c(cov_data_used_only$country)
cov_data_used_only <- cov_data_used_only[order(div,cov_data_used_only$year),]
# Making the X and Z matrices for the new data.

X <- model.matrix(formula_x,data=cov_data_used_only)
SE<- cbind(model.matrix(~Source,data=cov_data_used_only)[,-1])
Z <- make_Z_noPcov(cov_data_used_only,data_frame)

# Predicting the mean and variance for each new observation.
pred_mat_new <- pred_fun_noPcov(X,Z,Y=0,fitted_model,N,SE_var=SE,full_b = b, S = S)

pred_all  <- pred_mat_new$mu_T^2
sigma_Y_all <- diag(pred_mat_new$sigma_Y)*4*pred_all

pred_new_data <- data.frame(cov_data_used_only,pred_all, sigma_Y_all)
all_data <- left_join(all_data, pred_new_data)
plot(all_data$Point.Estimate_Both, all_data$pred_all)
lines(c(0,100),c(0,100))
plot(sqrt(all_data$Point.Estimate_Both), sqrt(all_data$pred_all) )
lines(c(0,100),c(0,100))


# Extracting countries that were not used in the original analysis
cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]

pred_all_noC <- 0
sigma_Y_all_noC <- 0
if(nrow(cov_data_not_used)>0){
  # Predicting for them
  Xi.full <- cbind(model.matrix(formula_x ,data=all_data))
  SE.full <- cbind(model.matrix(~Source,data=all_data)[,-1])
  X <- Xi.full[not_in(all_data$country,data_frame$country),]
  SE<- SE.full[not_in(all_data$country,data_frame$country),]
  
  N <- length(unique(cov_data_not_used$country))
  pred_mat_new_noC <- pred_fun_noPcov_noC(X,Y,fitted_model,N,SE_var=SE,full_b=b,S=S)
  pred_all_noC  <- pred_mat_new_noC$mu_T^2
  sigma_Y_all_noC <- diag(pred_mat_new_noC$sigma_Y)*4*pred_all_noC
  
  
  # Predicting the standard errors.
  pred_new_data_noC <- data.frame(cov_data_not_used,pred_all_noC, sigma_Y_all_noC)
  all_data <- left_join(all_data, pred_new_data_noC)
  plot(all_data$Point.Estimate_Both, all_data$pred_all_noC)
  lines(c(0,100),c(0,100))
}

### Final data for Female
surv_wide_Female <- surv_wide_B %>% 
  left_join( all_data) %>%   
  mutate(
    Sex = "Female",
    pred_all = case_when(
      is.na(pred_all) ~ pred_all_noC, 
      TRUE ~ pred_all
    ),
    sigma_Y_all = case_when(
      is.na(sigma_Y_all) ~ sigma_Y_all_noC, 
      TRUE ~ sigma_Y_all
    )
  ) %>% 
  dplyr::select(-c(pred_all_noC, sigma_Y_all_noC)) %>%   
  mutate(
    Point.Estimate_pred = case_when(
      is.na(Point.Estimate_Female) ~ pred_all, 
      TRUE ~ Point.Estimate_Female
    ),
    SE_pred = case_when(
      is.na(Point.Estimate_Female) ~ sqrt(SE_val_Female^2 + sigma_Y_all), 
      TRUE ~ SE_val_Female
    )
  ) 

plot(surv_wide_Female$Point.Estimate_Both, surv_wide_Female$Point.Estimate_pred)
lines(c(0,100),c(0,100))

if(FALSE){
  ggplot(data=surv_wide_Female, aes(y = SE_pred,
                                    x = SE_val_Female)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    facet_grid(~is.na(Point.Estimate_Female)) +
    geom_abline(intercept = 0, slope = 1)
  
  ggplot(data=surv_wide_Female, aes(y = Point.Estimate_pred,
                                    x = Point.Estimate_Both)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    facet_grid(~is.na(Point.Estimate_Female)) +
    geom_abline(intercept = 0, slope = 1)
  
  surv_wide_Female <- surv_wide_Female %>% 
    dplyr::select(
      c(ISO.code, country,  year, Sex,  UNICEFSurveyID,
        Point.Estimate_pred, SE_pred, ShortSource)
    )
}

test_merge <- surv_data %>% 
  left_join(surv_wide_Female, by = c(".imp","ISO.code", "year", "Sex", "country", "UNICEFSurveyID", "All_africa_HIW")) %>% 
  mutate(ShortSource = ShortSource.x) %>% 
  mutate(
    ShortSource = case_when(
      is.na(ShortSource) ~ ShortSource.y,
      TRUE ~ ShortSource
    )
  ) %>% 
  dplyr::select(-c(ShortSource.x,ShortSource.y))

















### Predicting for Males
# Fitting the model
fitted_model <- lme(formula,
                    random=~1|country,
                    weights = varIdent(form=~1|Source),
                    data = surv_wide, na.action = na.omit, 
                    control = cnt_list)

all_data <- surv_wide
data_frame <- surv_wide %>% 
  mutate(Y = sqrt(Point.Estimate_Male)) %>% 
  filter(!is.na(Y))

Xi.matrix <- cbind(model.matrix(formula_x ,data=data_frame))
SE.matrix <- cbind(model.matrix(~Source,data=data_frame)[,-1])
N <- length(unique(data_frame$country))
X_orig <- Xi.matrix
Z_orig <- make_Z_noPcov(data_frame,data_frame)

# Estimate EBLUPs, predicted means, and variances of the predicted means.
pred_mat <- pred_fun_noPcov(X=X_orig,Z=Z_orig,Y=data_frame$Y,fitted_model,N,SE_var=SE.matrix,full_b = 1)
S <- pred_mat$S
b <- pred_mat$b


# Extracting only those countries that were used in the original analysis
cov_data_used_only <- all_data[all_data$country %in% data_frame$country,]
# Ordering the data so the random effects align.
div=c(cov_data_used_only$country)
cov_data_used_only <- cov_data_used_only[order(div,cov_data_used_only$year),]

# Making the X and Z matrices for the new data.
X <- cbind(model.matrix(formula_x,data=cov_data_used_only))
SE<- cbind(model.matrix(~Source,data=cov_data_used_only)[,-1])
Z <- make_Z_noPcov(cov_data_used_only,data_frame)

# Predicting the mean and variance for each new observation.
pred_mat_new <- pred_fun_noPcov(X,Z,Y=0,fitted_model,N,SE_var=SE,full_b = b, S = S)
pred_all  <- pred_mat_new$mu_T^2
sigma_Y_all <- diag(pred_mat_new$sigma_Y)*4*pred_all

pred_new_data <- data.frame(cov_data_used_only,pred_all, sigma_Y_all)
all_data <- left_join(all_data,pred_new_data)
plot(all_data$Point.Estimate_Both, all_data$pred_all)
lines(c(0,100),c(0,100))


# Extracting countries that were not used in the original analysis
cov_data_not_used <- all_data[not_in(all_data$country,data_frame$country),]

pred_all_noC <- 0
sigma_Y_all_noC <- 0
if(nrow(cov_data_not_used)>0){
  # Predicting for them
  Xi.full <- cbind(model.matrix(formula_x,data=all_data))
  SE.full <- cbind(model.matrix(~Source,data=all_data)[,-1])
  X <- Xi.full[not_in(all_data$country,data_frame$country),]
  SE<- SE.full[not_in(all_data$country,data_frame$country),]
  
  N <- length(unique(cov_data_not_used$country))
  pred_mat_new_noC <- pred_fun_noPcov_noC(X,Y,fitted_model,N,SE_var=SE,full_b=b,S=S)
  pred_all_noC  <- pred_mat_new_noC$mu_T^2
  sigma_Y_all_noC <- diag(pred_mat_new_noC$sigma_Y)*4*pred_all_noC
  
  
  # Predicting the standard errors.
  pred_new_data_noC <- data.frame(cov_data_not_used,pred_all_noC, sigma_Y_all_noC)
  all_data <- merge(all_data,pred_new_data_noC,all.x=TRUE)
  plot(all_data$Point.Estimate_Both, all_data$pred_all_noC)
  lines(c(0,100),c(0,100))
}

surv_wide_Male <- surv_wide_B %>% 
  left_join( all_data) %>%   
  mutate(
    Sex = "Male",
    pred_all = case_when(
      is.na(pred_all) ~ pred_all_noC, 
      TRUE ~ pred_all
    ),
    sigma_Y_all = case_when(
      is.na(sigma_Y_all) ~ sigma_Y_all_noC, 
      TRUE ~ sigma_Y_all
    )
  ) %>% 
  dplyr::select(-c(pred_all_noC, sigma_Y_all_noC)) %>%   
  mutate(
    Point.Estimate_pred = case_when(
      is.na(Point.Estimate_Male) ~ pred_all, 
      TRUE ~ Point.Estimate_Male
    ),
    SE_pred = case_when(
      is.na(Point.Estimate_Male) ~ sqrt(SE_val_Female^2 + sigma_Y_all), 
      TRUE ~ SE_val_Female
    )
  ) 

plot(surv_wide_Male$Point.Estimate_Both, surv_wide_Male$Point.Estimate_pred)
lines(c(0,100),c(0,100))

if(FALSE){
  ggplot(data=surv_wide_Male, aes(y = SE_pred,
                                  x = SE_val_Male)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    facet_grid(~is.na(Point.Estimate_Male)) +
    geom_abline(intercept = 0, slope = 1)
  
  ggplot(data=surv_wide_Male, aes(y = Point.Estimate_pred,
                                  x = Point.Estimate_Both)) + 
    geom_point() +
    geom_smooth(method = "lm") +
    facet_grid(~is.na(Point.Estimate_Male)) +
    geom_abline(intercept = 0, slope = 1)
}
surv_wide_Male <- surv_wide_Male %>% 
  dplyr::select(
    c(.imp,ISO.code, country,  year, Sex, UNICEFSurveyID,  
      Point.Estimate_pred, SE_pred, ShortSource)
  )








### Final merge of all data
final_merge <- test_merge %>% 
  left_join(surv_wide_Male, by = c(".imp","ISO.code", "country",  "year", "Sex", "UNICEFSurveyID")) %>%   
  mutate(ShortSource = ShortSource.x) %>% 
  mutate(
    ShortSource = case_when(
      is.na(ShortSource) ~ ShortSource.y,
      TRUE ~ ShortSource
    )
  ) %>% 
  dplyr::select(-c(ShortSource.x,ShortSource.y)) %>% 
  mutate(
    Point.Estimate = case_when(
      !is.na(Point.Estimate) ~ Point.Estimate, 
      !is.na(Point.Estimate_pred.x) ~ Point.Estimate_pred.x, 
      !is.na(Point.Estimate_pred.y) ~ Point.Estimate_pred.y, 
      TRUE ~ NA_real_
    ),
    SE_pred = case_when(
      !is.na(SE_val) ~ SE_val, 
      !is.na(SE_pred.x) ~ SE_pred.x, 
      !is.na(SE_pred.y) ~ SE_pred.y, 
      TRUE ~ NA_real_
    )
  ) %>% 
  dplyr::select(
    -c(Point.Estimate_pred.x, SE_pred.x, Point.Estimate_pred.y, SE_pred.y, Region, "SEV.Z")
  ) %>% 
  rename(
    Point.Estimate.Orig = PointEstimate,
    Region = All_africa_HIW,
    Point.Estimate.Imp = Point.Estimate
  ) %>% 
  dplyr::select(
    c(.imp,ISO.code, country, Region, year, Sex, Point.Estimate.Orig, 
      Point.Estimate.Imp, Point.Estimate.NS, SE_pred, SE_val, ShortSource,
      SDI, MCI, MCI_5_yr, WB.Income.Group, SEV, UNICEFSurveyID,
      WHOSurveyID,	EstimateType,	LowerLimit,	UpperLimit,	weighted_N,	
      unweighted_N,	StandardFootnotes,	ReportAuthorProvider,	
      FullSourceTitle,	N,	Standard.Error
    )
  ) %>% 
  arrange(.imp,country,year,Sex)


## Checking for duplicates
final_merge_nodup <- final_merge %>% 
  distinct()
final_merge
final_merge_nodup

# Exporting
saveRDS(final_merge,file = paste("1- Sample Input Data/",marker,"_data_w_cov_",month,".rds",sep=""))





