
wd <- "~path to root directory"
setwd(wd)


library(nlme)
library(tidyverse)
library(readxl)

###################################################################
##### First choose the file to run, it works the same for all #####
###################################################################

marker <- "Stunting"


if(marker=="Stunting"){surv_data <- read_excel("1- Sample Input Data/JME_Country_Level_Input_Stunting_Mar23.xlsx")}
if(marker=="Overweight"){surv_data <- read_excel("1- Sample Input Data/JME_Country_Level_Input_Overweight_Mar23.xlsx")}
if(marker=="Underweight"){surv_data <- read_excel("1- Sample Input Data/JME_Country_Level_Input_Underweight_Mar23.xlsx")}
if(marker=="Severe_Wasting"){surv_data <- read_excel("1- Sample Input Data/JME_Country_Level_Input_Severe_Wasting_Mar23.xlsx")}
if(marker=="Wasting"){surv_data <- read_excel("1- Sample Input Data/JME_Country_Level_Input_Wasting_Mar23.xlsx")}

## Check to see all variables have the correct format and the number of estimates
str(surv_data)
dim(surv_data)
length(surv_data$PointEstimate[!is.na(surv_data$PointEstimate)])

## Export excluded estimates with really small N
surv_data_sm = surv_data %>% 
  filter(!is.na(unweighted_N) & unweighted_N < 6)
write_csv(surv_data_sm,paste0("1- Sample Input Data/Country_Level_Input_",marker,"_smallN_Mar23.csv"))
surv_data = surv_data %>% 
  filter(is.na(unweighted_N) | unweighted_N > 5)

## Number of estimates after exclusions
dim(surv_data)

# First we'll limit the number of sources by lumping many sources
# with not many survey's into an "other" category.  The sources
# lumped into the "other" category are CFSVA, NHANES, SDHS, IDHS, RHS
# VAD, Surveillance, CDHS, LSMS, PAPFAM, SMART, VAD Survey and EHIS.

table(surv_data$ShortSource)
surv_data$Source = surv_data$ShortSource
surv_data$Source[surv_data$Source %in% c("CDHS","CFSVA","CWIQ", "EHIS","GFHS","IDHS", "LSMS","NHANES", "PAPFAM","Priority Survey","RHS", "SDHS","Surveillance", "VAD Survey", "Priority Survey")] = "Other"  ### Note: no "Surveillance" surveys have an observed SSE

#Now merge like sources
surv_data$Source[surv_data$Source=="DHS-Style"] = "DHS"
surv_data$Source[surv_data$Source=="MICS-DHS"] = "MICS"
surv_data$Source[surv_data$Source=="MICS-Style"] = "MICS"
surv_data$Source <- factor(surv_data$Source,levels=unique(surv_data$Source),labels=unique(surv_data$Source))
table(surv_data$Source)

## Ordering data
surv_data <- surv_data[order(surv_data$Country,surv_data$Year),]

# Get the sample size
surv_data$N <- surv_data$unweighted_N
surv_data$Standard.Error <- surv_data$StandardError

#### Getting SE's from confidence intervals ####
## The following code will attempt to get SE's from confidence intervals.
## This is only used in situations where confidence intervals are available
## and the SE is not. Confidence intervals will usually be of the form:
## estimate +/- 1.96*SE. However, they are commonly back-transformed to be 
## on the (0,1) scale. Here, we test different transformations, find the 
## one that is the most symmetric, and use it to get an estimate of the SE.

# Trying 3 different transformation for getting symmetric confidence intervals.
lower_inter <- surv_data$PointEstimate/100 - surv_data$LowerLimit/100
upper_inter <- surv_data$UpperLimit/100 - surv_data$PointEstimate/100
lower_inter2 <- log(surv_data$PointEstimate/100) - log(surv_data$LowerLimit/100)
upper_inter2 <- log(surv_data$UpperLimit/100) - log(surv_data$PointEstimate/100)
logit <- function(p){log(p/(1-p))}
lower_inter3 <- logit(surv_data$PointEstimate/100) - logit(surv_data$LowerLimit/100)
upper_inter3 <- logit(surv_data$UpperLimit/100) - logit(surv_data$PointEstimate/100)

# Finding which transformation results in the most symmetric confidence interval.
diff_matrix <- abs(cbind(lower_inter-upper_inter,lower_inter2-upper_inter2,lower_inter3-upper_inter3))
diff_matrix[is.na(surv_data$LowerLimit) | is.na(lower_inter2)] <- 0
best_trans <- c(apply(t(diff_matrix),2,which.min))
best_min <- c(apply(t(diff_matrix),2,min))

# Distribution of which transformation is the best
table(best_trans[!is.na(surv_data$LowerLimit)])


# looking at how symmetric the best transformation is
summary(best_min[!is.na(surv_data$LowerLimit)])

# Seeing if any PointEstimate values are outside their confidence interval
surv_data %>% 
  filter((PointEstimate  < LowerLimit ) | (PointEstimate  > UpperLimit ) ) %>% 
  summarise(n())

# Estimating the SE's on the transformed scale
SE_matrix_trans <- cbind(lower_inter+upper_inter,lower_inter2+upper_inter2,lower_inter3+upper_inter3)/(2*1.96)

# Estimating the SE's on the non-transformed scale using the delta method
SE_matrix <- SE_matrix_trans
SE_matrix[,1] <- SE_matrix_trans[,1]*100
SE_matrix[,2] <- SE_matrix_trans[,2]*surv_data$PointEstimate
SE_matrix[,3] <- SE_matrix_trans[,3]*surv_data$PointEstimate*(1-surv_data$PointEstimate/100)

# Seeing if the values are relatively similar to the observed SE's
summary(cbind(SE_matrix,surv_data$Standard.Error))

# Getting the best non-transformed SE and comparing those to the observed SE's
SE_vec <- SE_matrix[cbind(1:(dim(SE_matrix)[1]),best_trans)]
summary(cbind(SE_vec,surv_data$Standard.Error,SE_vec/surv_data$Standard.Error))

# Changing surveys where confidence intervals are available and the SE is not
# to the estimated SE.
surv_data$Standard.Error[is.na(surv_data$Standard.Error) & !is.na(surv_data$LowerLimit)] <- SE_vec[is.na(surv_data$Standard.Error) & !is.na(surv_data$LowerLimit)] 

# Getting the number of observed SE's and total SE's. 
length(surv_data$PointEstimate[!is.na(surv_data$PointEstimate)])
dim(surv_data)

###############################################################
##### In this first part "outlying SE's are identified       ##
##### as SE with standard residuals that are less than       ##
##### -4. These are repeatedly removed until no standardized ##
##### residuals are less than -4 (only remove SE's if they   ##
##### are too low).                                          ##
###############################################################

# create the p(1-p) values, changing SE = 0 to NA. 
surv_data <- surv_data %>% 
  mutate(
    Point.Estimate = PointEstimate,
    P1mP  = PointEstimate/100*(1-PointEstimate/100),
         Standard.Error = case_when(
           Standard.Error == 0 ~ NA_real_,
           TRUE ~ Standard.Error
         ))
summary(surv_data$Point.Estimate)

# create a dataset with only complete values.
cc_surv_data <- surv_data %>% 
  filter(!is.na(Standard.Error) & !is.na(N))
nrow(cc_surv_data)

# what percentage of surveys have complete data.
nrow(cc_surv_data)/nrow(surv_data)

## Look at some plots of the data
par(mar = c(4,5,1,1))
plot(
  log(cc_surv_data$N),
  log(cc_surv_data$Standard.Error),
  pch=19,
  xlab = expression(log(n[ij])),
  ylab = expression(log(S[ij])),
  cex.lab=1.4,col=as.integer(cc_surv_data$Source)
  )
legend("topright",
       legend = unique(cc_surv_data$Source),
       col = unique(as.integer(cc_surv_data$Source)),
       pch=19)


removed <- cc_surv_data %>% filter(Point.Estimate == 0)
pred_rm <- NULL
resid_rm <- NULL
cc_surv_data_rm <- cc_surv_data %>% filter(Point.Estimate>0)
C <- -4
L <- 1
while (L>0){
  fit <- lme(
    log(Standard.Error)~log(P1mP) + log(N),
    data = cc_surv_data_rm,
    random=~1|Country,
    weights = varIdent(~as.factor(Source))
    )
  resid <- residuals(fit,type = "pearson")
  pred <- predict(fit, cc_surv_data_rm)
  pred_rm <- c(pred_rm, as.vector(pred[resid < C]))
  resid_rm <- c(resid_rm, resid[resid < C])
  
  ####Removing outliers####
  L <- dim(cc_surv_data_rm[resid < C,])[1]
  removed <- rbind(removed,cc_surv_data_rm[resid < C,])
  cc_surv_data_rm <- cc_surv_data_rm[resid > C,]
}

## SE's removed from outlier procedure.
View(removed)

# Now we create the SE variable to use for modeling and set "removed" SE's to NA
surv_data$SE <- surv_data$Standard.Error
surv_data$SE[surv_data$N %in% removed$N & surv_data$Point.Estimate %in% removed$Point.Estimate] <- NA

length(surv_data$PointEstimate[!is.na(surv_data$PointEstimate)])




#####Predicting with standard errors#######
##### Recode the reference level and rename some variables.
surv_data <-  within(surv_data, Source <- relevel(Source, ref = "SMART"))
surv_data <- surv_data %>% 
  rename(country = Country, year = Year)

##### Predicting for those with P1mP and N.
fitted_model <- lme(log(SE)~Source + log(P1mP) + log(N),random=~1|country,weights = varIdent(form=~1|Source),data = surv_data,na.action = na.omit)

miss_ind <- complete.cases(surv_data %>% select(c(Source, P1mP, N, country, Source)))
pred <- pred2 <- rep(NA, length(miss_ind))
pred[miss_ind] <- predict(fitted_model, newdata = surv_data, na.action = na.omit)
pred2 <- predict(fitted_model, newdata = surv_data, na.action = na.omit, level = 0)
pred[is.na(pred)] = pred2[is.na(pred)]
pred[surv_data$P1mP==0] <- NA
surv_data = surv_data %>% bind_cols(Pred_SE = exp(pred))


##### Predicting for those with P1mP, and no N but a weighted N #########

fitted_model <- lme(log(SE)~Source + log(P1mP) + log(weighted_N) + Source:log(P1mP),
                    random=~1|country,
                    weights = varIdent(form=~1|Source),
                    data = surv_data,
                    na.action = na.omit)

## predict for countries with other surveys in the data, then for those without.
miss_ind <- complete.cases(surv_data %>% select(c(Source, P1mP,weighted_N, country, Source)))
pred <- rep(NA, length(miss_ind))
pred[miss_ind] <- predict(fitted_model, newdata = surv_data, na.action = na.omit)
pred2 <- predict(fitted_model, newdata = surv_data, na.action = na.omit, level = 0)
pred[is.na(pred)] = pred2[is.na(pred)]
pred[surv_data$P1mP==0] <- NA
surv_data = surv_data %>% bind_cols(Pred_SE_noN = exp(pred))


##### Predicting for those with no N #########

fitted_model <- lme(log(SE)~Source + log(P1mP) + Source:log(P1mP),
                    random=~1|country,
                    weights = varIdent(form=~1|Source),
                    data = surv_data,
                    na.action = na.omit)

## predict for countries with other surveys in the data, then for those without.
miss_ind <- complete.cases(surv_data %>% select(c(Source, P1mP, country, Source)))
pred <- rep(NA, length(miss_ind))
pred[miss_ind] <- predict(fitted_model, newdata = surv_data, na.action = na.omit)
pred2 <- predict(fitted_model, newdata = surv_data, na.action = na.omit, level = 0)
pred[is.na(pred)] = pred2[is.na(pred)]
pred[surv_data$P1mP==0] <- NA
surv_data = surv_data %>% bind_cols(Pred_SE_zeroN = exp(pred))


##### Predicting for those with no log(P1mP) but with an N #########

fitted_model <- lme(log(SE)~Source+ log(N) + Source:log(N),random=~1|country,weights = varIdent(form=~1|Source),data = surv_data,na.action = na.omit)

## predict for countries with other surveys in the data, then for those without.
miss_ind <- complete.cases(surv_data %>% select(c(Source, N, country, Source)))
pred <- rep(NA, length(miss_ind))
pred[miss_ind] <- predict(fitted_model, newdata = surv_data, na.action = na.omit)
pred2 <- predict(fitted_model, newdata = surv_data, na.action = na.omit, level = 0)
pred[is.na(pred)] = pred2[is.na(pred)]
surv_data = surv_data %>% bind_cols(Pred_SE_noP1mP = exp(pred))

#### Setting the SE_val to the appropriate predicted value.
surv_data$SE_val <- surv_data$Standard.Error
surv_data$SE_val[is.na(surv_data$SE_val)] <- surv_data$Pred_SE[is.na(surv_data$SE_val)]
surv_data$SE_val[is.na(surv_data$SE_val)] <- surv_data$Pred_SE_noN[is.na(surv_data$SE_val)]
surv_data$SE_val[is.na(surv_data$SE_val)] <- surv_data$Pred_SE_zeroN[is.na(surv_data$SE_val)]
surv_data$SE_val[is.na(surv_data$SE_val)] <- surv_data$Pred_SE_noP1mP[is.na(surv_data$SE_val)]

#### Final look at the SE's
par(mar = c(4,5,1,1))
plot(
  log(surv_data$N),
  log(surv_data$SE_val),
  pch=19,
  xlab = expression(log(n[ij])),
  ylab = expression(log(S[ij])),
  cex.lab=1.4,col=as.integer(surv_data$Source)
)
legend("topright",
       legend = unique(surv_data$Source),
       col = unique(as.integer(surv_data$Source)),
       pch=19)

#### Final data
final_data <- surv_data %>% 
  select(
    -c("PointEstimate","StandardError","Source","SE",
       "Pred_SE","Pred_SE_noN","Pred_SE_noP1mP","Pred_SE_zeroN","P1mP")
    ) 


#View(final_data)

### Export data ###
if(marker=="Stunting"){saveRDS(final_data,"1- Sample Input Data/Global_stunt_SE_clean.rds")}
if(marker=="Overweight"){saveRDS(final_data,"1- Sample Input Data/Global_over_SE_clean.rds")}
if(marker=="Underweight"){saveRDS(final_data,"1- Sample Input Data/Global_under_SE_clean.rds")}
if(marker=="Severe_Wasting"){saveRDS(final_data,"1- Sample Input Data/Global_sv_waste_SE_clean.rds")}
if(marker=="Wasting"){saveRDS(final_data,"1- Sample Input Data/Global_waste_SE_clean.rds")}


