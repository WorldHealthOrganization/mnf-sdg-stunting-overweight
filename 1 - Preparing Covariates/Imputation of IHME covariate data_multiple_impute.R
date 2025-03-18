# Clear environment
rm(list = ls())

# Set-up environment
library(this.path)
wd <- dirname(this.path::here())
print(wd)
setwd(wd)

source("Utils/Programs_Feb_2020.R")
library(tidyverse)
library(mvtnorm)
library(readxl)
library(mice)
library(doParallel)  

###### 1. Read in and Merge Data   ############### 

## Read in IHME covariate data
last_year <- 2025
IHME_dat_2022 <- read_csv("Data/IHME_covs/GBD 2022 MCI and SDI.csv") %>% 
  select(c(location_name, year, MCI, SDI, ISO.code))


## Get all the countries that there are data for
## From JME
surv_data <- read_excel("Data/JME/2024/Raw/JME_Country_Level_Input_Stunting.xlsx")
surv_data <- surv_data %>% 
  group_by(ISO3Code) %>% 
  filter(row_number() == 1) %>% 
  select(c(ISO3Code, Country))
# From IHME data
IHME_countries <- IHME_dat_2022 %>% 
  rename(ISO3Code = ISO.code, Country = location_name) %>% 
  group_by(ISO3Code) %>% 
  filter(row_number() == 1) %>% 
  select(c(ISO3Code))

# Full merge
All_countries <- surv_data %>% 
  full_join(IHME_countries, by = "ISO3Code") %>% 
  bind_rows(
    tibble(ISO3Code = "XKX", Country = "Kosovo")
  )

## Merge with Fertility and GDP data
wpp_dat <- read_csv("Data/WPP/WPP2022_Demographic_Indicators_Medium.csv", 
                    show_col_types = FALSE) %>% 
  filter(!is.na(ISO3_code)) %>% 
  rename(year = Time, ISO3Code = ISO3_code) %>% 
  filter(year > 1983 & year <= last_year)
big_wpp_dat <- wpp_dat %>% 
  right_join(All_countries, by = c("ISO3Code")) %>% 
  left_join(IHME_dat_2022, by = c("ISO3Code" = "ISO.code", "year")) %>% 
  mutate(location_name = case_when(
    is.na(location_name) ~ Country,
    TRUE ~ location_name
  )
  )


GDP_wide <- read_csv("Data/WB/API_NY.GDP.MKTP.CD_DS2_en_csv_v2_4701247.csv", 
                     show_col_types = FALSE)
GDP_long <- GDP_wide %>% 
  pivot_longer(
    `1960`:`2021`, names_to = "year",
    names_transform = list(year = as.numeric),
    values_to = "GDP") %>% 
  filter(year > 1982 & year <= last_year) %>% 
  rename(ISO3Code = `Country Code`) %>% 
  select(c(ISO3Code, year, GDP))

all_cov <- big_wpp_dat %>% 
  left_join(GDP_long, by=c("ISO3Code", "year")) %>% 
  mutate(lGDP = log(GDP), 
         lpop = log(TPopulation1Jan),
         c_factor = factor(LocID)) %>% 
  select(c("ISO3Code","c_factor", "location_name", "year", "MCI", "SDI", "lGDP",
           "lpop","PopDensity", "MedianAgePop", "CBR", "TFR", 
           "Births1519", "CDR", "LEx", "IMR", "Q5")) %>% 
  mutate(location_name = case_when(
    year == 2023 ~ lag(location_name, n = 1),
    year == 2024 ~ lag(location_name, n = 2),
    year == 2025 ~ lag(location_name, n = 3),
    TRUE ~ location_name
  )
  ) %>%
  filter(!is.na(ISO3Code)) %>% 
  group_by(ISO3Code, year) %>% 
  filter(row_number() == 1) %>% 
  ungroup()

##### 2. Imputation to Create Country-Level Means by Imputation   ############### 
##### This first imputation will be used to generate the country level means by imputation.
B <- 20 # Number of imputations
ini <- mice(all_cov, maxit = 0)

pred <- ini$pred
meth <- ini$method
pred["MCI",] <- pred["SDI",] <- pred["lGDP",] <- c(0,0,0,
                                                   rep(1,length(pred["MCI",]) - 3)
                                                   )
diag(pred) <- 0

imp <- parlmice(data = all_cov, n.core = 5, n.imp.core = 4,
                pred = pred, meth = meth, print = TRUE,
                m = B, maxit = 20, cluster.seed = 238746)

all_cov_comp <- mice::complete(imp, action = "long")

all_cov_comp <- all_cov_comp %>% 
  as_tibble() %>% 
  group_by(.imp,ISO3Code) %>% 
  filter(year <= last_year) %>% 
  mutate(mn_MCI = mean(MCI),
         mn_SDI = mean(SDI),
         mn_lGDP = mean(lGDP)) %>% 
  ungroup() %>% 
  select(-c("MCI","SDI","lGDP"))



##### 3. Perform Imputation on Mean Centered Data  ############### 
##### Mean center the data by country and re-impute using wide format.
### First for MCI
all_cov_center_MCI_wide <- all_cov %>% 
  group_by(ISO3Code) %>% 
  mutate(Z_MCI = MCI - mean(MCI, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(c("ISO3Code","c_factor", "year", "Z_MCI")) %>%
  pivot_wider(names_from = year, values_from = Z_MCI, 
              names_glue = "{.value}_{year}",
              values_fn = {mean})

ini <- mice(all_cov_center_MCI_wide, maxit = 0)

pred <- ini$pred
meth <- ini$method

pred <- ini$pred
pred[1:3,] <- 0
pred[,1:3] <- 0

imp_mci <- parlmice(data = all_cov_center_MCI_wide, n.core = 7, n.imp.core = 3,
                    pred = pred, meth = meth, print = TRUE,
                    m = B, maxit = 20, cluster.seed = 238746)


all_cov_comp_wide <- mice::complete(imp_mci, action = "long")
all_cov_comp_wide_long <- all_cov_comp_wide %>% 
  pivot_longer(
    "Z_MCI_1984":paste0("Z_MCI_",last_year), names_prefix = "Z_MCI_",
    names_to = "year",
    names_transform = list(year = as.numeric),
    values_to = "Z_MCI")


### Now for SDI
all_cov_center_SDI_wide <- all_cov %>% 
  group_by(ISO3Code) %>% 
  mutate(Z_SDI = SDI - mean(SDI, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(c("ISO3Code","c_factor", "year", "Z_SDI")) %>%
  pivot_wider(names_from = year, values_from = Z_SDI, 
              names_glue = "{.value}_{year}",
              values_fn = {mean})

ini <- mice(all_cov_center_SDI_wide, maxit = 0)

pred <- ini$pred
pred[1:3,] <- 0
pred[,1:3] <- 0
meth <- ini$method
meth[names(meth)=="Z_SDI_2022"] = "pmm"

imp_sdi <- parlmice(data = all_cov_center_SDI_wide, n.core = 7, n.imp.core = 3,
                    pred = pred, meth = meth, print = TRUE,
                    m = B, maxit = 20, cluster.seed = 238746)

all_cov_center_SDI_wide <- mice::complete(imp_sdi, action = "long")
all_cov_comp_wide_long <- all_cov_center_SDI_wide %>% 
  pivot_longer(
    "Z_SDI_1984":paste0("Z_SDI_",last_year), names_prefix = "Z_SDI_",
    names_to = "year",
    names_transform = list(year = as.numeric),
    values_to = "Z_SDI") %>% 
  full_join(
    all_cov_comp_wide_long
    )

save.image("Data/IHME_covs/Imputation_1Aug24.RData")



############### 4. Perform Imputation for after 2022 for All Countries ############### 
### Run to impute for last year 
### First for MCI
P_all_data <- all_cov_comp %>% 
  left_join(
    all_cov_comp_wide_long %>% 
      select(".imp","ISO3Code","year", "Z_MCI","Z_SDI")
  ) %>% 
  mutate(
    imp_MCI = mn_MCI + Z_MCI,
    imp_SDI = mn_SDI + Z_SDI,
  ) %>% 
  mutate(
    country = ISO3Code, 
    Y= imp_MCI,
    SE_var = 1)




no_cores <- 4 # detectCores()
registerDoParallel(cores = no_cores)  
cl <- makeCluster(no_cores) 

settings_to_try <- 1:B
foreach(set_i = settings_to_try)  %dopar% {
  
  j <- set_i
  all_data <- P_all_data %>% 
    filter(.imp == j & year > 2013)
  data_w_out <- all_data %>% 
    select(c("country", "year", "Y", "SE_var")) 
  data_w_out$SE_pred <- data_w_out$SE_var
  
  DF_R <- quantile(data_w_out$year,probs = c(0.5))
  B.knots <- range(data_w_out$year[!is.na(data_w_out$Y)])
  B.knots[1] <- B.knots[1] - 1
  B.knots[2] <- B.knots[2] + 10
  DF_P <- seq(min(data_w_out$year),2022,2)
  cov_data <- as.matrix(data.frame(Sex = all_data[,c("lpop")], lpop2 = all_data[,c("lpop")]^2))
  colnames(cov_data)[1] <- "Sex"  
  zero_covs <- NULL
  cov_mat <- "VC"
  q.order <- 2
  
  ##### Covariate analysis with multiple penalized functions
  Estimation <- cmnpe(data_w_out, DF_P, DF_R, B.knots, q.order, 
                      cov_data=cov_data, Pcov_data = NULL, 
                      cov_mat = cov_mat, plots = FALSE, TRANS=FALSE,
                      zero_covs = zero_covs, slope = TRUE)
  
  t_plot_data <- Estimation$plot_data
  saveRDS(t_plot_data, file = paste0("Data/IHME_covs/Plot data for MCI imputation ",j,".rds"))
  
}

save.image("Data/IHME_covs/Imputation_1Aug24.RData")

plot_dat <- tibble()
for(j in 1:B){
  
  t_plot_data <- readRDS(file = paste0("Data/IHME_covs/Plot data for MCI imputation ",j,".rds"))
  t_plot_data <- t_plot_data %>% 
    mutate(.imp = j) %>% 
    filter(year >= 2023) %>% 
    mutate(ISO3Code = country) %>% 
    group_by(ISO3Code,year) %>% 
    filter(row_number()==1) %>% 
    ungroup() %>% 
    select(.imp,ISO3Code, year,pred,sigma_Y_est)
  
  t_plot_data$MCI_imp <- t_plot_data$pred + 
    rnorm(length(t_plot_data$sigma_Y_est),0,t_plot_data$sigma_Y_est)
  
  t_plot_data <- t_plot_data %>% 
    select(.imp,ISO3Code, year,MCI_imp)
  plot_dat <- plot_dat %>% 
    bind_rows(t_plot_data)
}


P_all_data <- P_all_data %>% 
  left_join(plot_dat) %>% 
  mutate(
    imp_MCI = case_when(
      year >= 2023 ~ MCI_imp,
      TRUE ~ imp_MCI
    ))



### Now for SDI
P_all_data <- P_all_data %>% 
  mutate(
    country = ISO3Code, 
    Y= log(imp_SDI/1.15/(1-imp_SDI/1.15)),
    SE_var = 1) %>% 
  group_by(.imp,ISO3Code,year) %>% 
  filter(row_number()==1) %>% 
  ungroup()

settings_to_try <- 1:B
foreach(set_i = settings_to_try)  %dopar% {
  
  j <- set_i
  all_data <- P_all_data %>% 
    filter(.imp == j & year > 2013)
  data_w_out <- all_data %>% 
    select(c("country", "year", "Y", "SE_var")) 
  data_w_out$SE_pred <- data_w_out$SE_var
  
  DF_R <- quantile(data_w_out$year,probs = c(0.5))
  B.knots <- range(data_w_out$year[!is.na(data_w_out$Y)])
  B.knots[1] <- B.knots[1] - 1
  B.knots[2] <- B.knots[2] + 10
  DF_P <- seq(min(data_w_out$year),2022,2)
  cov_data <- as.matrix(data.frame(Sex = all_data[,c("lpop")], lpop2 = all_data[,c("lpop")]^2))
  colnames(cov_data)[1] <- "Sex"  
  zero_covs <- NULL
  cov_mat <- "VC"
  q.order <- 2
  
  ######### Covariate analysis with multiple penalized functions
  Estimation <- cmnpe(data_w_out, DF_P, DF_R, B.knots, q.order, 
                      cov_data=cov_data, Pcov_data = NULL, 
                      cov_mat = cov_mat, plots=FALSE, TRANS=FALSE,
                      zero_covs = zero_covs, slope = TRUE)
  cat(j,2*Estimation$df - 2*c(summary(Estimation$result$model)$logLik) + summary(Estimation$result$model)$AIC+2*c(summary(Estimation$result$model)$logLik),"\n")
  
  t_plot_data <- Estimation$plot_data
  saveRDS(t_plot_data, file = paste0("Data/IHME_covs/Plot data for SDI imputation ",j,".rds"))
  
}

save.image("Data/IHME_covs/Imputation_1Aug24.RData")

plot_dat <- tibble()
for(j in 1:B){
  
  t_plot_data <- readRDS(file = paste0("Data/IHME_covs/Plot data for SDI imputation ",j,".rds"))
  t_plot_data <- t_plot_data %>% 
    mutate(.imp = j) %>% 
    filter(year >= 2023) %>% 
    mutate(ISO3Code = country) %>% 
    group_by(ISO3Code,year) %>% 
    filter(row_number()==1) %>% 
    ungroup() %>% 
    select(.imp,ISO3Code, year,pred,sigma_Y_est)
  
  t_plot_data$SDI_imp <- exp(t_plot_data$pred + 
                               rnorm(length(t_plot_data$sigma_Y_est),0,t_plot_data$sigma_Y_est))
  t_plot_data$SDI_imp <- 1.15*t_plot_data$SDI_imp/(1 + t_plot_data$SDI_imp)
  
  t_plot_data <- t_plot_data %>% 
    select(.imp,ISO3Code, year,SDI_imp)
  plot_dat <- plot_dat %>% 
    bind_rows(t_plot_data)
}


P_all_data <- P_all_data %>% 
  left_join(plot_dat) %>% 
  mutate(
    imp_SDI = case_when(
      year >= 2023 ~ SDI_imp,
      TRUE ~ imp_SDI
    ))





##### 5. Merge Data and Create Final Variables  ############### 
### Final data 
fin_all_cov <- P_all_data %>% 
  select(c(".imp","ISO3Code","year","imp_MCI","imp_SDI", "location_name")) %>% 
  left_join(
    all_cov
  ) %>% 
  select(
    c(".imp","ISO3Code","year","MCI","SDI","imp_MCI","imp_SDI","lGDP", "location_name")
  ) 

### Checking model fit
samp_coun <- fin_all_cov %>% filter(ISO3Code %in% sample(unique(fin_all_cov$ISO3Code),20))

ggplot(samp_coun, aes(x = year, y = imp_SDI, group = .imp)) + 
  geom_line(show.legend = FALSE) + 
  facet_wrap(~location_name)

ggplot(samp_coun, aes(x = year, y = imp_MCI, group = .imp)) + 
  geom_line(show.legend = FALSE) + 
  facet_wrap(~location_name)


here_coun <- fin_all_cov %>% 
  filter(ISO3Code == "XKX" | ISO3Code == "TCA" )

ggplot(here_coun, 
       aes(x = year, y = imp_SDI, group = .imp, color = .imp)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~ISO3Code)
ggplot(here_coun, 
       aes(x = year, y = imp_MCI, group = .imp, color = .imp)) +
  geom_line(show.legend = FALSE) +
  facet_grid(~ISO3Code)

all_cov <- fin_all_cov %>% 
  mutate(
    MCI = imp_MCI,
    SDI = imp_SDI) %>% 
  select(
    c(".imp","ISO3Code","year","MCI","SDI","lGDP", "location_name")
  ) %>% 
  group_by(.imp,ISO3Code,year) %>% 
  filter(row_number()==1) %>% 
  ungroup()



cov_data <- readRDS("Data/IHME_covs/Tidy_World_Devel_Indic.rds")

all_cov_f <- all_cov %>% 
  mutate(Sex = "Female")
all_cov_m <- all_cov %>% 
  mutate(Sex = "Male")
all_cov <- all_cov %>% 
  mutate(Sex = "Both") %>% 
  bind_rows(all_cov_m) %>% 
  bind_rows(all_cov_f) %>% 
  arrange(ISO3Code, year, Sex)

P_cov_data <- cov_data %>% 
  full_join(all_cov, by = c("ISO.code"="ISO3Code", "year", "Sex")) %>% 
  mutate(location_name = case_when(
    !is.na(location_name.x) ~ location_name.x,
    !is.na(location_name.y) ~ location_name.y,
    TRUE ~ as.character(Country)
  )) %>% 
  select(-c("location_name.x", "location_name.y"))

# Create the 5 year average for each country
cov_data <- P_cov_data %>% 
  arrange(.imp,ISO.code, Sex, year) %>% 
  mutate(lag1 = lag(MCI), 
         lag2 = lag(MCI, 2),
         lag3 = lag(MCI, 3),
         lag4 = lag(MCI, 4),
         lag5 = lag(MCI, 5)
  )  %>% 
  mutate(MCI_5_yr = (lag1+lag2+lag3+lag4+lag5)/5) %>% 
  select(-c("lag1","lag2","lag3","lag4","lag5")) %>% 
  filter(year > 1989) %>% 
  group_by(.imp,ISO.code, Sex) %>% 
  fill(Region:Country) %>% 
  ungroup()

here_coun <- cov_data %>% 
  filter(ISO.code == "XKX" | ISO.code == "TCA" )

ggplot(here_coun, 
       aes(x = year, y = MCI_5_yr, group = .imp)) +
  geom_line() +
  facet_grid(~ISO.code)

samp_coun <- cov_data %>% filter(ISO.code %in% sample(unique(cov_data$ISO.code),36))

ggplot(samp_coun, aes(x = year, y = MCI_5_yr, group = .imp)) + 
  geom_line(show.legend = FALSE) + 
  facet_wrap(~location_name)
ggplot(samp_coun, aes(x = year, y = MCI, group = .imp)) + 
  geom_line(show.legend = FALSE) + 
  facet_wrap(~location_name)

saveRDS(cov_data,"Data/IHME_covs/Multiple_Imputed_Mar2023.rds")

cov_data_mean <- cov_data %>% 
  group_by(ISO.code, Sex,year, Region, Country, location_name) %>% 
  summarise(
    MCI_5_yr = mean(MCI_5_yr),
    SDI = mean(SDI),
    MCI = mean(MCI),
    lGDP = mean(lGDP, na.rm=TRUE),
    year = mean(year),
    SEV = mean(SEV)
  ) %>% 
  ungroup()



saveRDS(cov_data,"Data/IHME_covs/Single_Impute_Mar2023.rds")


# Clear environment
rm(list = ls())

















