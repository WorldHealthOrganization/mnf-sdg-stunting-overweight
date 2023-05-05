

library(tidyverse)
library(mvtnorm)
library(readxl)
library(mice)

IHME_dat_2022 <- read_csv("1- Sample Input Data/GBD 2022 MCI and SDI.csv") %>% 
  select(c(location_name, year, MCI, SDI, ISO.code))


## Get all the countries that there are data for
## From JME
surv_data <- read_excel("1- Sample Input Data/JME_Country_Level_Input_Stunting_Mar23.xlsx")
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
wpp_dat <- read_csv("1- Sample Input Data/WPP2022_Demographic_Indicators_Medium.csv", 
                    show_col_types = FALSE) %>% 
  filter(!is.na(ISO3_code)) %>% 
  rename(year = Time, ISO3Code = ISO3_code) %>% 
  filter(year > 1983 & year < 2023)
big_wpp_dat <- wpp_dat %>% 
  right_join(All_countries, by = c("ISO3Code")) %>% 
  left_join(IHME_dat_2022, by = c("ISO3Code" = "ISO.code", "year")) %>% 
  mutate(location_name = case_when(
    is.na(location_name) ~ Country,
    TRUE ~ location_name
  ))

GDP_wide <- read_csv("1- Sample Input Data/API_NY.GDP.MKTP.CD_DS2_en_csv_v2_4701247.csv", 
                    show_col_types = FALSE)
GDP_long <- GDP_wide %>% 
  pivot_longer(
    `1960`:`2021`, names_to = "year",
    names_transform = list(year = as.numeric),
    values_to = "GDP") %>% 
  filter(year > 1982) %>% 
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
  filter(!is.na(ISO3Code))



##### This first imputation will be used to generate the country level means by imputation.
## Number of imputations
B <- 20
ini <- mice(all_cov, maxit = 0)

pred <- ini$pred
meth <- ini$method
pred["MCI",] <- pred["SDI",] <- pred["lGDP",] <- c(0,0,0,rep(1,14))
diag(pred) <- 0

imp <- mice(all_cov, pred = pred, meth = meth, print = TRUE,
            m = B, maxit = 20, seed = 238746)


all_cov_comp <- mice::complete(imp, action = "long")

all_cov_comp <- all_cov_comp %>% 
  as_tibble() %>% 
  group_by(.imp,ISO3Code) %>% 
  mutate(mn_MCI = mean(MCI),
         mn_SDI = mean(SDI),
         mn_lGDP = mean(lGDP)) %>% 
  ungroup() %>% 
  select(-c("MCI","SDI","lGDP"))

##### This one will mean center the data by country and re-impute using wide format.

### First for MCI
all_cov_center_MCI_wide <- all_cov %>% 
  group_by(ISO3Code) %>% 
  mutate(Z_MCI = MCI - mean(MCI, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(c("ISO3Code","c_factor", "location_name", "year", "Z_MCI")) %>%
  pivot_wider(names_from = year, values_from = Z_MCI, 
              names_glue = "{.value}_{year}",
              values_fn = {mean})

ini <- mice(all_cov_center_MCI_wide, maxit = 0)

pred <- ini$pred
meth <- ini$method
imp_mci <- mice(all_cov_center_MCI_wide, 
                pred = pred, meth = meth, print = TRUE,
            m = 20, maxit = 5, seed = 238746)

all_cov_comp_wide <- mice::complete(imp_mci, action = "long")
all_cov_comp_wide_long <- all_cov_comp_wide %>% 
  pivot_longer(
    `Z_MCI_1984`:`Z_MCI_2022`, names_prefix = "Z_MCI_",
    names_to = "year",
    names_transform = list(year = as.numeric),
    values_to = "Z_MCI")

### Now for SDI
all_cov_center_SDI_wide <- all_cov %>% 
  group_by(ISO3Code) %>% 
  mutate(Z_SDI = SDI - mean(SDI, na.rm = TRUE)) %>% 
  ungroup() %>% 
  select(c("ISO3Code","c_factor", "location_name", "year", "Z_SDI")) %>%
  pivot_wider(names_from = year, values_from = Z_SDI, 
              names_glue = "{.value}_{year}",
              values_fn = {mean})

ini <- mice(all_cov_center_SDI_wide, maxit = 0)

pred <- ini$pred
pred["Z_SDI_2022",] <- c(0,0,0, rep(1,length(pred["Z_SDI_2022",]) - 4), 4)
meth <- ini$method
meth[length(meth)] = "pmm"
imp_sdi <- mice(all_cov_center_SDI_wide, pred = pred, meth = meth, 
                print = TRUE,
            m = 20, maxit = 5, seed = 238746)

all_cov_center_SDI_wide <- mice::complete(imp_sdi, action = "long")
all_cov_comp_wide_long <- all_cov_center_SDI_wide %>% 
  pivot_longer(
    `Z_SDI_1984`:`Z_SDI_2022`, names_prefix = "Z_SDI_",
    names_to = "year",
    names_transform = list(year = as.numeric),
    values_to = "Z_SDI") %>% 
  full_join(all_cov_comp_wide_long)



### Final data 
fin_all_cov <- all_cov_comp %>% 
  left_join(
    all_cov_comp_wide_long %>% 
      select(".imp","ISO3Code","year", "Z_MCI","Z_SDI")
  ) %>% 
  mutate(
    imp_MCI = mn_MCI + Z_MCI,
    imp_SDI = mn_SDI + Z_SDI,
  ) %>% 
  select(c(".imp","ISO3Code","year","imp_MCI","imp_SDI", "location_name")) %>% 
  left_join(all_cov) %>% 
  select(c(".imp","ISO3Code","year","MCI","SDI","imp_MCI","imp_SDI","lGDP", "location_name")) 

miss_coun <- fin_all_cov %>% 
  mutate(
    mci_diff = MCI - imp_MCI, 
    sdi_diff = SDI - imp_SDI
    ) %>% 
  filter(is.na(SDI))
#ggplot(miss_coun, aes(x = year, y = imp_SDI, group = ISO3Code)) + geom_line() + facet_wrap(~.imp)
#ggplot(miss_coun, aes(x = year, y = imp_MCI, group = ISO3Code)) + geom_line() + facet_wrap(~.imp)


all_cov <- fin_all_cov %>% 
  mutate(
    MCI = imp_MCI,
    SDI = imp_SDI) %>% 
  select(
    c(".imp","ISO3Code","year","MCI","SDI","lGDP", "location_name")
    ) 



cov_data <- readRDS("IHME_covs/Mean_Imputed_Tidy_World_Devel_Indic_Nov2022.rds")
cov_data <- cov_data %>% 
  distinct(
    ISO.code, year, location_name, sex_name,
    .keep_all = TRUE) %>% 
  rename("Sex" = sex_name) %>% 
  mutate(
    SEV = SEV_MI
  ) %>% 
  select(-c("SEV_MI","MCI","SDI"))

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
  group_by(.imp,ISO.code, Sex) %>% 
  mutate(lag1 = lag(MCI), 
         lag2 = lag(MCI, 2),
         lag3 = lag(MCI, 3),
         lag4 = lag(MCI, 4),
         lag5 = lag(MCI, 5)
  )  %>% 
  mutate(MCI_5_yr = (lag1+lag2+lag3+lag4+lag5)/5) %>% 
  ungroup() %>% 
  select(-c("lag1","lag2","lag3","lag4","lag5")) %>% 
  filter(year > 1989)

here_coun <- cov_data %>% 
  filter(ISO.code == "XKX" | ISO.code == "TCA" )

ggplot(here_coun, 
       aes(x = year, y = SDI, group = .imp)) +
  geom_line() +
  facet_grid(~ISO.code)
ggplot(here_coun, 
       aes(x = year, y = MCI_5_yr, group = .imp)) +
  geom_line() +
  facet_grid(~ISO.code)

saveRDS(cov_data,"IHME_covs/Imputed_IHME_Mar2023.rds")



















