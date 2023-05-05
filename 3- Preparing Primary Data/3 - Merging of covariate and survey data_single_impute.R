
wd <- "~path to root directory"
setwd(wd)

library(tidyverse)

marker <- "Stunting"
month = "Mar23" #for appending filename

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

## Reading in covariate and Cleaned (by Age and SE) data.
cov_data <- read_rds("1- Sample Input Data/Mean_Imputed_IHME_Mar2023.rds")

stunt_res <- read_rds(paste0("1- Sample Input Data/Global_",marker,"_SE_Age_clean.rds")) %>% 
  as_tibble() %>% 
  rename("ISO.code" = "ISO3Code") %>% 
  mutate(ISO.code = as.character(ISO.code),
         Point.Estimate = Adj_PointEstimate*100)%>% 
  filter(Age_range == "National") %>% 
  select(-c(Age_range,Adj_PointEstimate)) %>% 
  filter(year > 1989)

### Get the total number of primary data sources before merging.
length(stunt_res$Point.Estimate[!is.na(stunt_res$Point.Estimate)])

### Merge with covariate data
Stunt_data_w_cov <- cov_data %>% full_join(stunt_res) %>% 
  dplyr::select(-c("Country")) %>% 
  dplyr::rename("Country"="location_name") %>% 
  filter(year>1989) %>% 
  arrange(.imp, Country, year, Sex) %>% 
  select(-c(M49, country)) %>% 
  select(ISO.code, Country, year, Sex, everything())


### Ensuring only 1 observation per survey year.
Stunt_data_w_cov <- Stunt_data_w_cov %>% 
  arrange(ISO.code, Country, year, Sex, UNICEFSurveyID) %>% 
  mutate(ISO.code  = as.factor(ISO.code),
    year_ID = case_when(
      is.na(UNICEFSurveyID) ~ year,
      !is.na(UNICEFSurveyID) ~ year+UNICEFSurveyID
    ), 
    UNICEFSurveyID = case_when(
      is.na(UNICEFSurveyID) ~ 0,
      !is.na(UNICEFSurveyID) ~ UNICEFSurveyID
    )) %>% 
  group_by(ISO.code, Country, year, UNICEFSurveyID, year_ID, Sex,.imp) %>%
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  dplyr::select(-"year_ID") 

View(Stunt_data_w_cov)

### Get the total number of primary data sources and countries (should be 206) after merging.
B <- max(Stunt_data_w_cov$.imp)
data_one <- Stunt_data_w_cov
length(unique(data_one$ISO.code))
length(data_one$Point.Estimate[!is.na(data_one$Point.Estimate)])/B

### Load in and edit regional groupings ###
reg_groups <- read_csv("1- Sample Input Data/Region_groupings.csv", show_col_types = FALSE)
reg_groups <- reg_groups %>% dplyr::select(-c(country,Country))
HIC_list <- read_csv("1- Sample Input Data/HIC_list_Nov_2020.csv", show_col_types = FALSE)
HIC_list <- HIC_list %>% dplyr::select(-c(country)) %>% rename(ISO.code = ISO)
reg_groups_HIC <- left_join(reg_groups,HIC_list)
reg_groups_HIC$All_africa_HIW <- as.character(reg_groups_HIC$All_africa_HIW)
unique(reg_groups_HIC$All_africa_HIW[is.na(reg_groups_HIC$HIC2019)])
reg_groups_HIC$All_africa_HIW[!is.na(reg_groups_HIC$HIC2019)] <- "High-income Countries"
reg_groups <- reg_groups_HIC %>% select(c(ISO.code,All_africa_HIW))

### Merge regional groupings with data ###
data_one <- merge(data_one,reg_groups,by="ISO.code",all = TRUE)
data_miss <- data_one[is.na(data_one$All_africa_HIW),]
data_one$All_africa_HIW <- apply(as.matrix(tolower(as.character(data_one$All_africa_HIW))),1,simpleCap)
data_one$country <- as.character(data_one$Country)
### Create "East and Central Africa" region
data_one$All_africa_HIW[data_one$All_africa_HIW=="East Africa"] <- "East and Central Africa"
data_one$All_africa_HIW[data_one$All_africa_HIW=="Central Africa"] <- "East and Central Africa"

### Check data to see if n' are correct
View(data_one)
length(c(data_one$Point.Estimate[!is.na(data_one$Point.Estimate)]))/B
length(unique(data_one$ISO.code))
unique(data_one$All_africa_HIW)

### Doing some cleaning of the outcome data.
if(marker == "Overweight"){
  ## Changing datapoint for India consultation
  data_one$Point.Estimate[data_one$ISO.code=="IND" & data_one$year == 1997] <- NA
  
  
  # DEMOCRATIC PEOPLE'S REP. OF KOREA (THE) MICS 2009 has 0 for overweight.  Here putting in the 1/(2*n).
  data_one$SE_val[data_one$Point.Estimate ==0 &
                    !is.na(data_one$Point.Estimate) ] <- 0.1
  data_one$Point.Estimate[data_one$Point.Estimate ==0 &
                            !is.na(data_one$Point.Estimate) ] <- 
    min(data_one$Point.Estimate[data_one$Point.Estimate>0],na.rm=TRUE)
}

# Changing estimates=0 to the minimum value
min_val_pr <- min(data_one$Point.Estimate[data_one$Point.Estimate>0],na.rm=TRUE)
min_val_SE <- min(data_one$SE_val[data_one$SE_val>0],na.rm=TRUE)
# Standardizing SDI by Sex
# Setting the Region variable
data_one <- as_tibble(data_one) %>% 
  group_by(Sex) %>% 
  mutate(SEV.Z = scale(SEV)) %>% 
  ungroup(Sex) %>% 
  arrange(Country,year) %>% 
  mutate(
    SEV.Z = SEV.Z[,1],
    SE_val = case_when(
      SE_val == 0 ~ min_val_SE/100,
      TRUE ~ SE_val/100
    ),
    Point.Estimate = case_when(
      (Point.Estimate==0 & !is.na(Point.Estimate)) ~ min_val_pr/100,
      TRUE ~ Point.Estimate/100
    )
  ) 




### Output data
if(marker=="Stunting"){
  saveRDS(data_one,file = paste("1- Sample Input Data/Stunt_data_w_cov_",month,"_reg_all.rds",sep=""))
  
  }
if(marker=="Overweight"){
  
  saveRDS(data_one,file = paste("1- Sample Input Data/Over_data_w_cov_",month,"_reg_all.rds",sep=""))
  
  }










