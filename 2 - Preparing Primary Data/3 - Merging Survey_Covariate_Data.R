

library(this.path)
wd <- dirname(this.path::here())
print(wd)
setwd(wd)

library(tidyverse)
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

# marker can be "Stunting","Overweight"
marker <- "Stunting"

month <- "May" #for appending filename
year <- "2024"

#### Note: "Imputation of IHME covariate data_multiple_impute.R" or
#### "Imputation of IHME covariate data_single_impute.R" must 
#### be ran before you can commence (as appropriate).

### If doing for multiple imputation set as TRUE, otherwise FALSE
multiple_imputation_vec <- c(TRUE,FALSE)

for(multiple_imputation in multiple_imputation_vec){
  if(multiple_imputation){
    cov_data <- readRDS("Data/IHME_covs/Multiple_Imputed_Mar2023.rds")
  }else{
    cov_data <- readRDS("Data/IHME_covs/Single_Impute_Mar2023.rds") %>% 
      mutate(.imp = 0)
  }
  
  #### Keeping only distinct rows ####
  cov_data_dis <- cov_data %>% 
    group_by(Sex,year,ISO.code,.imp) %>% 
    filter(row_number()==1) %>% 
    ungroup() %>% 
    dplyr::select(-".imp") %>% 
    distinct(.keep_all = TRUE) 
  
  num <- cov_data_dis %>% 
    group_by(Sex,year,ISO.code) %>% 
    summarise(n = n())
  cov_data <- cov_data_dis %>% 
    left_join(num) %>% 
    group_by(Sex,year,ISO.code) %>% 
    mutate(.imp = row_number() - 1*I(n == 1)) %>% 
    ungroup() %>% 
    select(-c("n"))
  
  ## Reading in covariate and Cleaned (by Age and SE) data.
  stunt_res <- read_rds(
    paste0("Data/JME/",year,"/Cleaned/",marker,"_SE_Age_clean.rds")
  ) %>% 
    as_tibble() %>% 
    rename("ISO.code" = "ISO3Code") %>% 
    mutate(ISO.code = as.character(ISO.code),
           Point.Estimate = Adj_PointEstimate*100)%>% 
    filter(Age_range == "National") %>% 
    dplyr::select(-c(Age_range,Adj_PointEstimate)) %>% 
    filter(year > 1989)
  
  ### Get the total number of primary data sources before merging.
  length(stunt_res$Point.Estimate[!is.na(stunt_res$Point.Estimate)])
  
  ### Merge with covariate data
  Stunt_data_w_cov <- cov_data %>% 
    full_join(stunt_res, relationship = "many-to-many") %>% 
    dplyr::select(-c("Country")) %>% 
    dplyr::rename("Country"="location_name") %>% 
    filter(year>1989) %>% 
    arrange(.imp, Country, year, Sex) %>% 
    dplyr::select(-c(M49, country)) %>% 
    dplyr::select(ISO.code, Country, year, Sex, everything())
  
  
  ### Ensuring only 1 observation per survey/year/sex/.imp
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
  
  #View(Stunt_data_w_cov)
  
  ### Get the total number of primary data sources and countries (should be 206) after merging.
  B <- max(Stunt_data_w_cov$.imp)
  length(unique(Stunt_data_w_cov$ISO.code))
  
  
  
  ### Load in and edit regional groupings ###
  reg_groups <- read_xlsx("Data/Country/Crosswalk_Jan_2025.xlsx") %>% 
    select(-"CountryorArea") %>% 
    rename(
      ISO.code = MNF_Country_Code,
      All_africa_HIW = Modelling_Group
    )
  
  ### Merge regional groupings with data ###
  Stunt_data_w_cov <- merge(Stunt_data_w_cov,reg_groups,by="ISO.code",all = TRUE)
  data_miss <- Stunt_data_w_cov[is.na(Stunt_data_w_cov$All_africa_HIW),]
  Stunt_data_w_cov$All_africa_HIW <- apply(as.matrix(tolower(as.character(Stunt_data_w_cov$All_africa_HIW))),1,simpleCap)
  Stunt_data_w_cov$country <- as.character(Stunt_data_w_cov$Country)
  ### Shorten "Eastern Asia, South Eastern Asia And Oceania Excluding Australia & New Zealand" region
  Stunt_data_w_cov$All_africa_HIW[Stunt_data_w_cov$All_africa_HIW=="Eastern Asia, South Eastern Asia And Oceania Excluding Australia & New Zealand"] <- "Eastern Asia, South Eastern Asia And Oceania"
  
  
  ### Check data to see if n' are correct
  length(unique(Stunt_data_w_cov$ISO.code))
  unique(Stunt_data_w_cov$All_africa_HIW)
  
  ### Doing some cleaning of the outcome data.
  if(marker == "Overweight"){
    ## Changing datapoint for India consultation
    Stunt_data_w_cov$Point.Estimate[Stunt_data_w_cov$ISO.code=="IND" & Stunt_data_w_cov$year == 1997] <- NA
    
    
    # DEMOCRATIC PEOPLE'S REP. OF KOREA (THE) MICS 2009 has 0 for overweight.  Here putting in the 1/(2*n).
    Stunt_data_w_cov$SE_val[Stunt_data_w_cov$Point.Estimate ==0 &
                              !is.na(Stunt_data_w_cov$Point.Estimate) ] <- 0.1
    Stunt_data_w_cov$Point.Estimate[Stunt_data_w_cov$Point.Estimate ==0 &
                                      !is.na(Stunt_data_w_cov$Point.Estimate) ] <- 
      min(Stunt_data_w_cov$Point.Estimate[Stunt_data_w_cov$Point.Estimate>0],na.rm=TRUE)
  }
  
  # Changing estimates=0 to the minimum value
  min_val_pr <- min(Stunt_data_w_cov$Point.Estimate[Stunt_data_w_cov$Point.Estimate>0],na.rm=TRUE)
  min_val_SE <- min(Stunt_data_w_cov$SE_val[Stunt_data_w_cov$SE_val>0],na.rm=TRUE)
  # Standardizing SDI by Sex
  # Setting the Region variable
  Stunt_data_w_cov <- as_tibble(Stunt_data_w_cov) %>% 
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
  
  
  # Exporting data
  if(multiple_imputation){
    saveRDS(Stunt_data_w_cov,file = paste0("Data/Merged/",year,"/",marker,"_",month,"_reg_all_multiple_impute.rds"))
  }else{
    saveRDS(Stunt_data_w_cov,file = paste0("Data/Merged/",year,"/",marker,"_",month,"_reg_all_mean_multiple.rds"))
  }
  
}



