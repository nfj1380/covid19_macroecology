# prep Indian data set

pacman::p_load("tidyverse","brms", 'sp', 'spdep', 'spatialreg')

rm(list=ls())

source("00_functions.R")

theme_set(theme_classic())

# load and prepare data set-------------------------------
cvd19_india <- read_csv("data/India_covid19.csv") %>% 
  rename(cases = Cases_1July2021_cumulative, 
         deaths = Deaths_1July2021_cumulative,
         urban=`%urban`,
         pop=Population ) %>% 
  mutate(Province_State = str_replace(Province_State, " ", "_"))

#cvd19_india$Province_State <-str_replace(cvd19_india$Province_State, " ", "_")

# impute 10 data sets
cvd19_imputed <- mice::mice(cvd19_india ,m=10,maxit=50,method = "cart",seed=500)

# transform data and add spatial lag
cvd19_imp <- 
  mice::complete(cvd19_imputed ,"all") %>% 
  map(~ .x %>%# transform variables
        mutate(log_tests_pp= log(Tested_1July2021_cumulative)) %>% 
        mutate(across(where(is.character), as.factor)) %>% 
        mutate(log_healthcare = log(Health_expenditure),
               log_cases = log(cases))) %>% 
  map(function(df){ # compute and add spatial lag
    df_spat <- df
    coordinates(df_spat) <- ~Lat+Long
    nb <- tri2nb(df_spat)
    spatL <- nb2listw(nb)
    df$lag_rate<-lag.listw(x=spatL, var=(df$cases))
    return(df)
  })

#save one iteration

india_data <- as.data.frame(cvd19_imp[1])
names(india_data) <- gsub(pattern = "X1.", replacement = "", 
                        x = names(india_data))

write_csv(india_data,"data/cvd19_India.csv")
