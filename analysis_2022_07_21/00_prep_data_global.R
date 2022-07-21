#-----------------------------------------------------------------------------
# This file should take as in input the original raw csv files: 
# RawDataApr2021.csv and Pop_size_Data_WorldBank2.csv
# and perform all processing including the generating of spatial lag values. 
# It should output the final csv file used for model fits.
#-----------------------------------------------------------------------------

library(tidyverse)

# load data set with cases per million (replace thos with original generating code)
cvd19_pre <- read_csv("data/cvd19_2021_06_24.csv") %>% 
  rename(cases_per_mil = cases, deaths_per_mil = deaths)

pop <- read_csv("data/Pop_size_Data_WorldBank2.csv") %>% 
  rename(pop = Pop_size, country = Country)

cvd19_all <- left_join(cvd19_pre, pop, by = "country") %>% 
  mutate(cases = round(cases_per_mil*pop*(1e6)^-1,0),
         deaths = round(deaths_per_mil*pop*(1e6)^-1,0),
         log_tests= log(tests),
         log_cases = log(cases),
         log_pop = log(pop),
         tests = NULL) %>% 
         mutate(across(where(is.character), as.factor))

# inspect correlation of predictors
cvd19_all %>% dplyr::select(where(is.numeric)) %>% cor %>% corrplot::corrplot()

# remove strongly correlated predictors
cvd19 <- cvd19_all %>% select(-invSimp, -GDPpc, -IDBurd, -cases_per_mil, -deaths_per_mil)

#write_csv(cvd19,"data/cvd19_2022_07_04.csv")
