#-----------------------------------------------------------------------------
# This file should take as in input the original raw csv files: 
# RawDataApr2021.csv and Pop_size_Data_WorldBank2.csv
# and perform all processing including the generating of spatial lag values. 
# It should output the final csv file used for model fits.
#-----------------------------------------------------------------------------

library(tidyverse)
#install.packages("factoextra")
library(factoextra) #pca visualization

# load data set with cases per million (replace thos with original generating code)
cvd19_pre <- read_csv("data/cvd19_2021_06_24.csv") %>% 
  rename(cases_per_mil = cases, deaths_per_mil = deaths)

glimpse(cvd19_pre)

####feature engineering#####

 #is mean age/pathogen diversity a better predictor of 
#cases than mean age alone?

PathComp <- cvd19_pre %>% select(country, Lymphatic.filariasis,
                                 Schistosomiasis,Hookworm,Tuberculosis,
                                 Trichuriasis, Malaria, Ascariasis,
                                 Genital.herpes, HIV.AIDS ) %>%
  column_to_rownames( var = "country") %>%
  princomp(cor = TRUE, scores = TRUE)

res.var <- get_pca_var(PathComp)
res.var$coord          # Coordinates
res.var$contrib

PathComp  %>% fviz_pca_biplot( repel = TRUE,col.var = "#2E9FDF", # Variables color
                               col.ind = "#696969")  # Individuals color
                                            
PathComp %>% fviz_eig()

res.ind <-  get_pca_ind(PathComp )

coords_PathComp <- res.ind$coord %>%
         as.data.frame() %>%
        rownames_to_column() %>%
        select(rowname, Dim.1, Dim.2, Dim.3) %>%  # 3 pcs
        rename(country= rowname, PathCompPC1=Dim.1,
                       PathCompPC2=Dim.2, PathCompPC3=Dim.3 )                        

pop <- read_csv("data/Pop_size_Data_WorldBank2.csv") %>% 
  rename(pop = Pop_size, country = Country)

#####put it all together


cvd19_all <- left_join(cvd19_pre, pop, by = "country") %>%
  left_join( coords_PathComp, by = "country") %>%
  left_join( coords_ageDiv, by = "country") %>%
  mutate(cases = round(cases_per_mil*pop*(1e6)^-1,0),
                       deaths = round(deaths_per_mil*pop*(1e6)^-1,0),
                       log_tests= log(tests),
                       log_cases = log(cases),
                       log_pop = log(pop),
                       IDburd_perMil =IDBurd/(pop/1e6),# or just use log_pop?
                       tests = NULL) %>%
  mutate(across(where(is.character), as.factor))


# inspect correlation of predictors
cvd19_all %>% dplyr::select(where(is.numeric)) %>% cor %>% corrplot::corrplot()

# remove strongly correlated predictors
cvd19 <- cvd19_all %>% select(-invSimp, -GDPpc, -IDBurd,
                              -cases_per_mil, -deaths_per_mil, #pathCOmpPC1 stronggly correlated with ageDiv
                              -invSimp, -Lymphatic.filariasis,
                              -Schistosomiasis,-Hookworm,-Tuberculosis,
                              -Trichuriasis, -Malaria, -Ascariasis,
                              -Genital.herpes, -HIV.AIDS, -temp, -ageDivPC1)

#write_csv(cvd19,"data/cvd19_2022_10_17.csv")
