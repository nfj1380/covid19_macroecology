#-----------------------------------------------------------------------------
# This file should take as in input the original raw csv files: 
# RawDataApr2021.csv and Pop_size_Data_WorldBank2.csv
# and perform all processing including the generating of spatial lag values. 
# It should output the final csv file used for model fits.
#-----------------------------------------------------------------------------

pacman::p_load("tidyverse","factoextra", 'sp', 'spdep', 'spatialreg', 'mice')


# load data set with cases per million (replace thos with original generating code)
cvd19_pre <- read_csv("data/RawDataOct2023.csv") %>% 
  rename(cases_per_mil = cases_perMil, deaths_per_mil = deaths_perMil,
         country = Country, tests =Total_tests_per_thousand,
         cardioDR = Cardiovasc_death_rate,
         reg =WHO_Region)  #countries with no cases are isolated islands
#lots of missing data

#convert region names
cvd19_pre$reg <- recode_factor(cvd19_pre$reg,  'Eastern Mediterranean'='EMed',
                                Africa= 'Afri', Europe = 'Euro',
                               'South-East Asia'='SEAs',Americas= 'Amer',
                               'Western Pacific' = 'WPac')
#convert trans 
cvd19_pre$Transmission_Classification <- recode_factor(cvd19_pre$Transmission_Classification,
                                                       Clusters_of_cases = 'clust',
                                                       Community_transmission = 'comm',
                                                       Sporadic_cases='spor',
                                                       Pending = 'spor') #pending is similar to sporadic

cvd19_imputed <- mice::mice(cvd19_pre ,m=10,maxit=50,method = "cart",seed=500)

cvd19_imp <- mice::complete(cvd19_imputed ,"all") %>% 
  map(function(df){ # compute and add spatial lag
    df_spat <- df
    coordinates(df_spat) <- ~latitude+longitude
    nb <- tri2nb(df_spat)
    spatL <- nb2listw(nb)
    df$lag_rate<-lag.listw(x=spatL, var=(df$pscore.mean)) #cases or mortality?
    return(df)
  })

####feature engineering#####

#path data doesn't have any NAs

vars_endemic_indiv <- c("Lymphatic_filariasis",
                  'Schistosomiasis','Hookworm_disease','Tuberculosis',
                  'Trichuriasis', 'Malaria', 'Ascariasis',
                  'Genital.herpes', 'HIV.AIDS')

vars_endemic_nokeyPath <- c("Lymphatic_filariasis",
                        'Schistosomiasis','Hookworm_disease','Tuberculosis',
                        'Trichuriasis', 'Ascariasis',
                        'Genital.herpes')
                  
PathComp <- cvd19_imp[[1]] %>% select(country, vars_endemic_nokeyPath) %>%
  column_to_rownames( var = "country") %>%
  princomp(cor = TRUE, scores = TRUE)


# res.var <- get_pca_var(PathComp)
# res.var$coord          # Coordinates
# res.var$contrib

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


#Helminth combined prevalence. Not used in the end as strongly correlated with PCA2
# vars_worm <- c("Lymphatic_filariasis",
#                'Schistosomiasis','Hookworm_disease',
#                'Trichuriasis', 'Ascariasis')
# 
# worm_prev <- cvd19_imp[[1]] %>% select(country, vars_worm ) %>% 
#             mutate(HelminthPrev = rowSums(across(where(is.numeric)))) %>% 
#           select(country, HelminthPrev)


#prev <- cvd19_imp %>% select(country, vars_hiv_malaria )

pop <- read_csv("data/Pop_size_Data_WorldBank2.csv") %>% 
  rename(pop = Pop_size, country = Country)


#####put it all together


cvd19_all <- left_join(cvd19_imp[[1]], pop, by = "country") %>%
  left_join( coords_PathComp, by = "country") %>%
  left_join( worm_prev, by = "country") %>% 
  #left_join( prev, by = "country") %>% 
  mutate(cases = round(cases_per_mil*pop*(1e6)^-1,0),
                       deaths = round(deaths_per_mil*pop*(1e6)^-1,0),
                       log_tests= log(tests),
                       log_cases = log(cases),
                       log_pop = log(pop),
                       IDburd_perMil =inf_diseaseBurden/(pop/1e6),
                       Vaccinated_perMil = People_vaccinated/(pop/1e6),# or just use log_pop?
                       tests = NULL) %>%
  mutate(across(where(is.character), as.factor)) %>% 
  select(-HCexpend) #perfectly collinear with GDP


# inspect correlation of predictors
cvd19_all %>% dplyr::select(where(is.numeric),-vars_endemic_nokeyPath,
      -contains(c("cases","death",'excess',"mortality","log","pscore", "Rainfall", "Temp", "Vaccinated", 'vaccinated', 'Helminth', 'inv_')), -latitude, -longitude ) %>%
              cor %>% corrplot::corrplot()
#0.7 threshold
#id_burd_per_mill and pathogenDiv_PC1/ivSimp strongly correlated.
# remove strongly correlated predictors
cvd19 <- cvd19_all %>% select(-inv_simpsons , -mean_age, -inf_diseaseBurden,
                              -cases_per_mil, -deaths_per_mil, 
                               -vars_endemic_nokeyPath,  -Ascariasis,
                               -contains('.y'), -HelminthPrev, -People_vaccinated)
#glimpse(cvd19)
cvd19final <- cvd19 %>%  rename( temp = Avg_temp91.2016_c,  
                           rainfall = Avg_Rainfall91.2016_mm, trans = Transmission_Classification,
                           density = popDen)
glimpse(cvd19final)

#write_csv(cvd19final,"data/cvd19_2023_04_19.csv")
