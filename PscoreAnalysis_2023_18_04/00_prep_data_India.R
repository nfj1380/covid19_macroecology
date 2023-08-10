# prep Indian data set

pacman::p_load("vegan", "factoextra","tidyverse","brms", 'sp', 'spdep', 'spatialreg')

rm(list=ls())

source("00_functions.R")

theme_set(theme_classic())

# load and prepare data set-------------------------------
cvd19_india <- read_csv("data/India_covid19_2021_08_30.csv") %>% 
  rename(cases = Cases_1July2021_cumulative, 
         deaths = Deaths_1July2021_cumulative,
         pop=Population,
         IDBurd = InfDiseaseBurden) %>% 
  mutate(Province_State = str_replace(Province_State, " ", "_"))

#cvd19_india$Province_State <-str_replace(cvd19_india$Province_State, " ", "_")

# impute 10 data sets
cvd19_imputed <- mice::mice(cvd19_india ,m=10,maxit=50,method = "cart",seed=500)

# transform data and add spatial lag
cvd19_imp <- 
  mice::complete(cvd19_imputed ,"all") %>% 
  map(~ .x %>%# transform variables
        mutate(log_tests= log(Tested_1July2021_cumulative)) %>% 
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


#add endmeic disease variables

India_endemic <- india_data %>% select(Province_State, 
                                           Hookworm,Tubercluosis,
                                           Trichuriasis, Malaria, Ascariasis, 
                                           Herpes, HIV )  %>% 
                column_to_rownames( var = "Province_State")

#pathogen diversity
inverse_simp <- diversity(India_endemic, index = "invsimpson") %>% 
              as.data.frame()# %>% 
              

names(inverse_simp) <- 'invSimp'

inverse_simpFinal <- inverse_simp %>% 
  rownames_to_column() %>% 
  rename(Province_State= rowname)


vars_endemic_nokeyPath <- c('Hookworm','Tubercluosis',
                            'Trichuriasis','Ascariasis', 
                            'Herpes') 
#glimpse(india_data)
India_endemic_nkp <- india_data %>% select(Province_State,
                                           vars_endemic_nokeyPath) %>% 
   column_to_rownames( var = "Province_State")

#pathogen composition
PathComp_Indiapre <- India_endemic_nkp %>% 
  princomp(cor = TRUE, scores = TRUE) 

PathComp_Indiapre %>% fviz_eig() #quite a few pcs

res.var <- get_pca_var(PathComp_Indiapre)
res.var$coord          # Coordinates
res.var$contrib #first two should suffice (3 has similar loadings to 2)

PathComp_Indiapre %>% fviz_pca_biplot( repel = TRUE,
                               col.var = "#2E9FDF", # Variables color
                               col.ind = "#696969"  # Individuals color
)

res.ind <-  get_pca_ind(PathComp_Indiapre) 

coords_PathComp <- res.ind$coord %>% 
  as.data.frame() %>%
  rownames_to_column() %>% 
  select(rowname, Dim.1, Dim.2, Dim.3) %>%  # 2 pcs
  rename(Province_State= rowname, PathCompPC1=Dim.1,
         PathCompPC2=Dim.2, PathCompPC3=Dim.3 )

vars_worm <- c('Hookworm',
               'Trichuriasis', 'Ascariasis')

worm_prev <- india_data %>% select(Province_State, vars_worm ) %>% 
  mutate(HelminthPrev = rowSums(across(where(is.numeric)))) %>% 
  select(Province_State, HelminthPrev)

glimpse(india_data)

cvd19_India_all <- left_join(india_data, inverse_simpFinal, by = "Province_State") %>%
  left_join( coords_PathComp, by = "Province_State") %>%
  left_join( worm_prev , by = "Province_State") %>% 
  mutate(log_cases = log(cases),
         log_pop = log(pop),
         IDburd_perMil =IDBurd/(pop/1e6),
         health_spending_perMil = Health_expenditure/(pop/1e6)) %>% 
  mutate(across(where(is.character), as.factor))

# inspect correlation of predictors
cvd19_India_all %>% dplyr::select(where(is.numeric), -vars_endemic_nokeyPath, -Lat, -Long,-Health_expenditure,
                                  -contains(c("cases","death","log")) ) %>% 
                                  cor %>% corrplot::corrplot()

glimpse(cvd19_India_all)
cvd19_India_sub <- cvd19_India_all%>% select(-vars_endemic_nokeyPath,-PathCompPC1,
                                             -IDBurd,-IDburd_perMil , 
                                             -invSimp, -Diabetes, -Density,
                                             -Diabetes)  #correlated vars) #strong corr 
                                             



#write_csv(cvd19_India_sub,"data/cvd19_India1.csv")
