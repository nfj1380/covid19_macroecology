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

#pathogen composition
PathComp_Indiapre <- India_endemic  %>% 
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
  select(rowname, Dim.1, Dim.2) %>%  # 2 pcs
  rename(Province_State= rowname, PathCompPC1=Dim.1,
         PathCompPC2=Dim.2 )

cvd19_India_all <- left_join(india_data, inverse_simpFinal, by = "Province_State") %>%
  left_join( coords_PathComp, by = "Province_State") %>% 
  mutate(log_cases = log(cases),
         log_pop = log(pop),
         IDburd_perMil =IDBurd/(pop/1e6)) %>% 
  mutate(across(where(is.character), as.factor))

# inspect correlation of predictors
cvd19_India_all %>% dplyr::select(where(is.numeric)) %>% cor %>% corrplot::corrplot()

cvd19_India_sub <- cvd19_India_all%>% select(-names(India_endemic),
                                             -PathCompPC1) #strong corr with diversity
                                             



write_csv(cvd19_India_sub,"data/cvd19_India.csv")
