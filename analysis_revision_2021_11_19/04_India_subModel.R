pacman::p_load("tidyverse","brms", 'sp', 'spdep', 'spatialreg')

rm(list=ls())

# load and prepare data set-------------------------------
cvd19_india <- read_csv("data/India_covid19.csv") %>% 
  rename(cases = Cases_1July2021_cumulative, deaths = Deaths_1July2021_cumulative,
         urban=`%urban` , pop=Population )
  
cvd19_india$Province_State <-str_replace(cvd19_india$Province_State, " ", "_")

cvd19_imputed <- mice::mice(cvd19_india ,m=5,maxit=50,method = "cart",seed=500)
completed_covid_Data <- complete(cvd19_imputed ,1)

cvd19Ind <- completed_covid_Data  %>% 
  mutate(log_tests_pp= log(Tested_1July2021_cumulative)) %>% 
  mutate(across(where(is.character), as.factor)) %>% 
  mutate(log_healthcare = log(Health_expenditure))#needs to be scaled

#India data is too small for missForest

#data_prep <- select(cvd19, -Province_State) %>% 
  #missForest::missForest()

#-----------------------------------------------------------------------------------
#create a lag variable to account for neighbouring cases
#-----------------------------------------------------------------------------------

India_Spatial <- cvd19Ind
coordinates(India_Spatial ) <- ~Lat+Long

nb <- tri2nb(India_Spatial  ) #create neighbour list 
spatL <- nb2listw(nb)
plot(nb, coordinates(India_Spatial ))

cvd19Ind$lag_rate<-lag.listw(x=spatL, var=(cvd19Ind$cases))

#---------------------------------------------------------
#could imputation have an effect?

# complete_cvd19 <- cvd19_india[complete.cases(cvd19_india),]
# 
# cvd19_noNA <- complete_cvd19  %>% 
#   mutate(log_tests_pp= log(Tested_1July2021_cumulative)) %>% 
#   mutate(across(where(is.character), as.factor))
# 
# India_Spatial_noNA <- cvd19_noNA
# coordinates(India_Spatial_noNA) <- ~Lat+Long
# nb <- tri2nb(India_Spatial_noNA  ) #create neighbour list 
# spatL <- nb2listw(nb)
# plot(nb, coordinates(India_Spatial_noNA ))
# cvd19_noNA$lag_rate<-lag.listw(x=spatL, var=(cvd19_noNA$cases))

# Data exploration

#DataExplorer::create_report(cvd19_noNA)

# ggplot(cvd19_noNA, aes(y=log(deaths), x=Ascariasis))+
#   geom_point()
#tv<- lm(cases~Malaria, data=complete_cvd19)
#summary(tv)

#corrS <- cvd19Ind %>% dplyr::select(where(is.numeric)) %>% cor %>% corrplot::corrplot()
# 
# # prepare variables and formula constructors--------------
# vars_numeric <- cvd19 %>% 
#   select(where(is.numeric), -contains("cases"),
#          -contains("deaths"), -Lat, -Long, -Tested_1July2021_cumulative) %>% names; vars_numeric
# # 

# sampler settings----------------------------------------
seed <- 768021 # sample(1e6,1)
control = list(adapt_delta = 0.9)
cores = 4
chains = 4
iter = 20000
thin = 10 # to leave 1000 samples
init_r = 0
run_date <- "2021_12_07_Final"
#---------------------------------------------------------
############CASE MODEL############

f_cases_mu <- cases|rate(pop) ~ 1+
  #s(HIV, k = 6, bs = "tp") + prevalence too low (< 1%)
  s(Ascariasis, k = 6, bs = "tp") + 
  s(urban, k = 6, bs = "tp") + 
  s(log_tests_pp, k = 6, bs = "tp") + 
  s(lag_rate, k = 6, bs = "tp") + 
  s(Mean_age, k = 6, bs = "tp")

f_shape <- shape ~ 1
#not enough data for mean age - missing 15 states

form_cases <- bf(f_cases_mu, f_shape) + negbinomial()


# fit and save case models (< 3 mins)
fit_cases <- brm(form_cases, 
                 data = cvd19_noNA, 
                 cores = cores, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 file = paste0("results/fit_cases_India",run_date),
                 control = control)

summary(fit_cases)

resp <- "cases"

plot_ce(fit_cases, "Ascariasis")


#---------------------------------------------------------
############DEATH MODEL############

f_deaths_mu <- deaths|rate(pop) ~ 1+
  cases+ #doesn't seem to like this
 # s(Malaria, k = 6, bs = "tp") + 
  s(Trichuriasis , k = 6, bs = "tp") + 
  #s(Hookworm, k = 6, bs = "tp") + 
  s(log_tests_pp, k = 6, bs = "tp") + 
  s(lag_rate, k = 6, bs = "tp")  + 
 # s(log_healthcare , k = 6, bs = "tp")+ 
  s(Mean_age , k = 6, bs = "tp")


f_shape <- shape ~ 1
#not enough data for mean age - missing 15 states

form_deaths <- bf(f_deaths_mu, f_shape) + negbinomial()

run_date <- "2021_12_08_test3"

# fit and save case models (< 3 mins)
fit_deaths <- brm(form_deaths, 
                 data = cvd19Ind, 
                 cores = cores, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 file = paste0("results/fit_deaths_India",run_date),
                 control = control)

summary(fit_deaths)

resp <- "deaths"

plot_ce(fit_deaths, "Malaria")
plot_ce(fit_deaths, "Hookworm")
plot_ce(fit_deaths, "Trichuriasis")

