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


# sampler settings----------------------------------------
seed <- 768021 # sample(1e6,1)
control = list(adapt_delta = 0.95)
cores = 1
chains = 4
iter = 20000
thin = 0.5*iter*chains/1000 # to leave 1000 samples
init_r = 0
run_date <- "2022_03_21"
#---------------------------------------------------------


#---------------------------------------------------------
# cases submodel
#---------------------------------------------------------

f_cases_mu <- cases|rate(pop) ~ 1+
  #s(HIV, k = 6, bs = "tp") + prevalence too low (< 1%)
    s(Ascariasis, k = 5, bs = "tp") + 
  s(urban, k = 5, bs = "tp") + 
  s(log_tests_pp, k = 5, bs = "tp") + 
  s(lag_rate, k = 5, bs = "tp") + 
  s(Mean_age, k = 5, bs = "tp")


f_shape <- shape ~ 1

form_cases <- bf(f_cases_mu, f_shape) + negbinomial()
form_cases
future::plan("multisession", workers = 20)

# fit and save case models
fit_cases <- brm_multiple(form_cases, 
                 data = cvd19_imp, # list of imputed data sets
                 combine = FALSE, # combine posteriors
                 future = TRUE, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 #file = paste0("results/fit_cases_India_imp_",run_date),
                 control = control)

#saveRDS(fit_cases,"results/old_fits/fit_cases_mi_separate_2022_03_21.rds")

cases_comb <- combine_models(mlist = fit_cases, check_data = F)
resp <- "cases"

cases_comb %>% conditional_effects("Ascariasis")

vars_cases <- c("Ascariasis","urban","log_tests_pp","lag_rate","Mean_age")

resp <- "cases"
for(x in c("Ascariasis","urban","log_tests_pp","lag_rate","Mean_age")){
  fit_cases %>% 
    map(~.x %>% plot_ce(x)) %>% 
    ggpubr::ggarrange(plotlist = .) %>% 
    ggsave(paste0("plots/ce_mi_",x,".pdf"), plot = .)
}

c("Ascariasis","urban","log_tests_pp","lag_rate","Mean_age") %>% 
  map(~ plot_ce(cases_comb, .x)) %>% ggpubr::ggarrange(plotlist = .) %>% 
  ggsave("plots/ce_mi_combined.pdf",plot = .)


#---------------------------------------------------------
# deaths submodel
#---------------------------------------------------------

run_date <- "2022_03_30"

f_deaths_mu <- deaths|rate(pop) ~ 1+ 
  s(log_cases , k = 4, bs = "tp") + 
  s(Malaria, k = 4, bs = "tp") + 
  s(Trichuriasis , k = 4, bs = "tp") + 
  # s(Hookworm, k = 6, bs = "tp") + # prevalence too low
  #s(log_tests_pp, k = 6, bs = "tp") + #tests not impacting deaths here 
  s(lag_rate, k = 4, bs = "tp")  +  #doesn't look to be important
  s(log_healthcare , k = 4, bs = "tp")+ 
  s(Mean_age , k = 4, bs = "tp")


f_shape <- shape ~ 1

form_deaths <- bf(f_deaths_mu, f_shape) + negbinomial()

future::plan("multisession", workers = 20)
fit_deaths <- brm_multiple(form_deaths, 
                          data = cvd19_imp, # list of imputed data sets
                          combine = FALSE, # combine posteriors
                          future = TRUE, 
                          chains = chains,
                          init_r = init_r, 
                          iter = iter, 
                          thin = thin, 
                          seed = seed, 
                          control = control)

#saveRDS(fit_deaths,paste0("results/old_fits/fit_deaths_mi_separate_",run_date,".rds"))

# deaths plots
run_date <- "2022_03_30"
#fit_deaths <- readRDS(paste0("results/old_fits/fit_deaths_mi_separate_",run_date,".rds"))
deaths_comb <- combine_models(mlist = fit_deaths, check_data = F)

vars_deaths <- c("log_cases","Malaria", "Trichuriasis","lag_rate","log_healthcare","Mean_age")
resp <- "deaths"

# plots for all predictors for all ten imputed data sets
for(x in vars_deaths){
  fit_deaths %>% 
    map(~.x %>% plot_ce(x, exp = if_else(x == "log_cases",TRUE, FALSE))) %>% 
    ggpubr::ggarrange(plotlist = .) %>% 
    ggsave(paste0("plots/deaths_ce_mi_",x,"_",run_date,".pdf"), plot = .)
}

# plots for all predictors for the combined model
vars_deaths %>% 
  map(~ plot_ce(deaths_comb, .x, exp = if_else(.x == "log_cases",TRUE, FALSE))) %>% 
  ggpubr::ggarrange(plotlist = .) %>% 
  ggsave(paste0("plots/deaths_ce_mi_combined_",run_date,".pdf"),plot = .)

# plot for Malaria for combined model
run_date <- "2022_03_30"
resp <- "deaths"
deaths_comb <- readRDS(paste0("results/fits_deaths_mi_combined_",run_date))
deaths_comb %>% plot_ce("Malaria", exp = F)
ggsave(paste0("plots/deaths_India_Malaria_",run_date,".pdf"))
