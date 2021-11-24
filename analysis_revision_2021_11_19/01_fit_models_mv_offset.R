#-------------------------------------
# Model specification with reduced set of covariates
# Set k=6 as the maximal spline basis dimensions to reduce 
# model complexity while permitting sufficient non-linearity.
# We use the spline basis "ts" which applies regularisation to the
# linear basis functions in addition to the usual thin-plate 
# smoothing penalty which is applied by default to the non-smooth terms.
#-------------------------------------

pacman::p_load("tidyverse","brms")

rm(list=ls())

# load and prepare data set-------------------------------
cvd19_pre <- read_csv("data/cvd19_2021_06_24.csv") %>% 
  rename(cases_per_mil = cases, deaths_per_mil = deaths)
pop <- read_csv("data/Pop_size_Data_WorldBank2.csv") %>% 
  rename(pop = Pop_size, country = Country)

cvd19_pre$tests

cvd19 <- left_join(cvd19_pre, pop, by = "country") %>% 
  mutate(cases = round(cases_per_mil*pop*(1e6)^-1,0),
         deaths = round(deaths_per_mil*pop*(1e6)^-1,0),
         log_tests_pp= log(tests/pop),
         log_prop_cases = log(cases_per_mil/1e6),
         tests = NULL) %>% 
         mutate(across(where(is.character), as.factor))

#---------------------------------------------------------

#cvd19 %>% dplyr::select(where(is.numeric)) %>% cor %>% corrplot::corrplot()

# prepare variables and formula constructors--------------
vars_numeric <- cvd19 %>% 
  select(where(is.numeric), -contains("cases"),
         -contains("deaths"), -pop, -invSimp, -GDPpc, -IDBurd) %>% names; vars_numeric
vars_factor <- cvd19 %>% select(where(is.factor), -country) %>% names; vars_factor

# regional dependence: 
mu_spline <- function(response, vars_num, vars_fac, k = 6, basis = "tp"){
  paste(response," | rate(pop) ~ 1 + ", paste(vars_fac, collapse = " + "), " + ", 
        paste0("s(",vars_num,", k = ", k,", bs = '", basis, "')", collapse = " + ")) %>% 
    as.formula()
} 

mu_linear <- function(response,vars) paste(response," | rate(pop) ~ 1 + ",paste(vars, collapse = " + ")) %>% as.formula()

#---------------------------------------------------------

# sampler settings----------------------------------------
seed <- 768021 # sample(1e6,1)
control = list(adapt_delta = 0.9)
cores = 4
chains = 4
iter = 10000
thin = 10 # to leave 1000 samples
init_r = 0
run_date <- "2021_11_24"
#---------------------------------------------------------

# fit all vars to select best predictors: cases/deaths----
if(T){
  resp  <- "deaths" # cases/deaths
  type <- "spline" # "spline" or "linear"
  message(paste0("Fitting model: ",resp))
  if(resp=="deaths") vars_numeric <- c("log_prop_cases",vars_numeric)
  f_mu <- mu_spline(resp,vars_numeric, vars_factor, k = 6, basis = "ts")
  if(type == "linear") f_mu <- mu_linear(resp, c(vars_factor, vars_numeric))
  f_shape <- shape ~ 1 + 1|reg
  f_model <- bf(f_mu,f_shape) + negbinomial(); f_model
  f_mu
  brm(f_model, 
      data = cvd19, 
      cores = cores, 
      chains = chains,
      init_r = init_r, 
      iter = iter, 
      thin = thin, 
      seed = seed, 
      file = paste0("results/fit_vars_all_",resp,"_",type,"_",run_date),
      control = control)
  }

# fit cases submodel-------------------------------------
if(F){
  f_cases_mu <- cases | rate(pop) ~ 
    1 + 
    reg + 
    trans + 
    s(HIV.AIDS, k = 6, bs = "ts") + 
    s(Ascariasis, k = 6, bs = "ts") + 
    s(Malaria, k = 6, bs = "ts") + 
    s(PerUrb, k = 6, bs = "ts") + 
    s(log_tests, k = 6, bs = "ts") + 
    s(spatLag, k = 6, bs = "ts") 
  
  f_shape <- shape ~ 1 + 1|reg
  
  form_cases <- bf(f_cases_mu,f_shape) + negbinomial(); form_cases
  
  # fit and save case models (< 3 mins)
  fit_cases <- brm(form_cases, 
                   data = cvd19, 
                   cores = cores, 
                   chains = chains,
                   init_r = init_r, 
                   iter = iter, 
                   thin = thin, 
                   seed = seed, 
                   file = paste0("results/fit_cases_nb_sreg",run_date),
                   control = control)
  
  fit_nb_sreg <- readRDS("results/fit_cases_nb_sreg2021_11_22.rds")
  #fit_nb_s1 <- readRDS("results/fit_cases_nb_s12021_11_22.rds")
  #fit_pois <- readRDS("results/fit_cases_pois2021_11_22.rds")
  
  loo_fits <- list(pois = fit_pois, nb_reg = fit_nb_sreg, nb_s1 = fit_nb_s1) %>% 
    map(loo, cores = 20) %>% loo_compare()
  
  #model    elpd_diff  se_diff 
  #---------------------------
  #nb_reg        0.0        0.0   # nb a clear winner, good support for shape ~ 1 + 1|reg
  #nb_s1       -33.3       22.4
  #pois   -6501154.9   842676.1
}

# fit deaths submodel------------------------------------
if(F){
  
  # new shortlist
  # meanAge, 
  
  f_deaths_mu <- deaths|rate(pop) ~ 
    1 + 
    reg + 
    s(cases, k = 6, bs = "ts") + 
    s(meanAge, k = 6, bs = "ts") + 
    s(Malaria, k = 6, bs = "ts") + 
    s(Trichuriasis, k = 6, bs = "ts") + 
    s(popDen, k = 6, bs = "ts") + 
    s(Hookworm, k = 6, bs = "ts") + 
    s(spatLag, k = 6, bs = "ts") 
    
    f_shape <- shape ~ 1 + 1|reg
    
    form_deaths <- bf(f_deaths_mu, f_shape) + negbinomial()
    
# fit and save case models (< 6 mins)
fit_deaths <- brm(form_deaths, 
                 data = cvd19, 
                 cores = cores, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 file = paste0("results/fit_deaths_nb_reg",run_date),
                 control = control)

fit_deaths <- readRDS("results/fit_deaths_nb_reg2021_11_22.rds") 
fit_deaths %>% conditional_effects()
}
