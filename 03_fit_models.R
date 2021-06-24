pacman::p_load("tidyverse","brms")

rm(list=ls())

# load data set
cvd19 <- read_csv("data/cvd19_2021_06_24.csv")

run_date <- "2021_06_24"

# Model specification with reduced set of covariates
# Set k=6 as the maximal spline basis dimensions to reduce 
# model complexity while permitting sufficient non-linearity.
# We use the spline basis "ts" which applies regularisation to the
# linear basis functions in addition to the usual thin-plate 
# smoothing penalty which is applied by default to the non-smooth terms.

f_cases_mu <- cases ~ 1 + reg + 
  trans + 
  s(HIV.AIDS, k = 6, bs = "ts") + 
  s(Ascariasis, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(PerUrb, k = 6, bs = "ts") + 
  s(tests, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts") 

f_deaths_mu <- deaths ~ 1 + reg + 
  s(cases, k = 6, bs = "ts") + 
  s(meanAge, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(Trichuriasis, k = 6, bs = "ts") + 
  s(popDen, k = 6, bs = "ts") + 
  s(Hookworm, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts") 

f_shape <- shape ~ 1 + 1|reg

form_cases <- bf(f_cases_mu,f_shape) + negbinomial()
form_deaths <- bf(f_deaths_mu,f_shape) + negbinomial()

# sampler settings
seed <- 687756 # sample(1e6,1)
control = list(adapt_delta = 0.9)
cores = 4
chains = 4
iter = 8000
thin = 16 # to leave 1000 samples
init_r = 1

# fit and save case models (< 3 mins)
fit_cases <- brm(form_cases, 
                 data = cvd19, 
                 cores = cores, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 file = paste0("results/fit_cases_",run_date),
                 control = control)

# fit and save case models (< 6 mins)
fit_deaths <- brm(form_deaths, 
                 data = cvd19, 
                 cores = cores, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 file = paste0("results/fit_deaths_",run_date),
                 control = control)



