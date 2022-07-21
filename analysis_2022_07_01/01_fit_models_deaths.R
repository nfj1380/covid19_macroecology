#-------------------------------------
# Model specification with reduced set of covariates
# Set k=6 as the maximal spline basis dimensions to reduce 
# model complexity while permitting sufficient non-linearity.
# We use the spline basis "ts" which applies regularisation to the
# linear basis functions in addition to the usual thin-plate 
# smoothing penalty which is applied by default to the non-smooth terms.
#-------------------------------------

library(tidyverse)
library(brms)
library(pbmcapply)
library(ggpubr)

rm(list=ls())

theme_set(theme_classic())
source("99_functions.R")

# load data set and identify predictors
cvd19 <- read_csv("data/cvd19_2022_07_01.csv")

vars_numeric <- cvd19 %>% 
  select(where(is.numeric), -deaths, -cases, - spatLag, -log_pop ) %>% names; vars_numeric
vars_factor <- cvd19 %>% select(!where(is.numeric), -country, -reg) %>% names; vars_factor

# set sampler settings
seed <- 768021 # sample(1e6,1)
control = list(adapt_delta = 0.9)
cores = 4
chains = 4
iter = 2000
thin = max(chains*iter/2000,1) # save a maximum of 1000 samples
init_r = 0
run_date <- "2022_07_05"

# fit global model and select best predictors via regularised estimates
if(T){
  resp  <- "deaths" 
  f_mu <- mu_spline(resp,vars_numeric, vars_factor, k = 6, basis = "ts")
  f_shape <- shape ~ 1 + (1|a|reg)
  f_model <- bf(f_mu,f_shape) + negbinomial(); f_model

  brm(f_model, 
      data = cvd19, 
      cores = cores, 
      chains = chains,
      init_r = init_r, 
      iter = iter, 
      thin = thin, 
      seed = seed, 
      file = paste0("results/fit_all_",resp,"_",run_date),
      file_refit = "always",
      control = control)
  }

# manually inspect fits
if(F){
  fit_all_deaths <- readRDS("results/fit_all_deaths_2022_07_01.rds")
  
  # check convergence (Rhat<1.1): very good
  rhat(fit_all_deaths) %>% max # 1.01
  
  ce_plots <- pbmclapply(vars_numeric, function(X) plot_ce(fit_all_deaths,X,resp = "deaths"), mc.cores = 5)
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 5, ncol = 5) %>% 
    annotate_figure(top = text_grob("deaths", size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  ggsave("plots/ce_all_deaths_2022_07_01.pdf", width = 9, height = 9, device = cairo_pdf)
  
}

# fit cases submodel
if(F){
  resp <- "cases"
  vars_submodel <- 
    c("meanAge", "HIV.AIDS", "Ascariasis", 
      "Malaria", "PerUrb", "HCexpend", "cardioDR",
      "spatLag", "log_tests", "log_pop")
  
  f_cases_mu <- mu_spline(resp,vars_submodel, vars_factor, k = 6, basis = "tp")
  f_shape <- shape ~ 1 + (1|a|reg)
  form_cases <- bf(f_cases_mu,f_shape) + negbinomial() 
  
  # try different values of the scale for the sds prior: 1,2.6(default),5,10 
  pr <- set_prior("student_t(3, 0, 2.6)", class = "sds")
   
  # fit and save case models (< 3 mins)
  fit_cases <- brm(form_cases, 
                   data = cvd19, 
                   prior = pr,
                   cores = cores, 
                   chains = chains,
                   init_r = init_r, 
                   iter = iter, 
                   thin = thin, 
                   seed = seed, 
                   file = paste0("results/fit_submodel_",resp,"_",run_date),
                   control = control)
}

## loo for cases models
if(F){
  fit_all_cases <- readRDS("results/fit_all_cases_sreg_2022_07_01.rds")
  fit_sub_cases <- readRDS("results/fit_submodel_cases_sreg_2022_07_01.rds")
  
  future::plan("multisession", workers = 11)
  loo_all <- loo(fit_all_cases, future = T) # 11 problematic obs
  loo_sub <- loo(fit_sub_cases, future = T) # 8 problematic obs
  loo_compare(loo_all,loo_sub)
  
  #m.all_cases_sreg %>% loo(reloo = T, future = T) %>% saveRDS("results/loo_all_cases_sreg.rds") # 40 refits
  #m.all_cases_s1 %>% loo(reloo = T, future = T) %>% saveRDS("results/loo_all_cases_s1.rds") # 41 refits
  #m.sub_cases_sreg %>% loo(reloo = T, future = T) %>% saveRDS("results/loo_subset1_cases_sreg.rds") # 20 refits
  list(all_sreg = readRDS("results/loo_all_cases_sreg.rds"), 
       all_s1 = readRDS("results/loo_all_cases_s1.rds"), 
       sub_sreg = readRDS("results/loo_subset1_cases_sreg.rds")) %>%  
  loo_compare()
  
  #           elpd_diff se_diff
  # sub_sreg    0.0       0.0 
  # all_s1    -62.1      30.9 
  # all_sreg -127.2      73.1 
  
}


# conditional effects plots for submodel
if(F){
  fit_sub_cases <- readRDS("results/fit_submodel_cases_2022_07_01.rds")
  resp <- "cases"
  vars_submodel <- 
    c("meanAge", "HIV.AIDS", "Ascariasis", 
      "Malaria", "PerUrb", "HCexpend", "cardioDR",
      "spatLag")
  
    # check convergence (Rhat<1.1): very good
  rhat(fit_sub_cases) %>% max # 1.005
  
  ce_plots <- pbmclapply(vars_submodel, function(X) plot_ce(fit_sub_cases,X,resp = "cases"), mc.cores = 5)
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 2, ncol = 4) %>% 
    annotate_figure(top = text_grob("cases", size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  #ggsave("plots/ce_submodel_cases_2022_07_01.pdf", width = 9, height = 9, device = cairo_pdf)

  }


## compare sds priors for cases
if(F){
  # load fits with different sds priors
  m <- list()
  m[["s1"]] <- readRDS("results/fit_subset1_cases_sreg_2022_03_31_sds1.rds")
  m[["s2.6"]] <- readRDS("results/fit_subset1_cases_sreg_2022_03_31_sds2.6.rds")
  m[["s5"]] <- readRDS("results/fit_subset1_cases_sreg_2022_03_31_sds5.rds")
  m[["s10"]] <- readRDS("results/fit_subset1_cases_sreg_2022_03_31_sds10.rds")
  
  source("00_functions.R")
  theme_set(theme_classic())
  
  # no discernible differences in the estimated smooths for differing sds priors
  c("meanAge","Ascariasis","rain","PerUrb","log_tests_pp","spatLag") %>% 
    map(function(v){
      m %>% imap(~.x %>% plot_ce(v)+ labs(subtitle = .y)) %>% ggpubr::ggarrange(plotlist = .,nrow = 1)
    }) %>% ggpubr::ggarrange(plotlist = .,ncol = 1) %>% 
    ggsave("cases_prior_comp.pdf", height = 15, width = 6, plot = .)
    

}


