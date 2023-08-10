#-------------------------------------
# 2022_09_29
#-------------------------------------

#library(tidyverse)
library(brms)
library(pbmcapply)
library(ggpubr)
library(loo)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
rm(list=ls())

#nlist <- rstan::nlist
theme_set(theme_classic())
source("99_functions.R")


# load data set and identify predictors ----
cvd19 <- read_csv("data/cvd19_2022_10_17.csv")
vars_endemic <- c("PathCompPC1","PathCompPC2","PathCompPC3","IDburd_perMil"
                  ); vars_endemic 
vars_numeric <- cvd19 %>% select(where(is.numeric), -log_tests,  -contains(c("cases","death","pop")), -all_of(vars_endemic)) %>% 
  names; vars_numeric
vars_factor <- cvd19 %>% select(!where(is.numeric), -country, -reg) %>% names; vars_factor

# set sampler settings ----
seed <- 750122 # sample(1e6,1)
control = list(adapt_delta = 0.98)
cores = 4
chains = 4
iter = 9000
thin = max(chains*iter/2000,1) # save a maximum of 1000 samples
init_r = 0
run_date <- "2022_10_18"

# set predictors per response and model type ----
resp <- "deaths" # cases, deaths
mod <- "all" # "sub" or "all"
offset <- 'log(cases)'

if(mod == "sub") f_mu <- mu_spline_subModel(response=resp,offset= offset, vars_numeric=vars_numeric, vars_factor=vars_factor)

if(mod == "all") f_mu <- mu_spline_all(response = resp, offset=offset, vars_numeric=vars_numeric, vars_endemic=vars_endemic, vars_factor=vars_factor)
  

  f_shape <- shape ~ 1 + (1|a|reg)
  f_model <- bf(f_mu,f_shape) + negbinomial(); f_model

# fit models ----
if(T){
  

  brm(f_model, 
      data = cvd19, 
      cores = cores, 
      chains = chains,
      init_r = init_r, 
      iter = iter, 
      thin = thin, 
      seed = seed, 
      file = paste0("results/fit_",mod,"_",resp,"_",run_date,".rds"),
      file_refit = "always",
      control = control)
  }

# inspect fits and plot conditional effects ----
if(F){
  fit_ <- readRDS(paste0("results/fit_",mod,"_",resp,"_",run_date,".rds"))
  rhat(fit_) %>% max # 1.013, check convergence (Rhat<1.1)

  ce_type <- "endemic" # "numeric", "endemic" # choose vars
  
  if(ce_type == "numeric") vars_plot <- vars_numeric else vars_plot <- vars_endemic
  ce_plots <- pbmclapply(vars_plot, function(X) plot_ce(fit_,X,resp = "deaths",exp = ifelse(stringr::str_starts(X,"log"),T,F)), 
                          mc.cores =1) #1 core supported on windows
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 3, ncol = 4) %>% 
    annotate_figure(top = text_grob("cases", size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  ggsave(paste0("plots/ce_",resp,"_",mod,"_",ce_type,"_",run_date,".pdf"), 
         plot = main_plot, width = 9, height = 9, device = cairo_pdf)
}

# LOO ----
if(F){
  
  run_date <- "2022_09_02"
  resp <- "deaths"
  
  fit_all <- readRDS(paste0("results/fit_all_",resp,"_",run_date,".rds"))
  fit_sub <- readRDS(paste0("results/fit_sub_",resp,"_",run_date,".rds"))


  r2_bayes(fit_all)
  r2_bayes(fit_sub)
  a <- model_performance(fit_sub)
  b <- model_performance(fit_all)
  
  future::plan("multisession", workers = 4)
  loo_all <- loo(fit_all, future = T) # 17 problematic obs
  loo_sub <- loo(fit_sub, future = T) # 13 problematic obs
  loo_compare(loo_all,loo_sub)
  
  future::plan("multisession", workers = 34)
  #fit_all_cases %>% loo(reloo = T) %>% saveRDS("results/loo_all_cases_2022_07_21.rds") # 34 refits
  #fit_submodel_cases %>% loo(reloo = T) %>% saveRDS("results/loo_submodel_cases_2022_07_21.rds") # 13 refits
  loo_all <- readRDS("results/loo_all_cases_2022_07_21.rds")
  loo_sub <- readRDS("results/loo_submodel_cases_2022_07_21.rds")
  nlist(loo_all, loo_sub) %>% loo_compare()
  
  
  # loo by region plot ----
  cvd19 %>% select(country, reg, Ascariasis, HIV.AIDS, Malaria) %>% 
    mutate(loo_diff = -loo_all$pointwise[,1] + loo_sub$pointwise[,1]) %>% 
    group_by(reg) %>% 
    #filter(Ascariasis > 0.05 | HIV.AIDS >0.025 | Malaria > 0.05) %>% 
    summarise(loo_diff_est = sum(loo_diff), loo_diff_sd = sd(loo_diff)*sqrt(n())) %>% 
    mutate(sig = abs(loo_diff_est) > loo_diff_sd & loo_diff_est < 0) %>% 
    ggplot(aes(x = reg)) +
    geom_point(aes(y = loo_diff_est, col = sig), size = 3) +
    geom_linerange(aes(ymin=loo_diff_est-loo_diff_sd, ymax=loo_diff_est+loo_diff_sd, col = sig)) +
    geom_hline(aes(yintercept = 0)) +
    labs(y = "LOO difference", title = "Cases: LOO scores by region")
  
  #ggsave("plots/loo_cases_by_region.pdf")

}

# posterior predictive distribution ----
if(F){
  resp <- "deaths"
  mod <- "all"
  run_date <- "2022_09_01"
  cvd19 <- read_csv("data/cvd19_2022_07_04.csv")
  
  fit_ <- readRDS(paste0("results/fit_",mod,"_",resp,"_",run_date,".rds"))
  pop_set <- cvd19$pop # cvd19$pop or 1e6 etc.
  pp_raw <- posterior_predict(fit_, newdata = cvd19 %>% 
                                mutate(log_tests = log((exp(log_tests)/pop)*pop_set), 
                                       log_pop = log(pop_set))) %>% 
    as.data.frame()
  
  colnames(pp_raw) <- cvd19$country
  pp_raw %>% 
    as_tibble() %>% 
    mutate(across(everything(),log)) %>% 
    pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
    group_by(country) %>% 
    #filter(is.finite(y_pred)) %>% 
    summarise(q025 = quantile(y_pred, probs = 0.025, na.rm = F),
              q975 = quantile(y_pred, probs = 0.975, na.rm = F),
              q25 = quantile(y_pred, probs = 0.25, na.rm = F),
              q75 = quantile(y_pred, probs = 0.75, na.rm = F),
              #mean1 = log(mean(exp(y_pred), na.rm = F)),
              mean1 = median(y_pred[is.finite(y_pred)], na.rm = F),
              mean = median(y_pred, na.rm = F)
    ) %>% 
    left_join(cvd19 %>% select(country,y_obs = cases, reg, pop), by = "country") %>% 
    mutate(y_obs = (y_obs/pop)*pop_set,
           q975 = pmin(q975,log(pop_set)),
           q75 = pmin(q75,log(pop_set)),
           y_obs = log(y_obs)) %>% 
    #mutate(q975 = pmin(q975,log(pop)), y_obs = log(y_obs)) %>% 
    arrange(reg,mean) %>% 
    ggplot(aes(x = factor(country, levels = country))) +
    geom_linerange(aes(ymin = q025, ymax = q975, col = reg), alpha = 0.5, size = 0.5) +
    geom_linerange(aes(ymin = q25, ymax = q75, col = reg), size = 1.5, alpha = 1) +
    #geom_point(aes(y = mean, col = reg), alpha = 1, shape = 18, size = 2) +
    geom_line(aes(y = mean, group = reg), col = "white", size =0.6) +
    geom_line(aes(y = mean1, group = reg), col = "blue") +
    #geom_line(aes(y = log(pop), group = reg), col = "green") +
    geom_point(aes(y = y_obs), size =0.7) +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_brewer(type = "div", palette = "Dark2") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    ylim(c(3,NA)) +
    labs(title = "Global cases", 
         subtitle = "Posterior predictive distribution",
         x = "Countries (ordered by region and predicted mean cases)",
         y = "log(cases) per 1,000,000 people") # edit as per pop_set
  
  #ggsave("plots/global_cases_ppd_perMillion_2022_09_01withop.pdf", height = 6, width = 8)
  
  
}



