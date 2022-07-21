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
library(loo)
rm(list=ls())

nlist <- rstan::nlist
source("99_functions.R")

# load data set and identify predictors
cvd19 <- read_csv("data/cvd19_2022_07_04.csv")
vars_numeric <- cvd19 %>% select(where(is.numeric), -pop,-contains(c("cases","death"))) %>% names; vars_numeric
vars_factor <- cvd19 %>% select(!where(is.numeric), -country, -reg) %>% names; vars_factor

# set sampler settings
seed <- 768021 # sample(1e6,1)
control = list(adapt_delta = 0.9)
cores = 4
chains = 4
iter = 2000
thin = max(chains*iter/2000,1) # save a maximum of 1000 samples
init_r = 0
run_date <- "2022_07_01"

# fit global model and select best predictors via regularised estimates
if(F){
  resp  <- "cases" # cases/deaths
  vars_no_end <- vars_numeric[!vars_numeric %in% c("Ascariasis", "Malaria")]
  f_mu <- mu_spline(resp,vars_no_end, vars_factor, k = 6, basis = "ts")
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
      file = paste0("results/fit_all_no_end_",resp,"_",run_date),
      file_refit = "always",
      control = control)
  }

# manually inspect fits
if(F){
  theme_set(theme_classic())
  fit_all_cases <- readRDS("results/fit_all_cases_2022_07_01.rds")
  
  # check convergence (Rhat<1.1): very good
  rhat(fit_all_cases) %>% max # 1.007
  
  ce_plots <- pbmclapply(vars_numeric, function(X) plot_ce(fit_all_cases,X,resp = "cases"), mc.cores = 5)
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 4, ncol = 5) %>% 
    annotate_figure(top = text_grob("cases", size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  ggsave("ce_plots_cases_all_vars_2022_07_04.pdf", width = 9, height = 9, device = cairo_pdf)
  
}

# fit cases submodel
if(F){
  resp <- "cases"
  drop_endemic <- T
  vars_endemic <- c("Ascariasis", "Malaria")
  vars_submodel <- c("meanAge", "HIV.AIDS", "PerUrb", "HCexpend", 
                     "cardioDR", "spatLag", "log_tests", "log_pop")
  if(!drop_endemic) vars_submodel <- c(vars_submodel,vars_endemic)
  f_cases_mu <- mu_spline(resp,vars_submodel, vars_factor, k = 6, basis = "tp")
  f_shape <- shape ~ 1 + (1|a|reg)
  form_cases <- bf(f_cases_mu,f_shape) + negbinomial() 
  
  # try different values of the scale for the sds prior: 1,2.6(default),5,10 
  pr <- set_prior("student_t(3, 0, 2.6)", class = "sds")
   
  if(drop_endemic) no_end <- "_no_end" else no_end <- ""
  
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
                   file = paste0("results/fit_submodel_",resp,no_end,"_",run_date),
                   control = control)
}

# loo for cases models
if(F){
  fit_all_cases <- readRDS("results/fit_all_cases_2022_07_04.rds")
  fit_all_cases$model
  fit_all_cases_old <- readRDS("../analysis_revision_2021_11_19/results/fit_all_cases_sreg_2022_03_21.rds")
  fit_sub_cases <- readRDS("results/fit_submodel_cases_2022_07_01.rds")
  fit_sub_no_end_cases <- readRDS("results/fit_submodel_cases_no_end_2022_07_01.rds")
  
  bayes_R2(fit_all_cases_old)
  bayes_R2(fit_all_cases)
  bayes_R2(fit_sub_cases)
  bayes_R2(fit_sub_no_end_cases)
  loo_R2(fit_all_cases)
  
  
  future::plan("multisession", workers = 11)
  loo_all <- loo(fit_all_cases, future = T) # 11 problematic obs
  loo_sub <- loo(fit_sub_cases, future = T) # 8 problematic obs
  loo_sub_no_end <- loo(fit_sub_no_end_cases, future = T) #  problematic obs
  loo_compare(loo_all,loo_sub,loo_sub_no_end)
  
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
  fit_ <- readRDS("results/fit_submodel_cases_2022_07_01.rds")
  fit_ <- readRDS("results/fit_submodel_cases_no_end_2022_07_01.rds")
  
  resp <- "cases"
  vars_submodel <- c("meanAge", "HIV.AIDS", 
                     #"Ascariasis", "Malaria", 
                     "PerUrb", "HCexpend", 
                     "cardioDR", "spatLag",
                     "log_tests", "log_pop")
  
    # check convergence (Rhat<1.1): very good
  rhat(fit_) %>% max # 1.005
  
  ce_plots <- pbmclapply(vars_submodel, function(X) plot_ce(fit_,X,resp = "cases", 
                                                            exp = (X=="log_pop")), mc.cores = 5)
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 2, ncol = 4) %>% 
    annotate_figure(top = text_grob("cases", size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  #ggsave("plots/ce_submodel_cases_no_end_2022_07_01.pdf", width = 9, height = 9, device = cairo_pdf)

  }

# compare sds priors for cases
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

# posterior predictive distribution
if(F){
  fit_01 <- readRDS("results/fit_all_cases_2022_07_01.rds")
  fit_04 <- readRDS("results/fit_all_cases_2022_07_04.rds")

  fit_ <- readRDS("results/fit_all_cases_2022_07_01.rds")
  fit_$model
  pop_set <- 1e7
  pop_set <- cvd19$pop
  
  pp_raw <- posterior_predict(fit_,
                              newdata = cvd19 %>% mutate(pop = pop_set),
                              offset = T
  ) %>% as.data.frame()
  
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
              #mean1 = median(y_pred[is.finite(y_pred)], na.rm = F),
              mean = median(y_pred, na.rm = F)
    ) %>% 
    left_join(cvd19 %>% select(country,y_obs = cases, reg, pop), by = "country") %>% 
    mutate(q975 = pmin(q975,log(pop_set)),
           q75 = pmin(q75,log(pop_set)),
           y_obs = log(y_obs)) %>% 
    #mutate(q975 = pmin(q975,log(pop)), y_obs = log(y_obs)) %>% 
    arrange(reg,mean) %>% 
    ggplot(aes(x = factor(country, levels = country))) +
    geom_linerange(aes(ymin = q025, ymax = q975, col = reg), alpha = 0.5, size = 0.5) +
    geom_linerange(aes(ymin = q25, ymax = q75, col = reg), size = 1.5, alpha = 1) +
    #geom_point(aes(y = mean, col = reg), alpha = 1, shape = 18, size = 2) +
    geom_line(aes(y = mean, group = reg), col = "white", size =0.6) +
    #geom_line(aes(y = mean1, group = reg), col = "blue") +
    #geom_line(aes(y = log(pop), group = reg), col = "green") +
    geom_point(aes(y = y_obs), size =0.7) +
    theme_classic() +
    theme(legend.position = "top") +
    scale_color_brewer(type = "div", palette = "Dark2") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    ylim(c(4,NA))
}

# loo-predictive distributions
if(F){
  
  fit_ <- readRDS("results/fit_all_no_end_cases_2022_07_01.rds")
  #future::plan("multisession", workers = 20)
  #loo_ <- loo(fit_, save_psis = T, reloo = T, chains = 2, cores = 2)
  
  loo_all_no_end <- readRDS("results/loo_all_no_end_cases_2022_07_01.rds")
  loo_all_end <- readRDS("results/loo_all_cases_2022_07_01.rds")
  loo_sub_no_end <- readRDS("results/loo_submodel_no_end_cases_2022_07_01.rds")
  loo_sub_end <- readRDS("results/loo_submodel_cases_2022_07_01.rds")
  loo_compare(nlist(loo_sub_end,loo_sub_no_end))

  
  loo_compare(nlist(loo_all_no_end,loo_all_end))
  
  cvd19 %>% select(country,reg,Malaria,Ascariasis) %>% 
    mutate(loo_all_no_end = loo_all_no_end$pointwise[,"elpd_loo"],
           loo_all_end = loo_all_end$pointwise[,"elpd_loo"]) %>% 
    mutate(loo_diff = loo_all_no_end - loo_all_end) %>% 
    #filter(reg == "Afri") %>% 
    #filter(Ascariasis < 0.1 | Malaria < 0.1) %>% 
    summarise(ldiff = sum(loo_diff), se_diff = sd(loo_diff)*sqrt(n())) 
    
    
  
  #filter(Ascariasis > 0.01) 
    
    
    loo_1$pointwise[,"elpd_loo"]
  
  fit_ <- readRDS("results/fit_all_no_end_cases_2022_07_01.rds")
  fit_$formula
  ppred <- posterior_predict(fit_)
  loo_pred1 <- loo::E_loo(log(ppred),
                          psis_object = loo_all_end$psis_object, 
                          log_ratios = -1*log_lik(fit_), 
                          type = "quantile", 
                          probs = c(0.05,0.25,0.5,0.75,0.95))
  
  loo_pred_probs1 <- loo_pred1$value %>% t %>% as.data.frame() %>% as_tibble()
  colnames(loo_pred_probs1) <- c("ll","l","m","h","hh")
  
  
  bind_cols(cvd19 %>% select(country,reg,cases),loo_pred_probs1) %>% 
    mutate(y_obs = log(cases)) %>% 
    arrange(reg,m) %>% 
    ggplot(aes(x = factor(country, levels = country))) +
    geom_linerange(aes(ymin = ll, ymax = hh, col = reg), alpha = 0.5, size = 0.5) +
    geom_linerange(aes(ymin = l, ymax = h, col = reg), size = 1.5, alpha = 1) +
    geom_line(aes(y = m, group = reg), col = "white", size =0.6) +
    geom_point(aes(y = y_obs), size =0.7) +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_brewer(type = "div", palette = "Dark2") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    #ylim(c(4,NA)) +
    labs(title = "Leave-one-country-out posterior predictive distribution", 
         x = "Countries (ordered by region and predicted mean cases)",
         y = "log(cases)")
  
  #ggsave("plots/loo_pp_cases_2022_07_01.pdf", height = 6, width = 8)
  
  # coverage
  bind_cols(cvd19 %>% select(country,reg,cases),loo_pred_probs1) %>% 
    mutate(y_obs = log(cases)) %>% 
    mutate(cov90 = (hh >= y_obs) & (ll <= y_obs),
           cov50 = (h >= y_obs) & (l <= y_obs)) %>% 
    summarise(c90 = sum(cov90)/n(), c50 = sum(cov50)/n())
}


