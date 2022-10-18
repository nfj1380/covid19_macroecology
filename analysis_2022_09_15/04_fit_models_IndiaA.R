pacman::p_load("tidyverse","brms", 
               'sp', 'spdep', 'spatialreg',
               'pbmcapply', 'ggpubr')

rm(list=ls())

source("99_functions.R")

theme_set(theme_classic())

# load and prepare data set-------------------------------
cvd19_india <- read_csv("data/cvd19_India.csv")
summary(cvd19_india)

#density and urban

#Trich prev not corr with ascariasis in India
  vars_endemic <- c("PathCompPC1", "PathCompPC2",
                     "IDburd_perMil"); vars_endemic
  vars_numeric <- cvd19_india %>% 
    select(where(is.numeric),  -pop, -Lat , -Long,
           -Health_expenditure, -Density,
           -invSimp, -IDBurd, #PC1 stronger correlation with Malaria
           -contains(c("tested", "cases","death",'tests','pop')), -all_of(vars_endemic)) %>% 
  names; vars_numeric
  

vars_factor <- cvd19_india  %>% select(!where(is.numeric), -Province_State) %>% names; vars_factor

# set sampler settings ----
seed <- 750122 # sample(1e6,1)
control = list(adapt_delta = 0.97)
cores = 4
chains = 4
iter = 9000
thin = max(chains*iter/2000,1) # save a maximum of 1000 samples
init_r = 0
run_date <- "2022_10_18"

# set predictors per response and model type ----
resp <- "deaths" # cases, deaths
mod <- "all" # "sub" or "all"
offset <- 'log_cases' #or log_testspp

if (mod=='all') f_mu <- mu_splineIndia_all(resp, offset, vars_numeric, vars_endemic)
if (mod=='sub') f_mu <- mu_splineIndia_subModel(resp, offset, vars_numeric)
f_shape <- shape ~ 1
f_model <- bf(f_mu,f_shape) + negbinomial(); f_model

# fit models ----
if(T){

  
  #individual model
  brm(f_model, 
      data = cvd19_india, #first imputed data set
      cores = cores, 
      chains = chains,
      init_r = init_r, 
      iter = iter, 
      thin = thin, 
      seed = seed, 
      file = paste0("results/fitIndia_",mod,"_",resp,"_",run_date,".rds"),
      file_refit = "always",
      control = control)
}


# fit each imputed data ---- 

if(F){
future::plan("multisession", workers = 4)

# fit and save case models
fit_<- brm_multiple(f_model, 
                          data = cvd19_imp, # list of imputed data sets
                          combine = FALSE, # combine posteriors
                          future = TRUE, 
                          chains = chains,
                          init_r = init_r, 
                          iter = iter, 
                          thin = thin, 
                          seed = seed, 
                          #file = paste0("results/fitIndia_",mod,"_",resp,"_",run_date,".rds"),
                          file_refit = "always",
                          control = control)

saveRDS(fit_cases, 'fit_cases_India')
cases_comb <- combine_models(mlist = fit_cases, check_data = F)

}

# inspect fits and plot conditional effects ----
if(F){
  
  fit_<- readRDS(paste0("results/fitIndia_",mod,"_",resp,"_",run_date,".rds"))
  rhat(fit_) %>% max # 1.013, check convergence (Rhat<1.1)
  
  ce_type <- "numeric" # "numeric", "endemic" # choose vars
  
  if(ce_type == "numeric") vars_plot <- vars_numeric else vars_plot <- vars_endemic
  ce_plots <- pbmclapply(vars_plot, function(X) plot_ce(fit_,X,resp = 'deaths',exp = ifelse(stringr::str_starts(X,"log"),T,F)), 
                         mc.cores =1) #1 core supported on windows
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 3, ncol = 4) %>% 
    annotate_figure(top = text_grob(resp, size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  ggsave(paste0("plots/ceIndia_",resp,"_",mod,"_",ce_type,"_",run_date,".pdf"), 
         plot = main_plot, width = 9, height = 9, device = cairo_pdf)
}  


# posterior predictive distribution ----
if(F){
  resp <- "deaths"
  mod <- "all"
  fit_ <- readRDS(paste0("results/fitIndia_",mod,"_",resp,"_",run_date,".rds"))
  pop_set <- cvd19_india$pop # cvd19$pop or 1e6 etc.
  pp_raw <- posterior_epred(fit_, newdata = cvd19_india %>% 
                              mutate(log_tests_pp = log((exp(log_tests)/pop)*pop_set), 
                                     log_pop = log(pop_set))) %>% 
    as.data.frame()
  
  colnames(pp_raw) <- cvd19_india$Province_State
  pp_raw %>% 
    as_tibble() %>% 
    mutate(across(everything(),log)) %>% 
    pivot_longer(everything(), names_to = "state", values_to = "y_pred") %>% 
    group_by(state) %>% 
    #filter(is.finite(y_pred)) %>% 
    summarise(q025 = quantile(y_pred, probs = 0.025, na.rm = F),
              q975 = quantile(y_pred, probs = 0.975, na.rm = F),
              q25 = quantile(y_pred, probs = 0.25, na.rm = F),
              q75 = quantile(y_pred, probs = 0.75, na.rm = F),
              #mean1 = log(mean(exp(y_pred), na.rm = F)),
              mean1 = median(y_pred[is.finite(y_pred)], na.rm = F),
              mean = median(y_pred, na.rm = F)
    ) %>% 
    left_join(cvd19 %>% select(Province_State,y_obs = resp, reg, pop), by = "country") %>% 
    mutate(y_obs = (y_obs/pop)*pop_set,
           q975 = pmin(q975,log(pop_set)),
           q75 = pmin(q75,log(pop_set)),
           y_obs = log(y_obs)) %>% 
    #mutate(q975 = pmin(q975,log(pop)), y_obs = log(y_obs)) %>% 
    arrange(reg,mean) %>% 
    ggplot(aes(x = factor(Province_State, levels = Province_State))) +
    geom_linerange(aes(ymin = q025, ymax = q975, col = reg), alpha = 0.5, size = 0.5) +
    geom_linerange(aes(ymin = q25, ymax = q75, col = reg), size = 1.5, alpha = 1) +
    #geom_point(aes(y = mean, col = reg), alpha = 1, shape = 18, size = 2) +
    geom_line(aes(y = mean, group = reg), col = "white", size =0.6) +
    # geom_line(aes(y = mean1, group = reg), col = "blue") +
    #geom_line(aes(y = log(pop), group = reg), col = "green") +
    geom_point(aes(y = y_obs), size =0.7) +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_brewer(type = "div", palette = "Dark2") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    ylim(c(3,NA)) +
    labs(title = "Global deaths", 
         subtitle = "Posterior predictive distribution",
         x = "Countries (ordered by region and predicted mean cases)",
         y = "log(deaths) per 1,000,000 people") # edit as per pop_set
  
  ggsave("plots/global_deaths_ppd_perMillion_2022_08_21_noline.pdf", height = 6, width = 8)
  
}  
  

# #---------------------------------------------------------
# # cases submodel
# #---------------------------------------------------------
# 
# f_cases_mu <- cases|rate(pop) ~ 1+
#   #s(HIV, k = 6, bs = "tp") + prevalence too low (< 1%)
#     s(Ascariasis, k = 5, bs = "tp") + 
#   s(urban, k = 5, bs = "tp") + 
#   s(log_tests_pp, k = 5, bs = "tp") + 
#   s(lag_rate, k = 5, bs = "tp") + 
#   s(Mean_age, k = 5, bs = "tp")
# 
# 
# f_shape <- shape ~ 1
# 
# form_cases <- bf(f_cases_mu, f_shape) + negbinomial()
# form_cases
# future::plan("multisession", workers = 20)
# 
# # fit and save case models
# fit_cases <- brm_multiple(form_cases, 
#                  data = cvd19_imp, # list of imputed data sets
#                  combine = FALSE, # combine posteriors
#                  future = TRUE, 
#                  chains = chains,
#                  init_r = init_r, 
#                  iter = iter, 
#                  thin = thin, 
#                  seed = seed, 
#                  #file = paste0("results/fit_cases_India_imp_",run_date),
#                  control = control)
# 
# #saveRDS(fit_cases,"results/old_fits/fit_cases_mi_separate_2022_03_21.rds")
# 
# cases_comb <- combine_models(mlist = fit_cases, check_data = F)
# resp <- "cases"
# 
# cases_comb %>% conditional_effects("Ascariasis")
# 
# vars_cases <- c("Ascariasis","urban","log_tests_pp","lag_rate","Mean_age")
# 
# resp <- "cases"
# for(x in c("Ascariasis","urban","log_tests_pp","lag_rate","Mean_age")){
#   fit_cases %>% 
#     map(~.x %>% plot_ce(x)) %>% 
#     ggpubr::ggarrange(plotlist = .) %>% 
#     ggsave(paste0("plots/ce_mi_",x,".pdf"), plot = .)
# }
# 
# c("Ascariasis","urban","log_tests_pp","lag_rate","Mean_age") %>% 
#   map(~ plot_ce(cases_comb, .x)) %>% ggpubr::ggarrange(plotlist = .) %>% 
#   ggsave("plots/ce_mi_combined.pdf",plot = .)
# 
# 
# #---------------------------------------------------------
# # deaths submodel
# #---------------------------------------------------------
# 
# run_date <- "2022_03_30"
# 
# f_deaths_mu <- deaths|rate(pop) ~ 1+ 
#   s(log_cases , k = 4, bs = "tp") + 
#   s(Malaria, k = 4, bs = "tp") + 
#   s(Trichuriasis , k = 4, bs = "tp") + 
#   # s(Hookworm, k = 6, bs = "tp") + # prevalence too low
#   #s(log_tests_pp, k = 6, bs = "tp") + #tests not impacting deaths here 
#   s(lag_rate, k = 4, bs = "tp")  +  #doesn't look to be important
#   s(log_healthcare , k = 4, bs = "tp")+ 
#   s(Mean_age , k = 4, bs = "tp")
# 
# 
# f_shape <- shape ~ 1
# 
# form_deaths <- bf(f_deaths_mu, f_shape) + negbinomial()
# 
# future::plan("multisession", workers = 20)
# fit_deaths <- brm_multiple(form_deaths, 
#                           data = cvd19_imp, # list of imputed data sets
#                           combine = FALSE, # combine posteriors
#                           future = TRUE, 
#                           chains = chains,
#                           init_r = init_r, 
#                           iter = iter, 
#                           thin = thin, 
#                           seed = seed, 
#                           control = control)
# 
# #saveRDS(fit_deaths,paste0("results/old_fits/fit_deaths_mi_separate_",run_date,".rds"))
# 
# # deaths plots
# run_date <- "2022_03_30"
# #fit_deaths <- readRDS(paste0("results/old_fits/fit_deaths_mi_separate_",run_date,".rds"))
# deaths_comb <- combine_models(mlist = fit_deaths, check_data = F)
# 
# vars_deaths <- c("log_cases","Malaria", "Trichuriasis","lag_rate","log_healthcare","Mean_age")
# resp <- "deaths"
# 
# # plots for all predictors for all ten imputed data sets
# for(x in vars_deaths){
#   fit_deaths %>% 
#     map(~.x %>% plot_ce(x, exp = if_else(x == "log_cases",TRUE, FALSE))) %>% 
#     ggpubr::ggarrange(plotlist = .) %>% 
#     ggsave(paste0("plots/deaths_ce_mi_",x,"_",run_date,".pdf"), plot = .)
# }
# 
# # plots for all predictors for the combined model
# vars_deaths %>% 
#   map(~ plot_ce(deaths_comb, .x, exp = if_else(.x == "log_cases",TRUE, FALSE))) %>% 
#   ggpubr::ggarrange(plotlist = .) %>% 
#   ggsave(paste0("plots/deaths_ce_mi_combined_",run_date,".pdf"),plot = .)
# 
# # plot for Malaria for combined model
# run_date <- "2022_03_30"
# resp <- "deaths"
# deaths_comb <- readRDS(paste0("results/fits_deaths_mi_combined_",run_date))
# deaths_comb %>% plot_ce("Malaria", exp = F)
# ggsave(paste0("plots/deaths_India_Malaria_",run_date,".pdf"))
# 
