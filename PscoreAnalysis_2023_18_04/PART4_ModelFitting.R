#######################################################################
########## PART 4: Model fitting
#######################################################################
pacman::p_load("tidyverse","brms", 
               'sp', 'spdep', 'spatialreg',
               'pbmcapply', 'ggpubr', 'loo', 'performance', 'projpred'
  )

rm(list=ls())

#nlist <- rstan::nlist
theme_set(theme_bw())
source("utility_functions.R")


# load data set and identify predictors ----
cvd19 <- read_csv("data/cvd19_2023_04_19.csv") 
W <- as.matrix(read.csv("data/neighbor.matrix.csv", row.names=1))

#spatial autocorrelation
W[lower.tri(W)]<-t(W)[lower.tri(W)] #make dure it is symmetric

vars_endemic <- c("Malaria","HIV.AIDS", 'PathCompPC1', "PathCompPC3",'PathCompPC2', "IDburd_perMil")

#cvd19 %>% select(where(is.numeric)) %>% cor()
#PCA1 - pos values correlate to L.filariasis, hookworm, schisto,Malaria
#PCA2 _negative values herpes and HIV
#PCA3 - negative values with high schisto/Malaria. Pos with tb.
vars_numeric <- cvd19 %>% select(where(is.numeric),  -contains(c("cases","death",'tests','excess','mortality',
                                                                 'pop','pscore', 'lat', 'long', 'temp', 'rainfall', 'Vaccinated')), -all_of(vars_endemic)) %>% 
  names; vars_numeric
#rainfall and temp were not found to be important in the last paper.
#removed vaccinated as 
vars_factor <- cvd19 %>% select(!where(is.numeric), -country, -reg) %>% names; vars_factor

# vars_numeric_enviro <- cvd19 %>% select( temp, rainfall ) %>%
#                                          
#   names; vars_numeric_enviro  
# 
# vars_numeric_demographic <- cvd19 %>% select(where(is.numeric),  -contains(c("cases","death", 'tests','excess','mortality', 'pop')), -all_of(vars_endemic),
#                                        -all_of(vars_numeric_enviro)) %>% 
#   names; vars_numeric_demographic

# set sampler settings ----
seed <- 750122 # sample(1e6,1)
control = list(adapt_delta = 0.999)
cores = 4
chains = 4
iter = 12000
thin = max(chains*iter/2000,1) # save a maximum of 1000 samples
init_r = 0
run_date <- "2023_27_05" #with vaccinated 2021 estimate

# set predictors per response and model type ----
resp <- "pscore.mean" # WHO excess mortality estimates
mod <- 'sub' # 'lm',"sub" or "all"

str(cvd19)

if(mod == "sub") f_mu <- mu_spline_subModel(response=resp, vars_numeric=vars_numeric, vars_factor=vars_factor)

if(mod == "all") f_mu <- mu_spline_all(response = resp,  vars_numeric=vars_numeric, vars_endemic=vars_endemic, k_end=5, k_num=5)

if(mod == "enviro") f_mu <- mu_spline_all(response=resp, vars_numeric=vars_numeric_enviro,  vars_endemic=vars_endemic)

if(mod == "lm") f_mu <- mu_linear(response=resp, vars_numeric=vars_numeric,  vars_endemic=vars_endemic)

if(mod == "interact") f_mu <- pscore.mean | mi(pscore.mean.se) ~ 1 + (1 | reg) + trans + s(density, k = 5, bs = "ts") + t2(GDPpc,PerUrb) + s(Diabetes_prevalence, k = 5, bs = "ts") + s(cardioDR, k = 5, bs = "ts") + s(Hospital_beds_per_thousand, k = 5, bs = "ts") + s(PropOver_70, k = 5, bs = "ts") + s(AirConnectivityIndex, k = 5, bs = "ts")+ s(lag_rate, k = 5, bs = "ts") + s(Malaria, k = 5, bs = "ts") + s(HIV.AIDS, k = 5, bs = "ts") + s(PathCompPC1, k = 5, bs = "ts") + s(PathCompPC3, k = 5, bs = "ts") + s(PathCompPC2, k = 5, bs = "ts") + s(IDburd_perMil, k = 5, bs = "ts") 


#f_model <- bf(f_mu) + gaussian(); f_model

   f_shape <- sigma ~ 1 + (1|a|reg)
   f_model <- bf(f_mu,f_shape) + gaussian(); f_model
   #f_model <- bf(f_mu,f_shape) + student(); f_model #doesnt change anything.Outliers no better predicted

# fit models ----
if(T){
  
  #got some 0s in log cases
  #cvd19$log_cases[cvd19$log_cases==0]=0.1

  brm(f_model, 
      data = cvd19, 
      #data2 = list(W = W),#spatial
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

  
  plot(conditional_effects(fit_ ), ask = FALSE)
  
  
cond <- conditional_effects(fit_, surface = TRUE)
plot(cond, rug = TRUE, stype='raster')[[4]]

ce_type <- "endemic" # "numeric", "endemic" # choose vars

  #if linear model
  #stanplot(fit_ , pars = "^b_")
  #reg_level = 'Afri'
  reg_level = 'Afri'
  Transmission_Classification = NA
  
  if(ce_type == "numeric") vars_plot <- vars_numeric else vars_plot <- vars_endemic
  ce_plots <- pbmclapply(vars_plot, function(X) plot_ce(fit_,X,resp = "pscore.mean",reg_level = reg_level),#,
                                                       # exp = ifelse(stringr::str_starts(X,"log"),T,F)), 
                          mc.cores =1) #1 core supported on windows
  main_plot <- ggarrange(plotlist = ce_plots, nrow = 3, ncol = 4) %>% 
    annotate_figure(top = text_grob("pscore", size = 16, face = "bold"), bottom = paste0("CI: 90%"))
  ggsave(paste0("plots/ce_",resp,"_",mod,"_",reg_level,"_", ce_type,"_",run_date,".pdf"), 
         plot = main_plot, width = 9, height = 9, device = cairo_pdf)
  
  
  #for each region
  conditions <- data.frame(reg = unique(cvd19$reg))
  rownames(conditions) <- unique(cvd19$reg)
  
  me_loss <- conditional_effects(
    fit_, conditions = conditions,
    re_formula = NULL, method = "predict")
  
  plot(me_loss, ncol = 5, points = TRUE)
  
  
  #experimental decomposition
  
  devtools::install_github("jmgirard/varde")
  
  library(varde)
  
  res_1 <- varde(fit_)
  
  plot(res_1, type = "river")
  
  pp_samples <- posterior_predict(fit_, draws = 1000)
  glimpse(pp_samples)
  # Extract observed y values

}

 
  
  
# LOO ----
if(F){
  
  run_date <- "2023_19_05"
  mod <- 'sub' 

  # fit_all <- readRDS(paste0("results/fit_all_",resp,"_",run_date,".rds"))
  # #weird performance for fit_all
  # fit_sub <- readRDS(paste0("results/fit_sub_",resp,"_",run_date,".rds"))
  
  fit_all <- readRDS('results/fit_all_pscore.mean_2023_27_05.rds')
 
  fit_sub <- readRDS("results/fit_sub_pscore.mean_2023_27_05.rds")
  
  #vs <- cv_varsel( fit_all, cv_method='LOO')#doesnt work

  # r2_bayes(fit_sub)
  # r2_bayes(fit_all)
  # a <- model_performance(fit_sub)
  # b <- model_performance(fit_all)
  
 # future::plan("multisession", workers = 4)
  loo_all <- loo(fit_all, future = T) # 14 problematic obs
  loo_sub <- loo(fit_sub, future = T) # 15 problematic obs
  loo_compare(loo_all,loo_sub)
  
  #adding interaction term didnt improve model fit
  future::plan("multisession", workers = 4)
  loo_all_fixed <- fit_all %>% loo(reloo = T) %>% saveRDS("results/loo_all_2023_27_05.rds") #%>% saveRDS("results/loo_all_deaths_2022_11_08.rds") # 34 refits
  loo_all_sub <-  fit_sub %>% loo(reloo = T) %>% saveRDS("results/loo_submodel_2023_27_05.rds") # 13 refits
  
  loo_all <- readRDS("results/loo_all_2023_27_05.rds")
  loo_sub <- readRDS("results/loo_submodel_2023_27_05.rds")
 # nlist(loo_all, loo_sub) %>% loo_compare()
  
  
  # # loo by region plot ----
  # cvd19 %>% select(country, reg) %>% 
  #   mutate(loo_diff = -loo_all$pointwise[,1] + loo_sub$pointwise[,1]) %>% 
  #   group_by(reg) %>% 
  #   #filter(Ascariasis > 0.05 | HIV.AIDS >0.025 | Malaria > 0.05) %>% 
  #   summarise(loo_diff_est = sum(loo_diff), loo_diff_sd = sd(loo_diff)*sqrt(n())) %>% 
  #   mutate(sig = abs(loo_diff_est) > loo_diff_sd & loo_diff_est < 0) %>% 
  #   ggplot(aes(x = reg)) +
  #   geom_point(aes(y = loo_diff_est, col = sig), size = 3) +
  #   geom_linerange(aes(ymin=loo_diff_est-loo_diff_sd, ymax=loo_diff_est+loo_diff_sd, col = sig)) +
  #   geom_hline(aes(yintercept = 0)) +
  #   labs(y = "LOO difference", title = " LOO scores by region")
  
  #ggsave("plots/loo_cases_by_region.pdf")

}

# posterior predictive distribution ----
if(F){
  resp <- "deaths"
  mod <- "all"
  run_date <- "2022_11_10"
  cvd19 <- read_csv("data/cvd19_2022_07_04.csv")
  
  fit_ <- readRDS(paste0("results/fit_",mod,"_",resp,"_",run_date,".rds"))
  pp_check(fit_) #check overall - not great 
  
  pop_set <- cvd19$pop # cvd19$pop or 1e6 etc.
  pp_raw <- posterior_predict(fit_, newdata = cvd19 %>% 
                                mutate(log_tests = log((exp(log_tests)/pop)*pop_set), 
                                       log_pop = log(pop_set))) %>% 
    as.data.frame()
  
  colnames(pp_raw) <- cvd19$country
  pp_raw %>% 
    as_tibble() %>% 
    #mutate(across(everything(),log)) %>% 
    pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
    group_by(country) %>% 
    filter(is.finite(y_pred)) %>% 
    summarise(q025 = quantile(y_pred, probs = 0.025, na.rm = F),
              q975 = quantile(y_pred, probs = 0.975, na.rm = F),
              q25 = quantile(y_pred, probs = 0.25, na.rm = F),
              q75 = quantile(y_pred, probs = 0.75, na.rm = F),
              #mean1 = log(mean(exp(y_pred), na.rm = F)),
              mean1 = median(y_pred[is.finite(y_pred)], na.rm = F),
              mean = median(y_pred, na.rm = F)
    ) %>% 
    left_join(cvd19 %>% select(country,y_obs =pscore.mean, reg, pop), by = "country") %>% 
    # mutate(y_obs = (y_obs/pop)*pop_set,
    #        q975 = pmin(q975,log(pop_set)),
    #        q75 = pmin(q75,log(pop_set)),
    #        y_obs = log(y_obs)) %>% 
    #mutate(q975 = pmin(q975,log(pop)), y_obs = log(y_obs)) %>% 
    arrange(reg,mean) %>% 
    ggplot(aes(x = factor(country, levels = country))) +
    geom_linerange(aes(ymin = q025, ymax = q975, col = reg), alpha = 1, size = 0.5) +
    geom_linerange(aes(ymin = q25, ymax = q75, col = reg), size = 1.5, alpha = 1) +
    #geom_point(aes(y = mean, col = reg), alpha = 1, shape = 18, size = 2) +
    geom_line(aes(y = mean, group = reg), col = "white", size =0.6) +
    #geom_line(aes(y = mean1, group = reg), col = "blue") +
   # geom_line(aes(y = tests, group = reg), col = "red") +
    geom_point(aes(y = y_obs), size =0.7) +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 90, size = 6, hjust = 1, vjust = 0.5),
          panel.grid.major.x = element_line()) + scale_x_discrete(label=abbreviate)+
    # theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    scale_color_brewer(type = "div", palette = "Dark2") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    ylim(c(-50,NA)) +
    labs(title = "Global pscore", 
         subtitle = "Posterior predictive distribution",
         x = "Countries (ordered by region and pscores)",
         y = "pscore") # edit as per pop_set
  
  #ggsave("plots/global_deaths_ppd_perMillion_2022_09_30.pdf", height = 6, width = 8)
  
  
}



