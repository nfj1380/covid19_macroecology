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


cvd19 <- left_join(cvd19_pre, pop, by = "country") %>% 
  mutate(cases = round(cases_per_mil*pop*(1e6)^-1,0),
         deaths = round(deaths_per_mil*pop*(1e6)^-1,0),
         log_tests= log(tests),
         log_cases = log(cases),
         log_pop = log(pop),
         tests = NULL) %>% 
         mutate(across(where(is.character), as.factor))

cvd19 %>% names
#---------------------------------------------------------

#cvd19 %>% dplyr::select(where(is.numeric)) %>% cor %>% corrplot::corrplot()

# prepare variables and formula constructors--------------
vars_numeric <- cvd19 %>% 
  select(where(is.numeric), -contains("cases"),
         -contains("deaths"), -pop, -invSimp, -GDPpc, -IDBurd) %>% names; vars_numeric
vars_factor <- cvd19 %>% select(where(is.factor), -country) %>% names; vars_factor

# regional dependence: 
mu_spline <- function(response, vars_num, vars_fac, k = 6, basis = "tp"){
  paste(response," ~ 1 + ", paste(vars_fac, collapse = " + "), " + ", 
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
iter = 9000
thin = 8 # to leave 1000 samples
init_r = 0
run_date <- "2022_06_24"
#---------------------------------------------------------

# fit all vars to select best predictors: cases/deaths----
if(F){
  resp  <- "cases" # cases/deaths
  type <- "spline" # "spline" or "linear"
  shape <- "reg" # "1" or "reg"
  message(paste("Fitting model: ",resp, type, shape))
  if(resp=="deaths") vars_numeric <- c("log_cases",vars_numeric)
  f_mu <- mu_spline(resp,vars_numeric, vars_factor, k = 6, basis = "ts")
  if(type == "linear") f_mu <- mu_linear(resp, c(vars_factor, vars_numeric))
  f_shape <- shape ~ 1 
  if(shape == "reg") f_shape <- shape ~ 1 + 1|reg
  f_model <- bf(f_mu,f_shape) + negbinomial(); f_model
  
  #make_stancode(f_model, data = cvd19) %>% cat(file = "brms_code.stan")
  st_data <- make_standata(f_model, data = cvd19)
  st_data$X
  st_data$knots_1
  st_data$Zs_1_1 %>% head
  
  
  library(mgcv)
  
s1 <- mgcv::smooth.con(s(meanAge, k = 6, bs = "ts"),cvd19,NULL)
s1 <- smoothCon(s(meanAge, k = 6, bs = "ts"),cvd19, NULL, absorb.cons = T, diagonal.penalty = T, modCon = 3)
s1[[1]]$X %>% head
  
  s1$UZ
  
  s(meanAge, k = 6, bs = "ts")
  
  brm(f_model, 
      data = cvd19, 
      cores = cores, 
      chains = chains,
      init_r = init_r, 
      iter = iter, 
      thin = thin, 
      seed = seed, 
      file = paste0("results/fit_all_",resp,"_s",shape,"_",run_date),
      control = control)
  }

# fit cases submodel-------------------------------------
# subset1: meanAge, HIV.AIDS, Ascariasis, rain, PerUrb, log_tests_pp, spatLag
if(F){
  resp <- "cases"
  shape <- "reg" # "1" or "reg"
  f_cases_mu <- cases ~ 
    1 + 
    reg + 
    trans + 
    s(meanAge, k = 6, bs = "tp") + 
    s(HIV.AIDS, k = 6, bs = "tp") + 
    s(Ascariasis, k = 6, bs = "tp") + 
    s(rain, k = 6, bs = "tp") + 
    s(PerUrb, k = 6, bs = "tp") + 
    s(log_tests, k = 6, bs = "tp") + 
    s(spatLag, k = 6, bs = "tp")  +
    s(log_pop, k = 6, bs = "tp") 
  
  f_shape <- shape ~ 1 
  if(shape == "reg") f_shape <- shape ~ 1 + 1|reg
  
  form_cases <- bf(f_cases_mu,f_shape) + negbinomial() # nb: negbinomial()
  #form_cases <- bf(f_cases_mu) + poisson()  # pois: poisson()
  form_cases
  
  # try different values of the scale for the sds prior: 1,2.6(default),5,10 
  #run_date <- "2022_03_31_sds2.6"
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
                   file = paste0("results/fit_subset1_",resp,"_s",shape,"_",run_date),
                   control = control)
  
  fit_cases %>% conditional_smooths()

}

## loo for cases models
if(F){
  m.all_cases_sreg <- readRDS("results/fit_all_cases_sreg_2022_03_21.rds")
  m.all_cases_s1 <- readRDS("results/fit_all_cases_s1_2022_03_21.rds")
  m.sub_cases_sreg <- readRDS("results/fit_subset1_cases_tp_nb_sreg2022_03_21.rds")
  
  m.sub_cases_sreg$data
  pp_raw %>% as.data.frame()
  #future::plan("multisession", workers = 21)
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


#----- fit deaths submodel------------------------------------
if(F){
  # subset1: log_prop_cases, meanAge, Malaria, Trichuriasis, HCexpend, Hookworm, Schistosomiasis, log_tests_pp, spatLag 
  resp <- "deaths"
  shape <- "reg" # "1" or "reg"
  
  f_deaths_mu <- deaths|rate(pop) ~ 
    1 + 
    reg + 
    s(log_prop_cases, k = 6, bs = "tp") + 
    s(meanAge, k = 6, bs = "tp") + 
    s(Malaria, k = 6, bs = "tp") + 
    s(Trichuriasis, k = 6, bs = "tp") + 
    s(HCexpend, k = 6, bs = "tp") + 
    s(Hookworm, k = 6, bs = "tp") + 
    s(Schistosomiasis, k = 6, bs = "tp") +
    s(log_tests_pp, k = 6, bs = "tp") + 
    s(spatLag, k = 6, bs = "tp") 
  
    f_shape <- shape ~ 1 + 1|reg
  
    form_deaths <- bf(f_deaths_mu, f_shape) + negbinomial()
    
    # fit and save deaths models (< 6 mins)
    fit_deaths <- brm(form_deaths, 
                 data = cvd19, 
                 cores = cores, 
                 chains = chains,
                 init_r = init_r, 
                 iter = iter, 
                 thin = thin, 
                 seed = seed, 
                 file = paste0("results/fit_subset1_",resp,"_s",shape,"_",run_date),
                 control = control)
}

# Model comparison
# sub1 vs full: comparable for cases, sub1 better for deaths
if(F){
  if(F){
    m.all_deaths_sreg <- readRDS("results/fit_all_deaths_sreg_2022_03_21.rds")
    m.all_deaths_s1 <- readRDS("results/fit_all_deaths_s1_2022_03_21.rds")
    m.sub_deaths_sreg <- readRDS("results/fit_subset1_deaths_sreg_2022_03_21.rds")
    
    #future::plan("multisession", workers = 22)
    #m.all_deaths_sreg %>% loo(reloo = T, future = T) %>% saveRDS("results/loo_all_deaths_sreg.rds") # 43 refits
    #m.all_deaths_s1 %>% loo(reloo = T, future = T) %>% saveRDS("results/loo_all_deaths_s1.rds") # 163 refits
    #m.sub_deaths_sreg %>% loo(reloo = T, future = T) %>% saveRDS("results/loo_subset1_deaths_sreg.rds") # 30 refits
    
    list(all_sreg = readRDS("results/loo_all_deaths_sreg.rds"), 
         all_s1 = readRDS("results/loo_all_deaths_s1.rds"), 
         sub_sreg = readRDS("results/loo_subset1_deaths_sreg.rds")) %>%  
      loo_compare()
    
    #           elpd_diff se_diff
    # sub_sreg    0.0       0.0 
    # all_sreg  -34.6      55.4 
    # all_s1   -131.8      61.1
    
  }
}

# posterior predictive plots
if(F){
  
  cvd19 %>% arrange(cases) %>% relocate(cases) %>% mutate(cases = log(cases))
  
  m.sub_cases_sreg$model
  
  m.sub_cases_sreg <- readRDS("results/fit_subset1_cases_tp_nb_sreg2022_03_21.rds")
  identical(m.sub_cases_sreg$data$cases, cvd19$cases)
  
  pop_set <- 1e6
  pop_set <- cvd19$pop
  pp_raw <- posterior_predict(fit_cases)
                              
  pp_raw <- posterior_predict(fit_cases,
                             newdata = cvd19 %>% mutate(log_pop = log(pop_set))
                             ) %>% 
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
              q50 = quantile(y_pred, probs = 0.50, na.rm = F),
              #mean1 = log(mean(exp(y_pred), na.rm = F)),
              #mean1 = median(y_pred[is.finite(y_pred)], na.rm = F),
              mean = median(y_pred, na.rm = F)
              ) %>% 
    left_join(cvd19 %>% select(country,y_obs = cases, reg, pop), by = "country") %>% 
    mutate(q975 = pmin(q975,log(max(pop))),
           q75 = pmin(q75,log(max(pop))),
           y_obs = log(y_obs)) %>% 
    #mutate(q975 = pmin(q975,log(pop)), y_obs = log(y_obs)) %>% 
    arrange(reg,mean) %>% 
    ggplot(aes(x = factor(country, levels = country))) +
    geom_linerange(aes(ymin = q025, ymax = q975, col = reg), alpha = 0.7) +
    geom_linerange(aes(ymin = q25, ymax = q75, col = reg), size = 1.5, alpha = 0.7) +
    #geom_point(aes(y = mean, col = reg), alpha = 1, shape = 18, size = 2) +
    #geom_line(aes(y = mean, group = reg), col = "red") +
    geom_line(aes(y = q50, group = reg), col = "white") +
    #geom_line(aes(y = log(pop), group = reg), col = "green") +
    geom_point(aes(y = y_obs), size =0.8) +
    theme_classic() +
    theme(legend.position = "top", axis.text.x = element_blank() ) +
    scale_color_brewer(type = "div", palette = "Dark2") +
    guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL)) +
    #ylim(c(0,NA)) +
    labs(x = "country", y = "log(cases)", title = "Posterior predictions")

}


