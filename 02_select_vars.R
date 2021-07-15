pacman::p_load("tidyverse",'brms')

rm(list=ls())


#make consistent dataframes
prepare_data_cases <- function(data){
  data %>% 
    select(-Deaths._perMil, - latitude, -longitude) %>% 
    rename(y = Cases._perMil, reg = WHO_Region, tests = Total_tests_per_thousand,
           rain = Avg_Rainfall91.2016_mm, temp = Avg_temp91.2016_c, 
           popDen = popDen_2018, trans = Transmission_Classification,
           meanAge = mean._age, country = X1,
           Diabetes = Diabetes_prevalence, Hookworm = Hookworm.disease,
           invSimp = inv_simpsons, cardioDR = Cardiovasc_death_rate,
           spatLag = lag_rate, IDBurd = inf_diseaseBurden) %>% 
    mutate(reg = reg %>% 
             str_replace(fixed("Eastern Mediterranean"), "EMed") %>% 
             str_replace(fixed("Western Pacific"), "WPac") %>% 
             str_replace(fixed("South-East Asia"), "SEAs") %>% 
             str_sub(1,4)) %>% 
    mutate(trans = trans %>% 
             str_replace(fixed("Clusters_of_cases"), "clust") %>% 
             str_replace(fixed("Community_transmission"), "comm") %>% 
             str_replace(fixed("Sporadic_cases"), "spor") %>% 
             str_replace(fixed("Pending"), "comm")) %>% 
    mutate(y = round(y)) %>% 
    mutate(across(where(is_character),factor))
}


prepare_data_deaths <- function(data){
  data %>% 
    select( -latitude, -longitude, -Total_tests_per_thousand, -Transmission_Classification) %>% 
    rename(y = Deaths._perMil, reg = WHO_Region,
           rain = Avg_Rainfall91.2016_mm, temp = Avg_temp91.2016_c, 
           popDen = popDen_2018, 
           meanAge = mean._age, country = X1,
           Diabetes = Diabetes_prevalence, Hookworm = Hookworm.disease,
           invSimp = inv_simpsons, cardioDR = Cardiovasc_death_rate,
           cases = Cases._perMil,
           spatLag = lag_rate, IDBurd = inf_diseaseBurden) %>% 
    mutate(reg = reg %>% 
             str_replace(fixed("Eastern Mediterranean"), "EMed") %>% 
             str_replace(fixed("Western Pacific"), "WPac") %>% 
             str_replace(fixed("South-East Asia"), "SEAs") %>% 
             str_sub(1,4)) %>% 
    mutate(y = round(y)) %>% 
    mutate(across(where(is_character),factor))
}

#import data

cvd19.raw <- read_csv("data/cvdDataPreprocessed")

############################################################################################
#####################CASE MODEL#####################
############################################################################################

cvd19cases <- prepare_data_cases(cvd19.raw); cvd19


# define covariates
vars_all <- setdiff(names(cvd19cases),c("country","y")); vars_all
vars_numeric <- cvd19cases %>% select(where(is.numeric),-y) %>% names; vars_numeric
vars_factor <- cvd19cases %>% select(where(is.factor), -country) %>% names; vars_factor

#remove all variables with 0 slopes based on the complete model

vars_reduced <- vars_numeric[-c(1,3,6:11, 13, 14:18,20:21)] 
#vars_reduced <- vars_numeric[-c(1,3,6:11, 13, 15:18,20:21)] #

#herpes, trich, tb, density (?), rainfall (?), temperature, GDP, HCspend, invSimp (?/ hinge), hookworm, schisto, Lmphatic,
#IDbud, diabetes, cardion (?, )


# functions to specify models
mu_linear <- function(vars) paste("y ~",paste(vars, collapse = " + ")) %>% as.formula()

mu_spline <- function(vars_num, vars_fac, k = 6, bas = "tp"){  # N.B. bas = "ts" adds shrinkage
  paste("y ~", paste(vars_fac, collapse = " + "), " + ", 
        paste0("s(",vars_num,", k = ", k,", bs = '", bas, "')", collapse = " + ")) %>% 
    as.formula()
} 

#regional splines
mu_hierachical_spline <- function(vars_num, vars_fac, k = 6, bas = "tp"){ 
  paste(paste("y ~", paste(vars_fac, collapse = " + "), " + ",
              paste0("s(",vars_num,", by = reg, k=",k,", bs='", bas, "')", collapse="+"))) %>% 
    as.formula()
}

# brms fitting options
future = F
cores = 4
chains = 4
iter = 6000
control = list(adapt_delta = 0.9)
init_r = 1

# specify model formulae
f.lin.s1.nbin<- bf(mu_linear(vars_all), shape ~ 1) + negbinomial()
f.lin.s2.nbin<- bf(mu_linear(vars_all), shape ~ 1 + reg) + negbinomial()
f.spl.s1.nbin.tp <- bf(mu_spline(vars_numeric, vars_factor), shape ~ 1) + negbinomial()
f.spl.s1.nbin.ts <- bf(mu_spline(vars_numeric, vars_factor, bas = "ts"), shape ~ 1) + negbinomial()
f.spl.s2.nbin.tp <- bf(mu_spline(vars_numeric, vars_factor), shape ~ 1 + reg) + negbinomial()
f.spl.s2.nbin.ts <- bf(mu_spline(vars_numeric, vars_factor, bas = "ts"), shape ~ 1 + reg) + negbinomial()
f.spl.s3.nbin.ts <- bf(mu_spline(vars_numeric, vars_factor, bas = "ts"), shape ~ 1 + reg+tests) + negbinomial()

#reduced sets. just for ts splines to start with. Including test in the error term led to model convergence issues.

f.spl.Reduced <- bf(mu_spline(vars_reduced, vars_factor, bas = "ts"), shape ~ 1 + reg ) + negbinomial()
f.spl.Reduced_tests <- bf(mu_spline(vars_reduced, vars_factor, bas = "ts"), shape ~ 1 + tests ) + negbinomial()

f.spl.Reduced_hierachical <- bf(mu_hierachical_spline (vars_reduced, vars_factor, bas = "ts"), shape ~ 1 ) + negbinomial()

#f.spl.Reduced <- bf(mu_spline(vars_reduced, vars_factor, bas = "ts"), shape ~ 1  + reg ) + negbinomial()
#test model doesn't work

# fit models in brms
#if(F){
#message(paste("Fitting model:", "f.spl.s2.nbin.ts"))


m1_full <-  brm(f.spl.s1.nbin.ts, data = cvd19cases,
                chains = chains, iter = iter, future = future, control = control, cores = cores, 
                init_r = init_r) 

m2_full_WHO <-  brm(f.spl.s2.nbin.ts, data = cvd19cases,
                    chains = chains, iter = iter, future = future, control = control, cores = cores, 
                    init_r = init_r) 

m3_reducedSet_WHO <-  brm(f.spl.Reduced , data = cvd19cases,
                          chains = chains, iter = iter, future = future, control = control, cores = cores, 
                          init_r = init_r) 

#m4_reducedSet_WHO_Tests <-  brm(f.spl.Reduced , data = cvd19,
# chains = chains, iter = iter, future = future, control = control, cores = cores, 
#init_r = init_r) #doesn't converge

m4_reducedSet_WHO <-  brm(f.spl.Reduced , data = cvd19cases,
                          chains = chains, iter = iter, future = future, control = control, cores = cores, 
                          init_r = init_r) 
m4_reducedSet_WHO %>% saveRDS("m4_reducedSetInSimp_WHO .rds")

m5_reducedSet_WHO_highK <-  brm(f.spl.Reduced , data = cvd19cases,
                                chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                init_r = init_r) 

m7_reducedSet_hierachical <-  brm(f.spl.Reduced_hierachical , data =cvd19cases,
                                  chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                  init_r = init_r)

summary(m7_reducedSet_hierachical )

m1loo <-  loo(m1_full  )

m2loo <-  loo(m2_full_WHO )

m3loo <-  loo(m4_reducedSet_WHO  )

m4loo <-  loo(m4_reducedSetInSimp_WHO ) 

m5loo <-  loo(m5_reducedSetInSimp_WHO_highK  )

m7loo <- loo(m7_reducedSet_hierachical)

loo_compare( m1loo, m2loo, m3loo )

loo_moment_match(m7_reducedSet_hierachical, m7loo )

#not converging with test error even for smaller set.
#seems like mean age is a better predictor than invSimp in the models

#Things to try: Loo comparison between mean age and diversity
#try removing tests as an error term
#high K doesn't change much - slightly worse fit (elpd_diff -0.2)

#fit %>% saveRDS("../results/f.spl.s2.nbin.ts_2021_04_30.rds")
#} else fit <- readRDS("../results/f.spl.s2.nbin.ts_2021_04_30.rds")

# check rhat diagnostics are close to 1: VERY GOOD
rhat(m7_reducedSet_hierachical) %>% unname %>% {c(min(.),max(.))}

# posterior predictive checks: NOT BAD
m7_reducedSet_hierachical%>% pp_check(nsamples = 1000) + xlim(c(NA,2e5)) + 
  labs(title = "Posterior predictive checks",
       subtitle = "Regularised spline with regional fixed effects on the shape (1000 draws)")

# example plot of smooths
conditional_smooths(m7_reducedSet_hierachical, spaghetti = T, nsamples = 500, smooths = c('s(spatLag,k=6,bs="ts")'), rug=T)

plot(p)+theme_bw()
#p+theme_bw() doesnt work yet. Have to manipulate the list

#all plots
m7_reducedSet_hierachical %>% conditional_smooths(spaghetti = T, nsamples = 500, rug=T)

#
m7_reducedSet_hierachical %>% conditional_effects(catergorical=T)

plot(fit )
summary(fit)

############################################################################################
#####################DEATH MODEL#####################
############################################################################################

source('prep_data_deaths.R')
cvd19deaths <- prepare_data_deaths (cvd19.raw ); cvd19deaths

# define covariates
vars_all <- setdiff(names(cvd19deaths),c("country","y")); vars_all
vars_numeric <- cvd19deaths %>% select(where(is.numeric),-y) %>% names; vars_numeric
vars_factor <- cvd19deaths %>% select(where(is.factor), -country) %>% names; vars_factor

#remove all variables with 0 slopes based on the complete model

#vars_reduced <- vars_numeric[-c(1,3,6:11, 13, 14:18,20:21)] 
vars_reduced <- vars_numeric[-c(3:5,8, 10:15, 17:21)] #

#herpes, trich, tb, density (?), rainfall (?), temperature, GDP, HCspend, invSimp (?/ hinge), hookworm, schisto, Lmphatic,
#IDbud, diabetes, cardion (?, )

# functions to specify models
mu_linear <- function(vars) paste("y ~",paste(vars, collapse = " + ")) %>% as.formula()

mu_spline <- function(vars_num, vars_fac, k = 6, bas = "tp"){  # N.B. bas = "ts" adds shrinkage
  paste("y ~", paste(vars_fac, collapse = " + "), " + ", 
        paste0("s(",vars_num,", k = ", k,", bs = '", bas, "')", collapse = " + ")) %>% 
    as.formula()
} 

#regional splines
mu_hierachical_spline <- function(vars_num, vars_fac, k = 6, bas = "tp"){ 
  paste(paste("y ~", paste(vars_fac, collapse = " + "), " + ",
              paste0("s(",vars_num,", by = reg, k=",k,", bs='", bas, "')", collapse="+"))) %>% 
    as.formula()
}

death_f.lin.s1.nbin<- bf(mu_linear(vars_all), shape ~ 1) + negbinomial() #not converging

f.lin.s2.nbin<- bf(mu_linear(vars_all), shape ~ 1 + reg) + negbinomial()

f.spl.s1.nbin.tp <- bf(mu_spline(vars_numeric, vars_factor), shape ~ 1) + negbinomial()

death_f.spl.s1.nbin.ts <- bf(mu_spline(vars_numeric, vars_factor, bas = "ts"), shape ~ 1) + negbinomial()


death_f.spl.s2.nbin.ts_reg <- bf(mu_spline(vars_numeric, vars_factor, bas = "ts"), shape ~ 1 + reg) + negbinomial()

death_m2.spl.Reduced <- bf(mu_spline(vars_reduced, vars_factor, bas = "ts"), shape ~ 1  ) + negbinomial()

death_m3.spl.Reduced.Reg <- bf(mu_spline(vars_reduced, vars_factor, bas = "ts"), shape ~ 1 +reg ) + negbinomial()

death_m3_spl.caseby.Reg <- as.formula(y ~ reg + s(cases, k = 6, bs = "ts", by=reg, id=1) + s(meanAge, k = 6, bs = "ts") + s(Malaria, k = 6, bs = "ts") + s(Trichuriasis, k = 6, bs = "ts") + s(popDen, k = 6, bs = "ts") + s(Hookworm, k = 6, bs = "ts") + s(spatLag, k = 6, bs = "ts"), 
                                      shape ~ 1 + 1|reg)

death_m4_spl.Reg_noCases <- as.formula(y ~ reg  + s(meanAge, k = 6, bs = "ts") + s(Malaria, k = 6, bs = "ts") + s(Trichuriasis, k = 6, bs = "ts") + s(popDen, k = 6, bs = "ts") + s(Hookworm, k = 6, bs = "ts") + s(spatLag, k = 6, bs = "ts"), 
                                       shape ~ 1 + 1|reg)
#SEM fit test

death_m4_spl.Reg_noCases <- as.formula(y ~ reg  + s(meanAge, k = 6, bs = "ts") + s(Malaria, k = 6, bs = "ts") + s(Trichuriasis, k = 6, bs = "ts") + s(popDen, k = 6, bs = "ts") + s(Hookworm, k = 6, bs = "ts") + s(spatLag, k = 6, bs = "ts"), 
                                       shape ~ 1 + 1|reg)

#--------------------------------------------------------
# brms fitting options
future = F
cores = 4
chains = 4
iter = 8000
control = list(adapt_delta = 0.99)
init_r = 1
#--------------------------------------------------------

#models

death_m1_full_lm <-  brm(death_f.lin.s1.nbin, data = cvd19deaths ,
                         chains = chains, iter = iter, future = future, control = control, cores = cores, 
                         init_r = init_r) #doesn't converge

death_m1_full_gam <-  brm(death_f.spl.s1.nbin.ts , data = cvd19deaths ,
                          chains = chains, iter = iter, future = future, control = control, cores = cores, 
                          init_r = init_r, save_all_pars=T) 
death_m1_full_gam  %>% saveRDS( 'death_m1_full_gam')

death_m1_full_gam_reg <-  brm(death_f.spl.s2.nbin.ts_reg , data = cvd19deaths ,
                              chains = chains, iter = iter, future = future, control = control, cores = cores, 
                              init_r = init_r, save_all_pars=T) #didnt converge
death_m1_full_gam_reg  %>% saveRDS( 'death_m1_full_gam')

death_m1_full_gam_reduced<-  brm(death_m2.spl.Reduced , data = cvd19deaths ,
                                 chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                 init_r = init_r) 
death_m1_full_gam_reduced %>% saveRDS('death_m1_full_gam_reduced')


death_m1_full_gam_reduced_reg<-  brm(death_m3.spl.Reduced.Reg  , data = cvd19deaths ,
                                     chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                     init_r = init_r) 
death_m1_full_gam_reduced_reg_by<-  brm(death_m3_spl.caseby.Reg , data = cvd19deaths ,
                                        chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                        init_r = init_r)#terrible fit
death_m1_full_gam_reduced_noCases<-  brm(death_m4_spl.Reg_noCases , data = cvd19deaths ,
                                         chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                         init_r = init_r) 
saveRDS(death_m1_full_gam_reduced_reg, 'death_m1_full_gam_reduced_reg')
load('death_m1_full_gam_reduced')


summary(death_m1_full_gam_reduced_reg )

# check rhat diagnostics are close to 1: VERY GOOD
rhat(death_m1_full_gam) %>% unname %>% {c(min(.),max(.))}

# posterior predictive checks: NOT BAD
death_m1_full_gam_reduced %>% pp_check(nsamples = 1000) + xlim(c(NA,2e5)) + 
  labs(title = "Posterior predictive checks",
       subtitle = "Regularised spline with regional fixed effects on the shape (1000 draws)")


