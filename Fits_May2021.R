pacman::p_load("tidyverse","brms") # mice, future

rm(list=ls())

# plan(multisession, workers = 3) # for use with future

# load functions and data
source("fns_2021_04_30.R") 
data.raw <- read_csv("cvd19_imp_2021_04_30.csv")
cvd19 <- prepare_data(data.raw); cvd19


############################################################################################
#####################CASE MODEL#####################

# define covariates
vars_all <- setdiff(names(cvd19),c("country","y")); vars_all
vars_numeric <- cvd19 %>% select(where(is.numeric),-y) %>% names; vars_numeric
vars_factor <- cvd19 %>% select(where(is.factor), -country) %>% names; vars_factor

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


m1_full <-  brm(f.spl.s1.nbin.ts, data = cvd19,
                chains = chains, iter = iter, future = future, control = control, cores = cores, 
                init_r = init_r) 
  
m2_full_WHO <-  brm(f.spl.s2.nbin.ts, data = cvd19,
                      chains = chains, iter = iter, future = future, control = control, cores = cores, 
                      init_r = init_r) 

m3_reducedSet_WHO <-  brm(f.spl.Reduced , data = cvd19,
                    chains = chains, iter = iter, future = future, control = control, cores = cores, 
                    init_r = init_r) 

#m4_reducedSet_WHO_Tests <-  brm(f.spl.Reduced , data = cvd19,
                         # chains = chains, iter = iter, future = future, control = control, cores = cores, 
                          #init_r = init_r) #doesn't converge

m4_reducedSet_WHO <-  brm(f.spl.Reduced , data = cvd19,
                                chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                init_r = init_r) 
m4_reducedSet_WHO %>% saveRDS("m4_reducedSetInSimp_WHO .rds")

m5_reducedSet_WHO_highK <-  brm(f.spl.Reduced , data = cvd19,
                                chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                init_r = init_r) 

#m6_reducedSetInSimp_tests <-  brm(f.spl.Reduced_tests  , data = cvd19,
                                      #chains = chains, iter = iter, future = future, control = control, cores = cores, 
                                      #init_r = init_r)  #not converging even without WHO

m7_reducedSet_hierachical <-  brm(f.spl.Reduced_hierachical , data = cvd19,
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

source('prep_data_deaths.R')
cvd19deaths <- prepare_data_deaths (data.raw); cvd19deaths

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




#################################################################################
#FINAL models from supercomputer
#################################################################################
library(brms)
case_FINAL <- readRDS('m.cases_2021_05_18.rds')
summary(case_FINAL)

multi_model <- readRDS('m_fit_1_2021_05_28.rds')
summary(multi_model)

#define limits for some variables
lim_cases <- quantile(cvd19deaths$cases, 0.95)
lim_age[1] <- quantile(cvd19deaths$meanAge, 0.95)
lim_age_lower <- quantile(cvd19deaths$meanAge, 0.05)
lim_den <- quantile(cvd19deaths$popDen, 0.95)
lim_spat <- quantile(cvd19deaths$spatLag, 0.95)
# example plot of smooths
#gg <- conditional_smooths(death_m1_full_gam_reduced_reg  , spaghetti = T, nsamples = 500, rug=T) #marginal is the same as conditional
gg <- conditional_effects(multi_model, effect=c('cases') ,resp='deaths', spaghetti = T, nsamples = 500) #%>% plot()
gg <- conditional_effects(multi_model,  spaghetti = T, nsamples = 500) #%>% plot()

plot(gg)
bayes_R2(case_FINAL)
bayes_R2(death_m1_full_gam_reduced_reg)
#gg<- conditional_effects(death_m1_full_gam_reduced_reg ,  int_conditions = list(reg = mean), select_points=10)

a_plot <- plot(gg, plot = FALSE, rug=T)[[1]] + 
  theme_bw()#+
 #ylim(0,20000)

b_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[2]] + 
  theme_bw()#+
  #xlim(0,lim_cases)+
#geom_rug(data=cvd19deaths$cases,sides="b")
#b_plot$layers$aes_params$colour <- 'black' 
 #ylim(0,20000)

c_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[3]] + 
  theme_bw()#+
  #xlim(lim_age_lower,lim_age)+
 #ylim(0,20000)


d_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[4]] + 
  theme_bw()#+
  #ylim(0,20000)


e_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[5]] + 
  theme_bw()#+
 #ylim(0,20000)

f_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[6]] + 
  theme_bw()#+
 # xlim(0,lim_den)+
  #ylim(0,20000)

g_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[7]] + 
  theme_bw()#+
  #ylim(0,20000)

h_plot <- plot(gg, plot = FALSE, rug=T, spaghetti_args=c(colour='gray85', alpha=0.01),
               line_args = c(colour = "black"))[[8]] + 
  theme_bw()#+
 # xlim(0, lim_spat)+
  #ylim(0,20000)

gridExtra::grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot,f_plot,g_plot,h_plot, nrow = 2)
gridExtra::grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot,f_plot, g_plot, nrow = 2)
gridExtra::grid.arrange(a_plot, b_plot, c_plot, d_plot, e_plot,f_plot,  nrow = 2)
  #geom_rug(data=data.frame(cvd19deaths$meanAge), sides="b") 

death_m1_loo <- loo(death_m1_full_gam) 
death_m1_loo_reg <- loo(death_m1_full_gam_reg) 
death_m1_loo_reduced <- loo(death_m1_full_gam_reduced, moment_match = TRUE) 
death_m1_loo_reduced_reg <- loo(death_m1_full_gam_reduced_reg) 
death_m1_loo_reduced_reg_by<- loo(death_m1_full_gam_reduced_reg_by, moment_match = TRUE)
#need to set 'save_all_pars'
plot(p)+theme_bw()

loo_compare( death_m1_loo , death_m1_loo_reduced_reg, death_m1_loo_reg)

f_epred <- case_FINAL  %>% rstantools::posterior_predict()
colnames(f_epred) <- cvd19$country
f <- f_epred %>% 
  as_tibble() %>% 
  pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
  group_by(country) %>% 
  summarise(q025 = quantile(y_pred, probs = 0.025),
            q975 = quantile(y_pred, probs = 0.975),
            mean = mean(y_pred)) %>% 
  arrange(country) %>% 
  mutate(y_obs = cvd19  %>% arrange(country) %>% pull(y)) %>% 
  arrange(y_obs) %>% 
  ggplot(aes(x = factor(country, levels = .$country))) +
  geom_point(aes(y = y_obs)) +
  geom_point(aes(y = mean), col = "blue") +
  geom_linerange(aes(ymin = q025, ymax = q975), col = "darkblue") +
  theme_classic() +
  labs(x = "Country (order by cases)", y = "cases/1000", 
       title = "Posterior predictive 95% credible intervals + data",
       subtitle = "data = black, model predictions = blue") +
 #scale_x_discrete(labels = )+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

saved_Preds <- f_epred %>% 
  as_tibble() %>% 
  pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
  group_by(country) %>% 
  summarise(q025 = quantile(y_pred, probs = 0.025),
            q975 = quantile(y_pred, probs = 0.975),
            mean = mean(y_pred)) %>% 
  arrange(country) %>% 
  mutate(y_obs = cvd19  %>% arrange(country) %>% pull(y)) %>% 
  arrange(y_obs)
################################################################################################
#IML testing
################################################################################################

fover <- f %>% filter(y_obs > q975)

library(iml)
X <- cvd19[which(names(cvd19) != "y")]

str(X)
predictor <- Predictor$new(case_FINAL, data = X, y = cvd19$y)

#effs <- FeatureEffects$new(predictor)

shapley <- Shapley$new(predictor, x.interest = X[172, ])
shapley$plot()



# extra stuff
if(F){
  # does not converge well - large rhat values
  m.lin.s1.nbin.no.5 <-  brm(f.lin.s1.nbin, data = cvd19, prior = set_prior("normal(0,0.5)", class = "b"),
                             chains = chains, iter = iter, future = future, control = control, cores = cores)
  
  # somewhat better using stronger regularisation
  m.lin.s1.nbin.hs1 <-  brm(f.lin.s1.nbin, data = cvd19, prior = set_prior("horseshoe(1)", class = "b"),
                            chains = chains, iter = iter, future = future, control = control, cores = cores, 
                            init_r = init_r)
  
  # will splines help? The problem may be too many modelled degrees of freedom relative to size of the data set.
  # Although mispecification may be the issue
  m.spl.s1.nbin.hs3 <-  brm(f.spl.s1.nbin, data = cvd19,
                            chains = chains, iter = iter, future = future, control = control, cores = cores, 
                            init_r = init_r)
  

}
