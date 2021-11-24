## imputation attempts
# imputation appears not work for splines
# works well for linear models with values missing from multiple predictors. 

pacman::p_load("tidyverse","brms")


rm(list=ls())

source("00_functions.R")

# load data set
cvd19_raw <- read_csv("data/cvd19_data_2021_11_19_no_imp.csv")
cvd19 <- cvd19_raw %>% prepare_data()

cvd19_raw
cvd19

run_date <- "2021_11_19"

missing_vars <- cvd19 %>% map_dbl(~ .x %>% is.na %>% sum) %>% {.[.!=0]} %>% names

vars_cases <- c("HIV.AIDS", "Ascariasis", "Malaria", "PerUrb", "logTests", "spatLag", "trans")


## try a linear model with variable imputation
f_c_linear_mu <- bf(paste("y_c ~ ",
                    paste(vars_cases[!vars_cases %in% missing_vars], collapse = " + "),
                    "+ ",
                    paste(paste0("mi(",vars_cases[vars_cases %in% missing_vars],")"),collapse = " + "))) + Beta()

f_c_linear_mi_1 <- paste0(vars_cases[vars_cases %in% missing_vars][1],"|mi()~ ", 
       paste(vars_cases[!vars_cases %in% missing_vars], collapse = " + ")) %>% bf + Beta()

f_c_linear_mi_2 <- paste0(vars_cases[vars_cases %in% missing_vars][2],"|mi()~ ", 
                          paste(vars_cases[!vars_cases %in% missing_vars], collapse = " + ")) %>% bf + lognormal()

f_c_linear <- f_c_linear_mu + f_c_linear_mi_1 + f_c_linear_mi_2 + set_rescor(FALSE)
  
f_c_linear

m.linear <- brm(f_c_linear, data = cvd19, cores = 4, control = list(adapt_delta=0.85), init = 0)
# terrible convergence






## try a weighted poisson regression
f_cases_mu <- cases ~ 1 + reg + 
  trans + 
  s(HIV.AIDS, k = 6, bs = "ts") + 
  s(Ascariasis, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(PerUrb, k = 6, bs = "ts") + 
  s(logTests, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts") 






f_c <- 
  bf(y_c ~ 1 + 
       + 
       mi(PerUrb)) + Beta() +
  bf(PerUrb|mi() ~ 1 + Ascariasis) + Beta()
  
m.2 <- brm(f_c, data = cvd19, cores = 4, control = list(adapt_delta=0.9))
m.2

trans + 
  s(HIV.AIDS, k = 6, bs = "ts") + 
  s(Ascariasis, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(PerUrb, k = 6, bs = "ts") + 
  s(tests, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts") 



# tests, transmission are missing 
f_cases_mu <- bf(y_c ~ 1 + reg + 
  #trans + 
  s(HIV.AIDS, k = 6, bs = "ts") + 
  s(Ascariasis, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(PerUrb, k = 6, bs = "ts") + 
  #s(tests, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts")) + Beta()





m.2$data %>% dim

fit <- m.2
data <- fit$data %>% mutate(country = cvd19$country) %>% rename(y = y_c)
f_epred <- fit %>% posterior_predict
f_epred 
fit$data

colnames(f_epred) <- data$country
p <- f_epred %>% 
  as_tibble() %>% 
  mutate(across(everything())) %>% 
  pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
  group_by(country) %>% 
  summarise(q025 = quantile(y_pred, probs = 0.025),
            q975 = quantile(y_pred, probs = 0.975),
            q25 = quantile(y_pred, probs = 0.25),
            q75 = quantile(y_pred, probs = 0.75),
            mean = median(y_pred)) %>% 
  arrange(country) %>% 
  mutate(y_obs = data %>% arrange(country) %>% pull(y)) %>% 
  mutate(reg = data %>% arrange(country) %>% pull(reg)) %>% 
  arrange(reg,mean) %>% 
  ggplot(aes(x = factor(country, levels = country))) +
  geom_linerange(aes(ymin = q025, ymax = q975, col = reg), alpha = 0.7) +
  geom_linerange(aes(ymin = q25, ymax = q75, col = reg), size = 1.5, alpha = 0.7) +
  #geom_point(aes(y = mean, col = reg), alpha = 1, shape = 18, size = 2) +
  geom_line(aes(y = mean, group = reg), col = "white") +
  geom_point(aes(y = y_obs), size =0.8) +
  theme_classic() +
  theme(legend.position = "top") +
  scale_color_brewer(type = "div", palette = "Dark2") +
  guides(color = guide_legend(nrow = 1, byrow = TRUE, title = NULL))
p

m.1 %>% pp_check(ndraws = 200)
m.1 %>% plot
zero_inflated_beta()