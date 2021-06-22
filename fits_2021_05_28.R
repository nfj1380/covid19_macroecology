pacman::p_load("tidyverse","brms") # mice, future

rm(list=ls())

# load functions and data
source("fns_2021_05_28.R") 
data.raw <- read_csv("../data/cvd19_imp_2021_04_30.csv")
cvd19 <- prepare_data(data.raw); cvd19

# reduced covariate set
f_cases_mu <- cases ~ 1 + reg + 
  trans + 
  s(HIV.AIDS, k = 6, bs = "ts") + 
  s(Ascariasis, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(PerUrb, k = 6, bs = "ts") + 
  s(tests, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts") 

f_death_mu <- deaths ~ 1 + reg + 
  s(cases, k = 6, bs = "ts") + 
  s(meanAge, k = 6, bs = "ts") + 
  s(Malaria, k = 6, bs = "ts") + 
  s(Trichuriasis, k = 6, bs = "ts") + 
  s(popDen, k = 6, bs = "ts") + 
  s(Hookworm, k = 6, bs = "ts") + 
  s(spatLag, k = 6, bs = "ts") 

f_shape <- shape ~ 1 + 1|reg

form_cases <- bf(f_cases_mu,f_shape) + negbinomial()
form_death <- bf(f_death_mu,f_shape) + negbinomial()

# sampler params
seed <- 687756 # sample(1e6,1)
control = list(adapt_delta = 0.9)

# fit and save
m_fit_1 <- brm(form_cases + form_death + set_rescor(F), data = cvd19, cores = 4, chains = 4, 
             init_r = 1, iter = 4000, thin = 8, seed = seed, control = control)
#saveRDS(m_fit_1,"m_fit_1_2021_05_28.rds")

m_fit_2 <- brm(form_cases + form_death + set_rescor(F), data = cvd19, cores = 4, chains = 4, 
             init_r = 1, iter = 4000, thin = 8, seed = seed, control = control)

m_fit_1 <- m_fit


#saveRDS(m.cases,"m.cases_2021_05_18.rds")


# check rhat diagnostics are close to 1: VERY GOOD
rhat(fit) %>% unname %>% {c(min(.),max(.))}

fit %>% plot(pars = "^b_")

# posterior predictive checks: NOT BAD
fit %>% pp_check(nsamples = 1000) + xlim(c(NA,2e5)) + 
  labs(title = "Posterior predictive checks",
       subtitle = "Regularised spline with regional fixed effects on the shape (1000 draws)")

# example plot of smooths
fit %>% conditional_smooths(spaghetti = T, nsamples = 500, smooths = c('s(meanAge,k=3,bs="ts")')) 

f_epred <- fit %>% posterior_epred()
colnames(f_epred) <- cvd19$country
f_epred %>% 
  as_tibble() %>% 
  pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
  group_by(country) %>% 
  summarise(q025 = quantile(y_pred, probs = 0.025),
            q975 = quantile(y_pred, probs = 0.975),
            mean = mean(y_pred)) %>% 
  arrange(country) %>% 
  mutate(y_obs = cvd19 %>% arrange(country) %>% pull(y)) %>% 
  arrange(y_obs) %>% 
  ggplot(aes(x = factor(country, levels = country))) +
  geom_point(aes(y = y_obs)) +
  geom_point(aes(y = mean), col = "blue") +
  geom_linerange(aes(ymin = q025, ymax = q975), col = "darkblue") +
  theme_classic() +
  labs(x = "Country (order by observed cases)", y = "cases/1000", 
       title = "Posterior predictive 95% credible intervals + data",
       subtitle = "data = black, model predictions = blue") +
  scale_x_discrete(labels = NULL)
  
  

