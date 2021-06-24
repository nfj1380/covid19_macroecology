pacman::p_load("tidyverse","brms","ggpubr","RColorBrewer","bayesplot")

rm(list=ls())

# load cases
fit_cases <- readRDS("results/fit_cases_2021_06_24.rds")
fit_deaths <- readRDS("results/fit_deaths_2021_06_24.rds")

# load data
cvd19 <- read_csv("data/cvd19_2021_06_24.csv")

# check rhat diagnostics are close to 1: VERY GOOD (Rhat < 1.007)
rhat(fit_cases) %>% unname %>% {c(min(.),max(.))}
rhat(fit_deaths) %>% unname %>% {c(min(.),max(.))}

# posterior predictive densities
pp <- list()
pp$cases <- pp_check(fit_cases, nsamples = 300) + xlim(c(0,2e5)) + labs(title = " ")
pp$deaths <- pp_check(fit_deaths, nsamples = 300) + xlim(c(0,3e3)) + labs(title = " ")
pp %>% ggarrange(plotlist = ., nrow = 1, labels = "AUTO", common.legend = T, legend = "bottom")
ggsave("plots/cvd19_pp_dens_2021_04_24.png") # Supp Fig XX

# posterior predictive intervals

# check data is in the correct order: YES
round(cvd19$cases) - fit_cases$data$cases
round(cvd19$deaths) - fit_deaths$data$deaths

plot_pred <- function(model, xlabels = F){
  if(model == "cases") fit <- fit_cases
  if(model == "deaths") fit <- fit_deaths
  data <- fit$data %>% mutate(country = cvd19$country) %>% rename(y = any_of(model))
  f_epred <- fit %>% posterior_predict()
  colnames(f_epred) <- data$country
  p <- f_epred %>% 
    as_tibble() %>% 
    mutate(across(everything(),log)) %>% 
    pivot_longer(everything(), names_to = "country", values_to = "y_pred") %>% 
    group_by(country) %>% 
    summarise(q025 = quantile(y_pred, probs = 0.025),
              q975 = quantile(y_pred, probs = 0.975),
              q25 = quantile(y_pred, probs = 0.25),
              q75 = quantile(y_pred, probs = 0.75),
              mean = median(y_pred)) %>% 
    arrange(country) %>% 
    mutate(y_obs = data %>% arrange(country) %>% pull(y) %>% log) %>% 
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
  
  if(!xlabels) p <- p + theme(axis.ticks.x = element_blank()) + scale_x_discrete(labels = NULL)
  if(xlabels) p <- p + theme(axis.text.x = element_text(angle = 90, size = 7, hjust = 1) )
  
  p
}

# without country labels for main manuscript
ggarrange(plot_pred("cases") + labs(y = "log(cases)/million", x = NULL),
          plot_pred("deaths") + labs(y = "log(deaths)/million", 
                                     x = "Country (by increasing predicted response)"), nrow = 2,
          common.legend = T, legend = "bottom")
ggsave("plots/cvd19_plot_pp_predict_2021_06_24.png")


# with country labels for supplementary materials (unfinished)
plot_pred("cases", xlabels = T) # need to abbreviate country names and possible split plot into two.
plot_pred("deaths", xlabels = T)


# plot region-level shape estimates
# the estimated shape parameter is log(phi) where the associated
# variance is: var = mu + mu^2/phi
plot_shape_est <- function(fit){
  mcmc_shape_reg <- fit %>% as.data.frame() %>% select(contains("reg__sha")) %>% 
    mutate(across(everything(), ~ .x + sd_reg__shape_Intercept)) %>% 
    select(-sd_reg__shape_Intercept)
  reg_levels <- fit$data$reg %>% factor %>% levels
  colnames(mcmc_shape_reg) <- reg_levels
  mcmc_shape_reg %>% 
    #mutate(across(everything(), ~ .x %>% {1/.})) %>% 
    mcmc_intervals(pars = rev(reg_levels),
                   prob_outer = 0.95,
                   point_est = "median",
                   point_size = 5) + 
    coord_flip() +
    labs(x = "shape intercept (log-link scale) ")
}

ggarrange(plot_shape_est(fit_cases) + xlim(c(-2.1,3)),
          plot_shape_est(fit_deaths) + xlim(c(-2.1,3)) + labs(x = NULL), nrow = 1, labels= "AUTO")

ggsave("plots/cvd19_shape_effects_2021_06_24.png", width = 6, height = 3)




