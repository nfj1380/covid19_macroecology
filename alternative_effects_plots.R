
#---------------------------------
# alternative conditional effects
#---------------------------------

pacman::p_load("tidyverse","brms","ggpubr","RColorBrewer","bayesplot")

rm(list=ls())

# load cases
m.cases <- readRDS("results/fit_cases_2021_06_24.rds")
m.deaths <- readRDS("results/fit_deaths_2021_06_24.rds")

# load data
cvd19 <- read_csv("data/cvd19_2021_06_24.csv")

m.cases <- readRDS("results/fit_cases_2021_06_24.rds")

# default brms plot
plot_ce_brms <- m.cases %>% conditional_effects("Ascariasis")


# range
vals_Asc <- seq(0,0.35,length.out = 40)

# index of value within the range at which to fix conditional smooths
ref_index <- 1 # vals_Asc[ref_index] (6 ~ mean)

d_cond <- m.cases$data %>% 
  group_by(reg,trans) %>% 
  mutate(across(-all_of(c("cases")),mean)) %>% 
  ungroup() 
  mutate(reg = NA, trans = NA)

d_cond_Asc <- vals_Asc %>% map_dfr(~ mutate(d_cond, Ascariasis = .x)) %>% mutate(trans = NA)

# as per brms
d_cond <- m.cases$data %>% mutate(reg = factor(reg)) %>% mutate(across(-all_of(c("cases","reg","trans")),mean)) %>% 
  mutate(reg = levels(reg)[1], trans = levels(factor(trans))[1]) %>% slice(1)

d_cond_Asc <- vals_Asc %>% map_dfr(~ mutate(d_cond, Ascariasis = .x)) %>% 
  mutate(trans = factor(trans, levels = c("clust","comm"))) %>% 
  mutate(reg = factor(reg, levels = c("Afri","Amer")))

ce_Asc <- posterior_epred(m.cases, nsamples = 500, newdata = d_cond_Asc)


# brms-type plots
ce_Asc_data <- ce_Asc %>% t %>% as_tibble() %>% mutate(x = d_cond_Asc$Ascariasis) %>% 
  group_by(x) %>% summarise(across(everything(),mean)) %>% 
  ungroup() %>% 
  pivot_longer(-x, values_to = "y", names_to = "sample") 

y_min <- ce_Asc_data %>% group_by(x) %>% summarise(y = quantile(y,0.05)) %>% pull(y) %>% min
y_max <- ce_Asc_data %>% group_by(x) %>% summarise(y = quantile(y,0.95)) %>% pull(y) %>% max
y_lims <- c(y_min-(y_max - y_min)/10,y_max+(y_max - y_min)/10)


plot_ce_brms_spaghetti <-ce_Asc_data %>% 
  ggplot(aes(x,y)) +
  geom_line(aes(group = sample), alpha = 0.2, col = blues9[7]) +
  geom_line(size = 1, lty = "solid", data = ~ .x %>% group_by(x) %>% summarise(y = mean(y))) +
  geom_line(size = 1, lty = "dashed", data = ~ .x %>% group_by(x) %>% summarise(y = quantile(y,0.05))) +
  geom_line(size = 1, lty = "dashed", data = ~ .x %>% group_by(x) %>% summarise(y = quantile(y,0.95))) +
  theme_minimal() +
  ylim(y_lims) +
  #geom_vline(aes(xintercept = vals_Asc[6]), col = "grey20", lty = "dashed") +
  labs(title = "Ascariasis: conditional effects",  y = "cases", x = "Ascariasis")


plot_ce_brms <-ce_Asc_data %>% 
  ggplot(aes(x,y)) +
  geom_ribbon(aes(y = NULL, ymin = ymin, ymax = ymax),
              size = 1, fill = "grey70",
              data = ~ .x %>% group_by(x) %>% summarise(ymin = quantile(y,0.05), ymax = quantile(y,0.95))) +

  theme_minimal() +
  geom_line(size = 1, lty = "solid", data = ~ .x %>% group_by(x) %>% summarise(y = mean(y))) +
  ylim(y_lims) +
  #geom_vline(aes(xintercept = vals_Asc[6]), col = "grey20", lty = "dashed") +
  labs(title = "Ascariasis: conditional effects",  y = "cases", x = "Ascariasis")

# shifted plots
ce_Asc_data <- ce_Asc %>% t %>% as_tibble() %>% mutate(x = d_cond_Asc$Ascariasis) %>% 
  group_by(x) %>% summarise(across(everything(),mean)) %>% 
  ungroup() %>% 
  mutate(across(-x,~ .x - .x[ref_index])) %>% 
  pivot_longer(-x, values_to = "y", names_to = "sample") 
  
y_min <- ce_Asc_data %>% group_by(x) %>% summarise(y = quantile(y,0.05)) %>% pull(y) %>% min
y_max <- ce_Asc_data %>% group_by(x) %>% summarise(y = quantile(y,0.95)) %>% pull(y) %>% max
y_lims <- c(y_min-(y_max - y_min)/10,y_max+(y_max - y_min)/10)


plot_ce_shifted <-ce_Asc_data %>% 
  ggplot(aes(x,y)) +
  geom_hline(aes(yintercept = 0), col = "grey20", lty = "dashed") +
  geom_line(aes(group = sample), alpha = 0.2, col = blues9[7]) +
  geom_line(size = 1, lty = "solid", data = ~ .x %>% group_by(x) %>% summarise(y = mean(y))) +
  geom_line(size = 1, lty = "dashed", data = ~ .x %>% group_by(x) %>% summarise(y = quantile(y,0.05))) +
  geom_line(size = 1, lty = "dashed", data = ~ .x %>% group_by(x) %>% summarise(y = quantile(y,0.95))) +
  theme_minimal() +
  ylim(y_lims) +
  #geom_vline(aes(xintercept = vals_Asc[ref_index]), col = "grey20", lty = "dashed") +
  labs(title = "Ascariasis: shifted conditional effects",  y = "change in cases", x = "Ascariasis")

ggarrange(plot_ce_brms,plot_ce_brms_spaghetti,plot_ce_shifted, ncol =1)
ggsave("plots/ce_plots.pdf", width = 4, height = 8)
