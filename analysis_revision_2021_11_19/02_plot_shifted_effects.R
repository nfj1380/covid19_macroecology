library(tidyverse)
library(pbmcapply)
library(ggpubr)
library(brms)

rm(list=ls())

theme_set(theme_classic())


plot_ce <- function(fit, X, CI = 0.9, shift = T, n_points = 100, n_draws = 1000, reg_level = NA, trans_level = NA){
  dat_0 = fit$data %>% as_tibble %>% select(!any_of(c("cases","deaths"))) %>% map_dfr(mean)
  
  X_vals <- fit$data %>% pull(X) %>% {seq(min(.),max(.),length.out = n_points)}
  X_vals <- fit$data %>% pull(X) %>% {seq(min(.),quantile(.,0.95),length.out = n_points)}
  
  X_mean <- fit$data %>% pull(X) %>% mean
  X_ref <- abs(X_vals-X_mean) %>% which.min()
  
  dat_new = X_vals %>% 
    map_dfr( ~ dat_0 %>% mutate(across(all_of(X), function(v) v = .x))) %>% 
    mutate(reg = reg_level, trans = trans_level)
  ep <- fit %>% posterior_epred(newdata = dat_new, ndraws = n_draws, allow_new_levels = T) %>% t
  colnames(ep) <- paste0("s",1:ncol(ep))
  
  pred_tib = ep %>% as_tibble %>% mutate(x = X_vals) %>% 
    pivot_longer(-x, names_to = "sample", values_to = "y")
  
  
  if(shift) pred_tib = pred_tib %>% group_by(sample) %>% mutate(y = y - y[X_ref]) %>% ungroup()
  
  pred_data = pred_tib %>%
    group_by(x) %>% 
    summarise(med = quantile(y,0.5), lwr = quantile(y,(1-CI)/2), upr = quantile(y,(1+CI)/2))
  
  y_max = min(max(pred_data$upr),1e7)
  
  pred_data %>% 
    ggplot(aes(x=x)) +
    geom_line(aes(y=med), size = 1, col = blues9[8]) +
    geom_line(aes(y = upr), lty = "dashed", col = blues9[8]) +
    geom_line(aes(y = lwr), lty = "dashed", col = blues9[8]) +
    geom_hline(aes(yintercept = 0)) +
    labs(x = X, y = resp) +
    ylim(c(NA,y_max))
  
}


resp <- "cases"
run_date <- "2021_11_24"
fit <- readRDS(paste0("results/fit_vars_all_",resp,"_spline_",run_date,".rds"))
fit
#fit %>% pp_check(ndraws = 200) + xlim(c(0,1.5e6))

vars <- fit$data %>% select(where(is.numeric), -any_of(c("cases","deaths")), -pop) %>% names
ce_plots <- pbmclapply(vars, function(X) plot_ce(fit,X), mc.cores = length(vars))

main_plot <- ggarrange(plotlist = ce_plots, nrow = 5, ncol = 4) %>% 
  annotate_figure(top = text_grob(resp, size = 16, face = "bold"), bottom = paste0("CI: 90%"))
ggsave(paste0("plots/",resp,"_all_cond_effects_",run_date,".pdf"), width = 18, height = 12, device = cairo_pdf)

                     