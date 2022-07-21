library(tidyverse)
library(pbmcapply)
library(ggpubr)
library(brms)

rm(list=ls())

theme_set(theme_classic())

source("00_functions.R")

resp <- "cases"
run_date <- "2021_11_29"
fit <- readRDS(paste0("results/fit_vars_all_",resp,"_spline_",run_date,".rds"))
fit <- readRDS(paste0("results/fit_subset1_",resp,"_tp_",run_date,".rds"))
fit
#fit %>% pp_check(ndraws = 200) + xlim(c(0,1.5e6))

vars <- fit$data %>% select(where(is.numeric), -any_of(c("cases","deaths")), -pop) %>% names

plot_ce(fit_cases, "Hookworm")

ce_plots <- pbmclapply(vars, function(X) plot_ce(fit,X), mc.cores = length(vars))
ce_plots[[10]]
main_plot <- ggarrange(plotlist = ce_plots, nrow = 3, ncol = 3) %>% 
  annotate_figure(top = text_grob(resp, size = 16, face = "bold"), bottom = paste0("CI: 90%"))

ggsave(paste0("plots/",resp,"_sub1_cond_effects_",run_date,".pdf"), width = 9, height = 9, device = cairo_pdf)

                     