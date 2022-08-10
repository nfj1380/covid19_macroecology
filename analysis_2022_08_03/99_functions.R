
# construct spline-based model formula for the mean
mu_spline <- function(response, vars_num, vars_fac, vars_end, k_num = 6, k_end = 6){
  paste(response," | rate(pop) ~ 1 + (1|a|reg) + ", paste(vars_fac, collapse = " + "), " + ", 
        paste0("s(",vars_num,", k = ", k_num,", bs = 'ts')", collapse = " + ")," + ", 
        paste0("s(",vars_end,", k = ", k_end,", bs = 'ts')", collapse = " + ")) %>% 
    as.formula()
} 


# plot conditional effects
plot_ce <- function(fit, X, resp, CI = 0.9, shift = T, n_points = 100, n_draws = 1000, 
                    reg_level = NA, trans_level = NA, exp = F){
  dat_0 = fit$data %>% as_tibble %>% select(where(is.numeric)) %>% map_dfr(mean)
  dat_0 <- dat_0 %>% mutate(trans = "clust")
  
  #X_vals <- fit$data %>% pull(X) %>% {seq(min(.),max(.),length.out = n_points)}
  X_vals <- fit$data %>% pull(X) %>% {seq(quantile(.,0.05),quantile(.,0.95),length.out = n_points)}
  if(exp) X_vals <- fit$data %>% pull(X) %>% exp %>%  {seq(min(.),quantile(.,0.95),length.out = n_points)} %>% log

  X_mean <- fit$data %>% pull(X) %>% mean
  if(exp) X_mean <- fit$data %>% pull(X) %>% exp %>% mean %>% log
  X_ref <- abs(X_vals-X_mean) %>% which.min() # find closest value to the mean

  dat_new = X_vals %>% map_dfr( ~ dat_0 %>% mutate(across(all_of(X), function(v) v = .x)))
  
  # to predict to specific regions or transmission types 
  if(!is.na(reg_level)) dat_new %>% mutate(reg = reg_level)
  if(!is.na(trans_level)) dat_new %>% mutate(trans = trans_level)
  
  ep <- fit %>% posterior_epred(newdata = dat_new, ndraws = n_draws, re_formula = NA) %>% t
  colnames(ep) <- paste0("s",1:ncol(ep))
  
  pred_tib = ep %>% as_tibble %>% mutate(x = X_vals) %>% 
    pivot_longer(-x, names_to = "sample", values_to = "y")
  
  if(shift) pred_tib = pred_tib %>% group_by(sample) %>% mutate(y = y - y[X_ref]) %>% ungroup()
  
  pred_data = pred_tib %>%
    group_by(x) %>% 
    summarise(med = quantile(y,0.5), lwr = quantile(y,(1-CI)/2), upr = quantile(y,(1+CI)/2))
  
  y_max = min(max(pred_data$upr),1e7)
  
  if(exp) pred_data <- pred_data %>% mutate(x = exp(x))
  if(exp) X <- paste0("exp(",X,")")
  
  pred_data %>% 
    ggplot(aes(x=x)) +
    geom_line(aes(y=med), size = 1, col = blues9[8]) +
    geom_line(aes(y = upr), lty = "dashed", col = blues9[8]) +
    geom_line(aes(y = lwr), lty = "dashed", col = blues9[8]) +
    geom_hline(aes(yintercept = 0)) +
    labs(x = X, y = paste("Change in",resp)) +
    ylim(c(NA,y_max))
  
}
