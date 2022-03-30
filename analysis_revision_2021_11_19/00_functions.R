# plot conditional effects
plot_ce <- function(fit, X, CI = 0.9, shift = T, n_points = 100, n_draws = 1000, reg_level = NA, trans_level = NA, exp = F){
  dat_0 = fit$data %>% as_tibble %>% select(!any_of(c("cases","deaths"))) %>% map_dfr(mean)
  
  X_vals <- fit$data %>% pull(X) %>% {seq(min(.),max(.),length.out = n_points)}
  X_vals <- fit$data %>% pull(X) %>% {seq(min(.),quantile(.,0.95),length.out = n_points)}
  if(exp) X_vals <- fit$data %>% pull(X) %>% exp %>%  {seq(min(.),quantile(.,0.95),length.out = n_points)} %>% log

  
  X_mean <- fit$data %>% pull(X) %>% mean
  if(exp) X_mean <- fit$data %>% pull(X) %>% exp %>% mean %>% log
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
  
  if(exp) pred_data <- pred_data %>% mutate(x = exp(x))
  if(exp) X <- paste0("exp(",X,")")
  
  pred_data %>% 
    ggplot(aes(x=x)) +
    geom_line(aes(y=med), size = 1, col = blues9[8]) +
    geom_line(aes(y = upr), lty = "dashed", col = blues9[8]) +
    geom_line(aes(y = lwr), lty = "dashed", col = blues9[8]) +
    geom_hline(aes(yintercept = 0)) +
    labs(x = X, y = resp) +
    ylim(c(NA,y_max))
  
}
