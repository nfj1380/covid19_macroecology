
# construct spline-based model formula for the mean
#ts shrings linear and non-linear component, tp only shrings the non-linear
mu_spline_all <- function(response,  vars_numeric, vars_endemic, k_num = 6, k_end = 6){
  paste(response,"| se(pscore.mean.se, sigma = TRUE)  ~ 1 + (1|a|reg)+trans", " + ", 
        paste0("s(",vars_numeric,", k = ", k_num,", bs = 'ts')", collapse = " + ")," + ", 
        paste0("s(",vars_endemic,", k = ", k_end,", bs = 'ts')", collapse = " + ")) %>% 
    as.formula()
} 
mu_spline_all_sp <- function(response,  vars_numeric, vars_endemic, k_num = 6, k_end = 6){
  paste(response," | se(pscore.mean.se) ~ car(W, type='icar', gr = NA)", " + ", 
        paste0("s(",vars_numeric,", k = ", k_num,", bs = 'ts')", collapse = " + ")," + ", 
        paste0("s(",vars_endemic,", k = ", k_end,", bs = 'ts')", collapse = " + ")) %>% 
    as.formula()
} 

mu_spline_subModel <- function(response,  vars_numeric, vars_factor, k_num = 6){
  paste(response,"| se(pscore.mean.se, sigma = TRUE)  ~ 1 + (1|a|reg)+trans", " + ", 
        paste0("s(",vars_numeric,", k = ", k_num,", bs = 'ts')", collapse = " + ")) %>% 
    as.formula()
} 

mu_linear <- function(response,  vars_numeric, vars_endemic, k_num = 6, k_end = 6){
  paste(response,"| se(pscore.mean.se)  ~  trans", " + ", 
        #paste0(vars_endemic, collapse = " + ")," + ", 
        paste0(vars_numeric, collapse = " + ")) %>% 
    as.formula()
}
mu_spline_by <- function(response, offset, vars_numeric, vars_factor, vars_endemic, k_num = 6, k_end = 6){
  paste(response," | rate(",offset,") ~ 1 + (1|a|reg) + ", paste(vars_factor, collapse = " + "), " + ", 
        paste0("s(",vars_numeric,", k = ", k_num,", bs = 'ts')", collapse = " + ")," + ", 
        paste0("s(",vars_endemic,", k = ", k_end,", bs = 'ts')+s(HelminthPrev, k=6 ,bs = 'ts', by=reg)", collapse = " + ")) %>% 
    as.formula()
}



#interaction model

#India models

mu_splineIndia_all <- function(response, offset, vars_numeric, vars_endemic, k_num = 4, k_end = 4){
  paste(response," | rate(",offset,") ~ 1 + ", 
        paste0("s(",vars_numeric,", k = ", k_num,", bs = 'ts')", collapse = " + ")," + ", 
        paste0("s(",vars_endemic,", k = ", k_end,", bs = 'tp')", collapse = " + ")) %>% 
    as.formula()
} 


mu_splineIndia_subModel <- function(response, offset, vars_numeric,  k_num = 6){
  paste(response," | rate(",offset,") ~ 1 + ", 
        paste0("s(",vars_numeric,", k = ", k_num,", bs = 'ts')", collapse = " + ")) %>% 
    as.formula()
} 


mu_linearIndia <- function(response, offset, vars_numeric, vars_endemic){
  paste(response," | rate(",offset,") ~ 1 + ",
        paste(vars_factor, collapse = " + "), " + ",
        paste0(vars_endemic, collapse = " + ")," + ", 
        paste0(vars_numeric, collapse = " + ")) %>%
    as.formula()
}

mu_splineIndia2 <- function(response, vars_num,  vars_end, k_num = 5, k_end = 5){
  paste(response," | rate(log_tests_pp) ~ 1  + ",
        paste0("s(",vars_num,", k = ", k_num,", bs = 'ts')", collapse = " + ")," + ",
        paste0("s(",vars_end,", k = ", k_end,", bs = 'tp')", collapse = " + ")) %>%
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
    #ylim(c(NA,y_max))
    ylim(c(-10,10))
  
}
