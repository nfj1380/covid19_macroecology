prepare_data <- function(data){
  data %>% 
    select(-Deaths._perMil, - latitude, -longitude) %>% 
    rename(y = Cases._perMil, reg = WHO_Region, tests = Total_tests_per_thousand,
           rain = Avg_Rainfall91.2016_mm, temp = Avg_temp91.2016_c, 
           popDen = popDen_2018, trans = Transmission_Classification,
           meanAge = mean._age, country = Country,
           Diabetes = Diabetes_prevalence, Hookworm = Hookworm.disease,
           invSimp = inv_simpsons, cardioDR = Cardiovasc_death_rate,
           spatLag = Spatial_lag, IDBurd = inf_diseaseBurden) %>% 
    mutate(reg = reg %>% 
             str_replace(fixed("Eastern Mediterranean"), "EMed") %>% 
             str_replace(fixed("Western Pacific"), "WPac") %>% 
             str_replace(fixed("South-East Asia"), "SEAs") %>% 
             str_sub(1,4)) %>% 
    mutate(trans = trans %>% 
             str_replace(fixed("Clusters_of_cases"), "clust") %>% 
             str_replace(fixed("Community_transmission"), "comm") %>% 
             str_replace(fixed("Sporadic_cases"), "spor") %>% 
             str_replace(fixed("Pending"), "comm")) %>% 
    mutate(y = round(y)) %>% 
    mutate(across(where(is_character),factor))
}