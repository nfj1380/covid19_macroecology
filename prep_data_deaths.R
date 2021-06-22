prepare_data_deaths <- function(data){
  data %>% 
    select( -latitude, -longitude, -Total_tests_per_thousand, -Transmission_Classification) %>% 
    rename(y = Deaths._perMil, reg = WHO_Region,
           rain = Avg_Rainfall91.2016_mm, temp = Avg_temp91.2016_c, 
           popDen = popDen_2018, 
           meanAge = mean._age, country = Country,
           Diabetes = Diabetes_prevalence, Hookworm = Hookworm.disease,
           invSimp = inv_simpsons, cardioDR = Cardiovasc_death_rate,
           cases = Cases._perMil,
           spatLag = Spatial_lag, IDBurd = inf_diseaseBurden) %>% 
    mutate(reg = reg %>% 
             str_replace(fixed("Eastern Mediterranean"), "EMed") %>% 
             str_replace(fixed("Western Pacific"), "WPac") %>% 
             str_replace(fixed("South-East Asia"), "SEAs") %>% 
             str_sub(1,4)) %>% 
    mutate(y = round(y)) %>% 
    mutate(across(where(is_character),factor))
  }