prepare_data <- function(data){
  data %>% 
    select( -latitude, -longitude) %>% 
    rename(deaths = Deaths._perMil,
           cases =  Cases._perMil,
           reg = WHO_Region,
           rain = Avg_Rainfall91.2016_mm, 
           temp = Avg_temp91.2016_c, 
           popDen = popDen_2018, 
           meanAge = mean._age, 
           country = Country,
           diabetes = Diabetes_prevalence, 
           hookworm = Hookworm.disease,
           invSimp = inv_simpsons, 
           cardioDR = Cardiovasc_death_rate,
           cases = Cases._perMil,
           spatLag = Spatial_lag, 
           dBurd = inf_diseaseBurden,
           tests = Total_tests_per_thousand,
           trans = Transmission_Classification) %>% 
    mutate(reg = reg %>% 
             str_replace(fixed("Eastern Mediterranean"), "EMed") %>% 
             str_replace(fixed("Western Pacific"), "WPac") %>% 
             str_replace(fixed("South-East Asia"), "SEAs") %>% 
             str_sub(1,4)) %>% 
    mutate(y_d = deaths/1e6, 
           y_c = cases/1e6,
           PerUrb = PerUrb/100.1,
           logTests = log(tests*1000)) %>% 
    mutate(across(where(is_character),factor))
}
