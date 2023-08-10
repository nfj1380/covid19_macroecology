##############################################################################################################
#----------------------------------------Part 3: Mapping Pscores-------------------------------------------------------
##############################################################################################################

# install.packages(c("cowplot", "googleway",  "ggrepel", 
# "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))
# 
# 
# install.packages(c( "ggspatial", "libwgeom", "rgeos",
#                     "rnaturalearth", "rnaturalearthdata"))
# 
library("sf")
library('tidyverse')
library("rnaturalearth")
library("rnaturalearthdata")
library('rgeos')
library('RColorBrewer')


# load and prepare data set-------------------------------
cvd19_pre <- read.csv("data/cvd19_2023_04_19.csv") %>% 
  rename(name_long = country)

cvd19_pre <- read.csv("data/RawDataOct2023.csv") %>% 
  rename(name_long = Country)


world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)



mapData <- world %>% arrange(name_long) %>% 
  mutate(name_long= str_replace(name_long, 'United Kingdom', 'The United Kingdom')) %>% 
  mutate(name_long= str_replace(name_long, 'United States', 'United States of America')) %>%
  mutate(name_long= str_replace(name_long, 'Tanzania', 'United Republic of Tanzania')) %>% 
  mutate(name_long= str_replace(name_long, 'Vietnam', 'Viet Nam')) %>%  
  mutate(name_long= str_replace(name_long, 'Venezuela', 'Venezuela (Bolivarian Republic of)')) %>% 
  mutate(name_long= str_replace(name_long, 'Iran', 'Iran (Islamic Republic of)')) %>%
  mutate(name_long= str_replace(name_long, 'Bolivia', 'Bolivia (Plurinational State of)')) %>%
  mutate(name_long= str_replace(name_long, 'Lao PDR', "Lao People's Democratic Republic"))

AllData_combined <- mapData %>% left_join(cvd19_pre, by='name_long')

str(AllData_combined )
str(world)


a <- ggplot(data= AllData_combined )+
  geom_sf(aes(fill= pscore.mean))+
  scale_fill_gradientn(colours=brewer.pal(7,"Reds"))+
  theme_bw()


#add schisto

library(gridExtra)
grid.arrange(a,b,c,d,e,f,g, h, i,j, nrow=3)

