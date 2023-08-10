#######################################################################
########## PART 1: Extracting excess deaths annd P score estimates
#######################################################################

library(tidyr)
library(readr)
library(dplyr)
library(brms)
library(stringr)
# https://www.nature.com/articles/

# P score = 100 * (acm- expected)/expected

#note that the raw data used here is too large for github. 
#Please Find here: https://github.com/WHOexcessc19/Codebase

acm_raw <- read_csv("data/acm.csv")
expected_raw <- read_csv("data/expected.csv")

acm <-  
  acm_raw %>% 
  select(-WHO_region) %>% 
  rename(country = Country,
         month = months) %>% 
  pivot_longer(starts_with("acm"), names_to = "sample", values_to = "acm") %>% 
  group_by(country,month) %>% 
  mutate(sample = str_replace(sample,"acm","s"))


expected <-  
  expected_raw %>% 
  select(-WHO_region) %>% 
  rename(country = Country,
         month = months) %>% 
  pivot_longer(starts_with("expec"), names_to = "sample", values_to = "expected") %>% 
  group_by(country,month) %>% 
  mutate(sample = str_replace(sample,"expec","s"))
  

data <- full_join(acm,expected,by = c("country","month","sample"))

data_final <- data %>% 
  mutate(pscore = 100 * (acm- expected)/expected) %>% 
  group_by(country,sample) %>% 
  summarise(pscore = mean(pscore)) %>% # average over months per sample
  group_by(country) %>% 
  summarise(pscore.mean = mean(pscore),  # average over posterior samples
            pscore.mean.se = sd(pscore)/sqrt(1000), # standard error of the posterior mean
            pscore.lo = quantile(pscore,0.025),
            pscore.hi = quantile(pscore,0.975))%>% 
  mutate(plo = pscore.mean - 2*pscore.mean.se,
         plo = pscore.mean + 2*pscore.mean.se) %>% 
  filter(country=="Azerbaijan")

write.csv(data_final, 'pscores.csv')





# random stuff for spatial autocorrelation models

east <- north <- 1:10
Grid <- expand.grid(east, north)
K <- nrow(Grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1
W %>% View


# generate some spatial data
east <- north <- 1:10
Grid <- expand.grid(east, north)
K <- nrow(Grid)

# set up distance and neighbourhood matrices
distance <- as.matrix(dist(Grid))
W <- array(0, c(K, K))
W[distance == 1] <- 1

# generate the covariates and response data
x1 <- rnorm(K)
x2 <- rnorm(K)
theta <- rnorm(K, sd = 0.05)
phi <- rmulti_normal(
  1, mu = rep(0, K), Sigma = 0.4 * exp(-0.1 * distance)
)
eta <- x1 + x2 + phi
prob <- exp(eta) / (1 + exp(eta))
size <- rep(50, K)
y <- rbinom(n = K, size = size, prob = prob)
dat <- data.frame(y, size, x1, x2)

# fit a CAR model
fit <- brm(y | trials(size) ~ x1 + x2 + car(W),
           data = dat, data2 = list(W = W),
           family = binomial(), cores =4)
summary(fit)

fit$model

#calculating neighbours

# install.packages(c("rgeos", "rgdal"), type="source")
library(rgdal)
library(rgeos)
library(spdep)


# Variables for getting map shapefiles
map.url <-  "http//https://www.naturalearthdata.com/downloads/50m-cultural-vectors/50m-admin-0-countries-2/.zip"

map.path <- file.path("Raw data", "maps")
map.zip.name <- basename(map.url)
map.name <- tools::file_path_sans_ext(map.zip.name)

# Download Natural Earth shapefiles
# download.file(url=map.url, file.path(map.path, map.zip.name), "auto")

library('sf')
countries <- read_sf(dsn='ne_50m_admin_0_countries.shp')

# Extract the ISO codes and map them to the numeric row names
country.names <- countries$NAME_LONG

# Determine which countries are neighbors
# Adapted from http://stackoverflow.com/a/32318128/120898
#
# spdep::poly2nb/nb2mat method is faster and more accurate than rgeos::gTouches
#
# gTouches gives wrong results; doesn't see Russia-North Korea border; is suuuuuuuper slow
#   neighbor.matrix <- gTouches(countries, byid=TRUE)
#
neighbor.list <- poly2nb(countries)
neighbor.matrix <- nb2mat(neighbor.list, style="B", zero.policy=TRUE)
colnames(neighbor.matrix) <- country.names
rownames(neighbor.matrix) <- country.names

write.csv(neighbor.matrix, 'neighbor.matrix.csv') #need country names fixed
