pacman::p_load("tidyverse","brms", 'missForest', 'sp', 'spdep', 'spatialreg')

rm(list=ls())

################################################
#----------------Import raw data----------------
################################################

TotalData<-read.csv('data/RawDataApr2021.csv', row.names=1)

#Remove countries with 0 reported cases or deaths
TotalDataNoZeros <- filter(TotalData, Cases._perMil  > 0, Deaths._perMil >0)

#-----------------------------------------------------------------------------------
#create a lag variable to account for neighbouring cases
#-----------------------------------------------------------------------------------

TotalDataNoZerosSpatial <- TotalDataNoZeros
coordinates(TotalDataNoZerosSpatial) <- ~latitude+longitude

nb <- tri2nb(TotalDataNoZerosSpatial ) #create neighbour list 
spatL <- nb2listw(nb)
plot(nb, coordinates(TotalDataNoZerosSpatial ))

TotalDataNoZeros$lag_rate<-lag.listw(x=spatL, var=(TotalDataNoZeros$Cases._perMil))

#-----------------------------------------------------------------------------------
#impute missing data
#-----------------------------------------------------------------------------------

catData <- Filter(is.character, TotalDataNoZeros) #all catergorical data complete. No need to impute

missForestData <- Filter(is.numeric, TotalDataNoZeros)

TotalDataNoNA <- missForest(TotalDataNoZeros [4:29], ntree = 1000, variablewise = TRUE)

# check error
TotalDataNoNA$OOBerror 

interpData <- cbind(data.frame(TotalDataNoZeros[1:3]), TotalDataNoNA$ximp)

#-----------------------------------------------------------------------------------
#Check for correlations
#-----------------------------------------------------------------------------------
cvd19data <- read.csv('data/cvdDataPreprocessed', row.names=1)

#check out correlations

ggpairs(cvd19data[4:27])

#this is a little hard to read. We focussed on numeric data

cvd19dataNoCat <- Filter(is.numeric, cvd19data)

ggcorr(cvd19dataNoCat[3:25], label = TRUE, label_size = 3,  hjust = 0.9, size = 3 ) #pearson correlations are default

write.csv(interpData , "data/cvdDataPreprocessed", row.names=TRUE)








