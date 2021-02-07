#---------------------------------------------------------------------------
# Script: Archaeological radiocarbon dates analysis
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------

# Load libraries for functions used
library(rcarbon) #calibrate 14C dates
library(tidyverse) #maipulate dataframes

#Read radiocarbon archaeological sites from Ziolkowski 1994 
arch_ecuador <- read.csv("data/archaeologicalsites/archaeological-andes-Ziolkowski.csv") %>%
  filter(!C14_date==0) %>%
  filter(!C14_date>=4000) %>%
  filter(!Error>=200) %>%
  filter(!Material=="unkwown")

#running calibration with Southern Hemisphere cc
arch.caldates <- calibrate(x=arch_ecuador$C14_date, errors=arch_ecuador$Error,calCurves='shcal13') 

#plot caibrated radiocarbon ages
plot(arch.caldates,HPD=TRUE)

#spd() aggregates (sums) calibrated radiocarbon dates within a defined chronological range
Ecuador.spd <- spd(arch.caldates,timeRange=c(3000,100)) 
plot(Ecuador.spd, title="Test") 

#extract and write radiocarbon calibrated ages
ecuadorspd_ts <- Ecuador.spd$grid
write.csv(ecuadorspd_ts, "data/archaeologicalsites/ecuadorspd_ts.csv")



   