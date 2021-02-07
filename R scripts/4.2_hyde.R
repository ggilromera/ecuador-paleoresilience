#---------------------------------------------------------------------------
# Script: Historical human and cropland time series analysis 
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------

# Load libraries for the functions used
library(dplyr)
library(tidyverse)
library(cowplot)

## Load HYDE time series
#for all information on the HYDE 3.2 data used, see  ftp://ftp.pbl.nl/hyde/hyde3.2/readme_release_HYDE3.2.1.txt
# Read in Historical Human population density
load("data/HYDE/sitesHydes-human-1700-2000ADLong.Rdata")

# Average over the lower, upper and baseline estimations
humanPop <- sitesHydeL %>% 
  group_by(site, time,scenario) %>% 
  summarise(pop_avg=mean(pop,na.rm=T)) %>% 
  spread(time, pop_avg)

humanPop17002000 <- humanPop
load("data/HYDE/sitesHydes-human-0-1700ADLong.Rdata")
humanPop0_1700 <- humanPop[,-ncol(humanPop)] #remove 1700 variable

humanPopdf <- merge(humanPop0_1700, humanPop17002000, by="site")
colnames(humanPopdf) <- paste(colnames(humanPopdf), "AD", sep = "_")
#write.csv(croplanddf, "data/HYDE/HumanDensity_baseline.csv")


# Read in Historical Cropland area
load("data/HYDE/sitesHydes-cropland-0-1700ADLong.Rdata")

# Average over the lower, upper and baseline estimations
cropland <- sitesHydeL %>% 
  group_by(site, time,scenario) %>% 
  summarise(cropland_avg=mean(cropland,na.rm=T)) %>% 
  spread(time, cropland_avg)

cropland0_1700 <- cropland[,-ncol(cropland)]

load("data/HYDE/sitesHydes-cropland-1700-2000ADLong.Rdata")
colnames(sitesHydeL)[3] <- c("cropland")
cropland1700_2000 <- cropland #remove 1700 variable

croplanddf <- merge(cropland0_1700, cropland1700_2000, by="site")
colnames(croplanddf) <- paste(colnames(croplanddf), "AD", sep = "_")
#write.csv(croplanddf, "data/HYDE/Cropland_baseline.csv")


### HYDE human population density
humanPopdf <- read.csv("data/HYDE/HumanDensity_baseline.csv", row.names = 1)
colnames(humanPopdf) <- gsub("X", "", colnames(humanPopdf))
colnames(humanPopdf) <- gsub("_AD", "", colnames(humanPopdf))

select.lakes <- paste(c("Piñan1", "Yahurcch", "Fondocch1", "Llaviucu"), collapse = '|')
test <- humanPopdf %>%
  filter(str_detect(sites, select.lakes))
humanPopLakes <- as.data.frame(t(test))[-1,]
colnames(humanPopLakes) <- c("Yahuarcocha", "Fondococha", "Llaviucu", "Pinan")
humanPopLakes$variable <- "HumanDensity"


## HYDE cropland
cropland <- read.csv("data/HYDE/Cropland_baseline.csv", row.names=1)
colnames(cropland) <- gsub("X", "", colnames(cropland))
test <- cropland %>%
  filter(str_detect(sites, select.lakes))
croplandLakes <- as.data.frame(t(test))[-1,]
colnames(croplandLakes) <- c("Yahuarcocha", "Fondococha", "Llaviucu", "Pinan")
croplandLakes$variable <- "cropland"

#Rbind historical human and cropland data
humanfootprint <- rbind(humanPopLakes, croplandLakes)
indx <- sapply(humanfootprint, is.factor)
humanfootprint[indx] <- lapply(humanfootprint[indx], function(x) as.numeric(as.character(x)))
humanList <- split(humanfootprint, humanfootprint$variable)

years <- as.numeric(rownames(humanPopLakes))
humanfootprint <- cbind(humanfootprint, years)

# plot non-binned humanfootprint data
long_matrixHumanHist<-gather(humanfootprint, key=lake, value=value, -variable, -years)
head(long_matrixHumanHist)

plt_human_ts <- ggplot(long_matrixHumanHist, aes(x = years, y = value, colour=lake)) +
  geom_line()+
  #scale_y_continuous(trans='log10') +
  facet_wrap(~ variable, ncol = 1, scales = "free") +
  theme_bw()
plt_human_ts

# Write Human and cropland historical time series
#write.csv(long_matrixHumanHist, "data/HYDE/long_matrixHumanHist.csv")


## bin Human predictors
binhumanData <- function(i, cores, ...) {
  core <- cores[[i]]
  years <- rownames(humanList$HumanDensity) #this is to get ages_AD
  core <- core[ , -which(colnames(core) %in% c("variable"))] 
  
  # Function for binning samples in each core (Alistair Seddon)
  source("R scripts/functions/functions.R")
  
  #make the age categories
  diatomBin1 <- binFunc(as.data.frame(core), as.numeric(years), 100, 0, 1700) ##
  diatomBin2 <- binFunc(as.data.frame(core), as.numeric(years), 30, 1700, 2000) ##
  
  #diatomBin = binFunc(diatoms, as.numeric(rownames(diatoms)), 40, 20, 2660)
  
  #merge the two binning dataframes
  rbind.data.frame(diatomBin1, diatomBin2)
}

#wrap up the function
binnedHuman <- lapply(seq_along(humanList), binhumanData, cores = humanList)
names(binnedHuman) <- names(humanList)

#save the list of historical footprint predictors
saveRDS(binnedHuman, "data/binnedHuman.rds")



