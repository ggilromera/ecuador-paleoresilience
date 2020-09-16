## COMMUNITY SYNCHRONY ANALYSIS
setwd("/Volumes/xbenitogranell-data/0_project/data/training set")
# setwd("/nfs/xbenitogranell-data/0_project/data/training set")

#load packages
library(synchrony)
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(TeachingDemos)
library(maps)
library(analogue)
library(ggplot2)
library(plotrix)

# community synchrony function from Pedersen et al 2017 
community.sync.window <- function (d, twin=NULL, nrands=50) {
  years = as.numeric(rownames(d))
  starts=seq(from=1, to=length(years)-twin+1, by=1)
  ends=starts+twin-1
  mid=starts+(twin-1)/2
  
  results=matrix(nrow=length(starts), ncol=6, NA)
  for (i in 1:length(starts)) {
    locs=which(rownames(d) %in% years[starts[i]:ends[i]])
    c = community.sync(d[locs,], nrands=nrands, method = "spearman")
    results[i,] = c(years[starts[i]], years[ends[i]],years[mid[i]], c$meancorr, c$obs, c$pval)
  }
  colnames(results)=c("start", "end","mid", "meancorr", "sync", "pval")
  results=as.data.frame(results)
  return(results)
}


# read diatom Ecuadorean core datasets 30-yr binned 
diat_data <- read.csv("diat_binned_data.csv", row.names=1)
diat_data[is.na(diat_data)] <- 0

  llav <- diat_data %>% filter(lake=="Llaviucu")
  rownames(llav) <- llav$age
  llav <- llav[,!names(llav) %in% c("age", "lake")]
  llav <- llav[, colSums(llav) > 0]
  llav <- decostand(llav, method = "hellinger")


yah <- diat_data %>% filter(lake=="Yahuarcocha")
pin <- diat_data %>% filter(lake=="Pinan")
fondo <- diat_data %>% filter(lake=="Fondococha")

#comb1 <- analogue::join(llav, yah)

comb1 <- merge(pin, fondo, by="age")
ages <- comb1$age 
comb1 <- comb1[,3:ncol(comb1)]
comb1 <- comb1[,!names(comb1) %in% c("lake.y")]
comb1 <- comb1[, colSums(comb1) > 0]
rownames(comb1) <- ages
comb1 <- decostand(comb1, method = "hellinger")

#twindow = 10
sync_allsp = community.sync.window(llav, twin=30)

plot.ts(sync_allsp$sync)





#read diatom core datasets
mergedCores <- read.csv("mergedCores_counts3.csv") #read dataframe with diatom absolute counts including Fondococha

agedepth <- mergedCores[, names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- mergedCores[, !names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat[is.na(diat)] <- 0
diatoms_save <- cbind(agedepth, diat)

coresList <- split(diatoms_save, diatoms_save$lake)


# this is function to calculate relative aundance from counts data
RA <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[ , -which(names(core) %in% c("depth","upper_age", "lower_age", "lake"))] # drop year & depths vars
  core <- tran(core, method="percent")
  core[is.na(core)] <- 0
  depth <- coresList[[i]]$depth
  upper_age <- coresList[[i]]$upper_age
  lower_age <- coresList[[i]]$lower_age
  lake <- coresList[[i]]$lake
  cbind.data.frame(depth,upper_age, lower_age, lake, core) #combine extracted columns and remove first row to match with scd
  
}

## apply RA function function to each core
coresRA <- lapply(seq_along(coresList), RA, cores=coresList)
names(coresRA) <- names(coresList)

#extract dataframes from list
merged <- plyr::ldply(coresRA, data.frame)
merged <- merged[,-1] #remove .id variable


#Select most common species
agedepth <- merged[, names(merged) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- merged[, !names(merged) %in% c("depth", "upper_age", "lower_age", "lake")]

abund <- apply(diat, 2, max)
n.occur <- apply(diat>0, 2, sum)
diat_red <- diat[, n.occur >10 & abund>5] #more than 3% of RA and present in >2 samples

#check rows with NA
row.has.na <- apply(diat_red, 1, function(x){any(is.na(x))})
sum(row.has.na)

diatoms_save <- cbind(agedepth, diat_red)

#new1: ecological groups
#new2: harmonized taxonomic names

#transform dataframe to tidy format
#mergedCores <- read.csv("mergedCores_counts3.csv") #read dataframe with diatom absolute counts including Fondococha
mergedCores <- read.csv("mergedCores_counts4.csv")[-1] #with new Fondococha agedepth model

changes <- read.csv("old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)

diatoms_save <- mergedCores #save dataframe

new <- diatoms_save %>%
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  group_by(depth, taxa, lake, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == 0) %>% #this is to remove empty samples (rows)
  filter(!upper_age == 0) %>% #this is to remove ages == 0 (triumfo and fondodocha record)
  spread(key = taxa, value = count) %>%
  as.data.frame()

cores_merged <- new %>%
  mutate(AgeCE = upper_age*(-1)+1950) %>%
  filter(AgeCE >= 0)

cores_merged <- split(cores_merged, cores_merged$lake)


# Filter which lake I want, remove spp=0 and set ages
lake <- cores_merged %>%
  filter(lake == "yahuarcocha")

#agedepth <- lake[, 1:4]
agedepth <- lake[, names(lake) %in% c("depth", "lake", "upper_age", "lower_age", "AgeCE")]

diatoms <- lake[,5:ncol(lake)]
diatoms[is.na(diatoms)] <- 0
diatoms <- diatoms[, colSums(diatoms) > 0]
#rownames(diatoms) <- agedepth$upper_age
rownames(diatoms) <- agedepth$AgeCE
diatoms <- diatoms[,!names(diatoms) %in% c("AgeCE")]
diatoms <- decostand(diatoms, method = "hellinger")

#twindow = 10
sync_allsp = community.sync.window(diatoms, twin=20)
plot.ts(sync_allsp$sync)

## synchrony by ecological groups
new <- diatoms_save %>% 
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_1)) %>% #ecological grouping
  group_by(depth, taxa, lake, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == "0" ) %>% #this is to remove empty samples (rows)
  filter(!upper_age==0.0) %>% #this is to drop extraneous ages
  ungroup() %>%
  group_by(depth, lake) %>%
  mutate(relative_abundance_percent = count / sum(count) * 100) %>%
  dplyr::select(depth, lake, upper_age, taxa, relative_abundance_percent) %>%
  spread(key = taxa, value = relative_abundance_percent) %>%
  ungroup()

coresecologicalgroups <- split(new, new$lake)

## binndata
doBin <- function(i, cores, ...) {
  # Function for binning samples in each core (Alistair Seddon)
  source("/Volumes/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
  #source("/nfs/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
  
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "lake", "upper_age", "lower_age"))] # drop year & depths vars
  core[is.na(core)] <- 0
  core <- core[, colSums(core) > 0] #select only present species
  ages <- coresecologicalgroups[[i]]$upper_age #this is to get ages

  #make the age categories
  diatomBin1 <- binFunc(as.data.frame(core), as.numeric(ages), 30, -60, 2000) 
  #diatomBin2 <- binFunc(as.data.frame(core), as.numeric(ages), 100, 150, 2000) 

  #merge the two binning dataframes for North lakes
  rbind.data.frame(diatomBin1)
  
}

coresecologicalgroups_binned <- lapply(seq_along(coresecologicalgroups), doBin, cores = coresecologicalgroups) 
names(coresecologicalgroups_binned) <- names(coresecologicalgroups)

ages <- plyr::ldply(as.numeric(rownames(coresecologicalgroups_binned$llaviucu, data.frame)))
matrix <- plyr::ldply(coresecologicalgroups_binned, data.frame)

matrix <- plyr::rename(matrix,c(".id"="lake"))
diat <- cbind(ages,matrix)
colnames(diat)[1] <- "age"

llav <- diat %>% filter(lake=="llaviucu")
yah <- diat %>% filter(lake=="yahuarcocha")
pinan <- diat %>% filter(lake=="pinan")
fondo <- diat %>% filter(lake=="fondococha")

#merge dataframes
df <- plyr::join_all(list(llav, yah), by = 'age')
rownames(df) <- df$age
df <- df[,-which(names(df) %in% c("age", "lake"))]

head(df)
row.all.na <- apply(df, 1, function(x) all(is.na(x)))
sum(row.all.na)
df <- df[ !row.all.na, ]

#interpolate to fill empty bins 
library(zoo)
df_planktic <- df %>% dplyr::select(freshwater_planktic, freshwater_planktic.1,
                                  tycoplanktonic, tycoplanktonic.1)
df_benthic <- df %>% dplyr::select(benthic, benthic.1, ephiphytic, ephiphytic.1, saline, saline.1)

dfInter_planktic <- as.data.frame(na.approx(df_planktic, na.rm = TRUE)) #do interpolation between adjacent samples
dfInter_benthic <- as.data.frame(na.approx(df_benthic, na.rm = TRUE)) #do interpolation between adjacent samples

rownames(dfInter_planktic) <- rownames(df_planktic)
rownames(dfInter_benthic) <- rownames(df_benthic)

# perform synchcrony
twindow = 3 #30 years
sync_planktic<- community.sync.window(dfInter_planktic, twin=twindow)
sync_benthic<- community.sync.window(dfInter_benthic, twin=twindow)

#sync_allsp <- community.sync.window(dfInter[as.numeric(rownames(dfInter))<500,], twin=twindow)

# community.sync(dfInter[as.numeric(rownames(dfInter))<250,], nrands = 20)
# community.sync(dfInter[as.numeric(rownames(dfInter))>250,], nrands = 20)

# plot synchrony plots
png(filename ="sync_Llav_Yah.png",
     res = 300,                                             
     width = 10, height = 8, units = 'in')

alpha=0.05
plot(sync_planktic$end, sync_planktic$sync, type="l", lwd=2, lty=1, ylim=c(0,1),
     xlab="Cal yr BP", ylab="Community synchrony")
lines(sync_benthic$end, sync_benthic$sync, type="l", lwd=1, lty=2)

for(i in 1:4) abline(v=zones_llaviucu[i], col=i+1, lty=2, lwd=2)
for(i in 1:4) abline(v=zones_yahuarcocha[i], col=i+1, lty=2, lwd=2)
# for(i in 1:4) abline(v=zones_pinan[i], col=i+1, lty=2, lwd=2)
# for(i in 1:4) abline(v=zones_fondococha[i], col=i+1, lty=2, lwd=2)

points(sync_planktic$end, sync_planktic$sync, col="black", pch=21, 
       bg=ifelse(sync_planktic[,"pval"] < alpha, 'black', 'white'))

points(sync_benthic$end, sync_benthic$sync, col="black", pch=21, 
       bg=ifelse(sync_benthic[,"pval"] < alpha, 'black', 'white'))

legend("topright",legend = c("Planktic", "Benthic"),
       lwd=c(2,1), lty=c(1,2), col="black",
       box.col="#FF003300",inset=0.005,cex=0.9)

dev.off()


#BC dissimilarity time series

## BC dissimilarity from consecutive and cumulative samples
#successive samples
## non-binned data
ts_BCdismilarity_north_lakes <- read.csv("ts_BCdismilarity_north_lakes.csv", row.names = 1)
tsref_BCdismilarity_north_lakes <- read.csv("tsref_BCdismilarity_north_lakes.csv", row.names = 1)


llav <- ts_BCdismilarity_north_lakes %>% filter(id=="llaviucu")
yah <- ts_BCdismilarity_north_lakes %>% filter(id=="yahuarcocha")
pin <- ts_BCdismilarity_north_lakes %>% filter(id=="pinan")
fondo <- ts_BCdismilarity_north_lakes %>% filter(id=="fondococha")

llav_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="llaviucu")
yah_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="yahuarcocha")
pin_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="pinan")
fondo_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="fondococha")

#read in fossil distance to modern analogues
dist_to_analogues_all_lakes<- read.csv("dist_to_analogues_all_lakes.csv",row.names = 1) 

axis.V<-1.1
label.V<-0.9
alpha<-0.05
seq_palette <- viridis(4)


png(filename ="BC_timeseries.png",
    res = 300,                                             
    width = 6, height = 8, units = 'in')

#par(mfrow=c(5,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)
par(mfrow=c(8,1),mar=c(2,5.5,0,2),oma=c(3,1,1,1),las=1)

plot(pin$age, pin$timeseriesBC, type="n", ylab="BC\ndissimilarity", xlab=NA,xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), xlim=c(-60,2000))
abline(v=c(764:783), col="grey", lty=1, lwd=3)
abline(v=c(946:1024), col="grey", lty=1, lwd=3)
lines(pin$age, pin$timeseriesBC, type="l", lwd=2, ylab="BC\ndissimilarity", xlab=NA,xaxt='n',
    cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), xlim=c(-60,2000), col=seq_palette[1])
lines(pin_ref$age, pin_ref$timeseriesBC_fromreference, lwd=1, col=seq_palette[1])
for(i in 1:4) abline(v=zones_pinan[i], col=i+1, lty=2, lwd=2)
dist_to_analogues_pinan <- dist_to_analogues_all_lakes %>% 
    filter(lake=="Piñan") %>%
    filter(quality=="bad")
points(dist_to_analogues_pinan$ages,rep(max(dist_to_analogues_pinan$dist_to_analogues)*1.03,
                     length(dist_to_analogues_pinan$ages)),pch=25,bg="darkorange",col="black",cex=1.5)
legend("topright",legend = c("Piñan"),inset=0.010,cex=0.9, bg="transparent", bty="n")

plot(yah$age, yah$timeseriesBC, type="l", lwd=2, ylab="BC\ndissimilarity", xlab=NA,
     xaxt='n',cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), xlim=c(-60,2000),col=seq_palette[2])
lines(yah_ref$age, yah_ref$timeseriesBC_fromreference, lwd=1, col=seq_palette[2])
for(i in 1:4) abline(v=zones_yahuarcocha[i], col=i+1, lty=2, lwd=2)
dist_to_analogues_yah <- dist_to_analogues_all_lakes %>% 
  filter(lake=="Yahuarcocha") %>%
  filter(quality=="bad")
points(dist_to_analogues_yah$ages,rep(max(dist_to_analogues_yah$dist_to_analogues)*1.10,
                                        length(dist_to_analogues_yah$ages)),pch=25,bg="darkorange", col="black", cex=1.5)
legend("topright",legend = c("Yahuarcocha"), bg="transparent", bty="n",inset=0.010)

plot(fondo$age[1:79], fondo$timeseriesBC[1:79], type="n", ylab="BC\ndissimilarity", xlab=NA,xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), xlim=c(-60,2000))
abline(v=c(1996:2063), col="grey", lty=1, lwd=3)
lines(fondo$age[1:79], fondo$timeseriesBC[1:79], type="l", lwd=2, ylab="BC\ndissimilarity",
     xlab=NA, xaxt='n', cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1),xlim=c(-60,2000),
     col=seq_palette[3])
lines(fondo_ref$age[1:79], fondo_ref$timeseriesBC_fromreference[1:79], lwd=1, col=seq_palette[3])
for(i in 1:5) abline(v=zones_fondococha[i], col=i+1, lty=2, lwd=2)
dist_to_analogues_fondo <- dist_to_analogues_all_lakes %>% 
  filter(lake=="Fondococha") %>%
  filter(quality=="bad")
points(dist_to_analogues_fondo$ages,rep(max(dist_to_analogues_fondo$dist_to_analogues)*0.96,
                                      length(dist_to_analogues_fondo$ages)),pch=25,bg="darkorange", col="black", cex=1.5)
legend("topright",legend = c("Fondococha"), bg="transparent", bty="n",inset=0.010)

plot(llav$age[1:68], llav$timeseriesBC[1:68], type="l", lwd=2, ylab = "BC\ndissimilarity",
     xlab=NA,xaxt='n', cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), xlim=c(-60,1900),col=seq_palette[4])
lines(llav_ref$age[1:68], llav_ref$timeseriesBC_fromreference[1:68], type="l", lwd=1, ylab = "BC dissimilarity from reference", xlab=NA, xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), col=seq_palette[4])
for(i in 1:4) abline(v=zones_llaviucu[i], col=i+1, lty=2, lwd=2)
dist_to_analogues_llaviucu <- dist_to_analogues_all_lakes %>% 
  filter(lake=="Llaviucu") %>%
  filter(quality=="bad")
points(dist_to_analogues_llaviucu$ages,rep(max(dist_to_analogues_llaviucu$dist_to_analogues)*1.06,
                                      length(dist_to_analogues_llaviucu$ages)),pch=25,bg="darkorange", col="black", cex=1.5)
legend("topright",legend = c("Llaviucu"), bg="transparent", bty="n",inset=0.010,cex=0.9)


# Here add human and climatic historical timeseries
library(zoo)
long_matrixHumanHist <- read.csv("long_matrixHumanHist.csv", row.names=1) %>%
  mutate(age=1950-years)
var <- "cropland"
#var <- "HumanDensity"

llav_human <- long_matrixHumanHist %>% 
  filter(lake=="Llaviucu", variable==var)
llav_human$value[36:38] <- NA
llav_human$value[31:41] <- na.approx(llav_human$value[31:41], na.rm = TRUE) #do interpolation between adjacent samples

yah_human <- long_matrixHumanHist %>% filter(lake=="Yahuarcocha", variable==var)
pin_human <- long_matrixHumanHist %>% filter(lake=="Pinan", variable==var)
fondo_human <- long_matrixHumanHist %>% filter(lake=="Fondococha", variable==var)

#plot
plot(yah_human$age, yah_human$value, type="l", lwd=2, ylab=bquote(atop("Cropland area",
                                                                       ~ Km^2)), 
     xlab=NA,xaxt='n', cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col=seq_palette[2])
lines(pin_human$age, pin_human$value, type="l", lwd=2, col=seq_palette[1])
lines(fondo_human$age, fondo_human$value, type="l", lwd=2, col=seq_palette[3])
lines(llav_human$age, llav_human$value, type="l", lwd=2, col=seq_palette[4])
legend("topright",legend = c("Piñan", "Yahuarcocha", "Fondococha", "Llaviucu"),
       lwd=2,col=c(seq_palette[1], seq_palette[2],seq_palette[3],seq_palette[4]),
       inset=0.005,cex=0.9,bg="transparent", bty="n")

    #read in simulated climatic data from Paleoview
    binnedClimateSimulated <- readRDS("/Volumes/xbenitogranell-data/0_project/data/historic/CommonEra/binnedClimateSimulated.rds")
    matrixClimateHist <- plyr::ldply(binnedClimateSimulated, data.frame)
    
    #this is to extract binned ages; only needs to do it once because binned ages are the same for all the cores
    ages <- as.data.frame(as.numeric(row.names(binnedClimateSimulated$PrecipSeason)))
    ages <- plyr::rename(ages,c("as.numeric(row.names(binnedClimateSimulated$PrecipSeason))"="age"))
    
    matrixClimHist <- cbind(matrixClimateHist, ages)
    
    long_matrixClimHist<-gather(matrixClimHist, key=lake, value=value, -age, -.id)

#read pages2k temperature data
long_SAtempdata <- read.csv("/Volumes/xbenitogranell-data/0_project/data/climate/long_SAtempdata.csv", row.names = 1) %>%
  filter(method=="LNA_anomalies")

plot(long_SAtempdata$age_calBP, long_SAtempdata$Temp, type="l", lwd=2, ylab="Temperature\n\u00B0C",
     cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col="black", xlab=NA,xaxt='n')
legend("topright",legend = c("PAGES2K Temperature anomalies reconstruction"), 
       bg="transparent", bty="n",inset=0.010,cex=0.9)

### read South America hydroclimate records
pumacocha <- read.csv("/Volumes/xbenitogranell-data/0_project/data/climate/SouthAmerica_hydroclimate_records.csv", row.names = 1) %>%
  filter(core %in% c("pumacocha")) %>% filter(!d13C==-9999.00)

# triumfo <- read.csv("/Volumes/xbenitogranell-data/0_project/data/climate/SouthAmerica_hydroclimate_records.csv", row.names = 1) %>%
#   filter(core %in% c("triumfo")) %>% filter(!d13C==0.0)
# ubaque <- read.csv("/Volumes/xbenitogranell-data/0_project/data/climate/SouthAmerica_hydroclimate_records.csv", row.names = 1) %>%
#   filter(core %in% c("ubaque")) %>% filter(!d13C==-9999.00)

plot(pumacocha$age_calBP, pumacocha$d18O, type="l", lwd=2, ylab=expression(paste (delta^18, "O \u2030")),
     cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col="black", xlab=NA,xaxt='n')
legend("topright",legend = c("Pumacocha lake"), bg="transparent", bty="n",inset=0.010,cex=0.9)

#read archaecological time series
ecuadorspds <- read.csv("/Volumes/xbenitogranell-data/0_project/data/archaeology/ecuadorspd_ts.csv", row.names = 1)

plot(ecuadorspds$calBP, ecuadorspds$PrDens, type="h", lwd=2, ylab="SPD",
     cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col="grey")
legend("topright", legend=c("Ecuador radiocarbon dates"), bg="transparent", bty="n",inset=0.010,cex=0.9)

mtext("Cal yr BP",side=1,outer=T,line=1,  adj=0.525,cex=1)
dev.off()


##

plot(altN$age, altN$tau, type="l", lwd=2, ylab = "Spearman rho",
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1), col=seq_palette[1])
abline(h=0, lwd=1,lty=2, col=8)
points(altN$age, altN$tau, pch=21, bg=ifelse(altN$pval < alpha, seq_palette[1], 'white'))
lines(altS$age, altS$tau, lwd=2, col=seq_palette[2])
points(altS$age, altS$tau, pch=21, bg=ifelse(altS$pval < alpha, seq_palette[2], 'white'))
lines(latS$age, latS$tau, lwd=2, col=seq_palette[3])
points(latS$age, latS$tau, pch=21, bg=ifelse(latS$pval < alpha, seq_palette[3], 'white'))
lines(latN$age, latN$tau, lwd=2, col=seq_palette[4])
points(latN$age, latN$tau, pch=21, bg=ifelse(latN$pval < alpha, seq_palette[4], 'white'))
lines(pair5$age, pair5$tau, lwd=2, col=seq_palette[5])
points(pair5$age, pair5$tau, pch=21, bg=ifelse(pair5$pval < alpha, seq_palette[5], 'white'))
lines(pair6$age, pair6$tau, lwd=2, col=seq_palette[6])
points(pair6$age, pair6$tau, pch=21, bg=ifelse(pair6$pval < alpha, seq_palette[6], 'white'))

legend("topright",legend = c("Llaviucu-Yahuarcocha", "Llaviucu-Fondococha", 
                            "Llaviucu-Pinan", "Yahuarcocha-Fondococha", "Yahuarcocha-Pinan",
                            "Fondococha-Pinan"),
       lwd=2,col=c(seq_palette[1],seq_palette[2], seq_palette[3], seq_palette[4], seq_palette[5], seq_palette[6]),
       box.col="#FF003300",inset=0.005,cex=0.7)

mtext("Cal yr BP",side=1,outer=T,line=1,  adj=0.525,cex=1)

###
#Altitude
par(mfrow=c(3,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)

plot(alt$age[1:57], alt$tau[1:57], type="l", lwd=2, ylab = "Mann-Kendall Tau", xlab=NA, xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1), col=seq_palette[1])
abline(h=0, lwd=1,lty=2, col=8)
text(-0.4,-0.9,paste("Community synchrony across elevation"),pos=4)
lines(alt$age[58:nrow(alt)], alt$tau[58:nrow(alt)], lwd=2, col=seq_palette[2])
points(alt$age[1:57], alt$tau[1:57], pch=21, 
       bg=ifelse(alt$pval < alpha, seq_palette[1], 'white'))
points(alt$age[58:nrow(alt)], alt$tau[58:nrow(alt)], pch=21, 
       bg=ifelse(alt$pval < alpha, seq_palette[2], 'white'))
legend("topleft",legend = c("Llaviucu-Yahuarcocha", "Pinan-Fondococha"),
       lwd=2,col=c(seq_palette[1], seq_palette[2]),box.col="#FF003300",inset=0.005,cex=0.9)

#Latitude
plot(lat$age[1:55], lat$tau[1:55], type="l", lwd=2, ylab = "Mann-Kendall Tau", xlab=NA, xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1), col=seq_palette[3])
abline(h=0, lwd=1,lty=2, col=8)
text(-0.4,-0.9,paste("Community sychrony across latitude"),pos=4)
lines(lat$age[56:nrow(lat)], alt$tau[56:nrow(lat)], lwd=2, col=seq_palette[4])
points(lat$age[1:55], lat$tau[1:55], pch=21, 
       bg=ifelse(lat$pval < alpha, seq_palette[3], 'white'))
points(lat$age[56:nrow(lat)], alt$tau[56:nrow(lat)], pch=21, 
       bg=ifelse(lat$pval < alpha, seq_palette[4], 'white'))
legend("topleft",legend = c("Pinan-Yahuarcocha", "Llaviucu-Fondococha"),
       lwd=2,col=c(seq_palette[3], seq_palette[4]),box.col="#FF003300",inset=0.005,cex=0.9)


# Compare across altitude and latitude
plot(corr_alt$age, corr_alt$mean, type="l", lwd=2, ylab="Mann-Kendall Tau", xlab="cal yr BP",
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1))
abline(h=0, lwd=1,lty=2, col=8)
plotCI(corr_alt$age,corr_alt$mean,uiw=corr_alt$se,add=T,pch=NA,sfrac=0)
lines(corr_lat$age, corr_lat$mean, col="blue",lwd=2)
plotCI(corr_lat$age,corr_lat$mean,uiw=corr_lat$se,add=T,pch=NA,sfrac=0, col="blue")
legend("topright",legend = c("Across elevation", "Across latitude"),
       lwd=2,col=c("black", "blue"),box.col="#FF003300",inset=0.005,cex=0.9)
mtext("Cal yr BP",side=1,outer=T,line=1,  adj=0.525,cex=1)

dev.off()


## BC dissimilarity to baseline
#Altitude
par(mfrow=c(3,1),mar=c(2,5,0,2),oma=c(3,1,1,1),las=1)

plot(alt$age[1:58], alt$tau[1:58], type="l", lwd=2, ylab = "Mann-Kendall Tau", xlab=NA, xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1), col=seq_palette[1])
abline(h=0, lwd=1,lty=2, col=8)
text(-0.4,-0.9,paste("Community synchrony across elevation to baseline"),pos=4)
lines(alt$age[59:nrow(alt)], alt$tau[59:nrow(alt)], lwd=2, col=seq_palette[2])
points(alt$age[1:58], alt$tau[1:58], pch=21, 
       bg=ifelse(alt$pval < alpha, seq_palette[1], 'white'))
points(alt$age[59:nrow(alt)], alt$tau[59:nrow(alt)], pch=21, 
       bg=ifelse(alt$pval < alpha, seq_palette[2], 'white'))
legend("topleft",legend = c("Llaviucu-Yahuarcocha", "Pinan-Fondococha"),
       lwd=2,col=c(seq_palette[1], seq_palette[2]),box.col="#FF003300",inset=0.005,cex=0.9)

#Latitude
plot(lat$age[1:56], lat$tau[1:56], type="l", lwd=2, ylab = "Mann-Kendall Tau", xlab=NA, xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1), col=seq_palette[3])
abline(h=0, lwd=1,lty=2, col=8)
text(-0.4,-0.9,paste("Community synchrony across latitude to baseline"),pos=4)
lines(lat$age[57:nrow(lat)], alt$tau[57:nrow(lat)], lwd=2, col=seq_palette[4])
points(lat$age[1:56], lat$tau[1:56], pch=21, 
       bg=ifelse(lat$pval < alpha, seq_palette[3], 'white'))
points(lat$age[57:nrow(lat)], alt$tau[57:nrow(lat)], pch=21, 
       bg=ifelse(lat$pval < alpha, seq_palette[4], 'white'))
legend("topleft",legend = c("Pinan-Yahuarcocha", "Llaviucu-Fondococha"),
       lwd=2,col=c(seq_palette[3], seq_palette[4]),box.col="#FF003300",inset=0.005,cex=0.9)

# Compare across altitude and latitude
plot(corr_alt$age, corr_alt$mean, type="l", lwd=2, ylab="Mann-Kendall Tau", xlab="cal yr BP",
     cex.lab=label.V,cex.axis=axis.V, ylim = c(-1,1))
abline(h=0, lwd=1,lty=2, col=8)
plotCI(corr_alt$age,corr_alt$mean,uiw=corr_alt$se,add=T,pch=NA,sfrac=0)
lines(corr_lat$age, corr_lat$mean, col="blue",lwd=2)
plotCI(corr_lat$age,corr_lat$mean,uiw=corr_lat$se,add=T,pch=NA,sfrac=0, col="blue")
legend("topright",legend = c("Across elevation", "Across latitude"),
       lwd=2,col=c("black", "blue"),box.col="#FF003300",inset=0.005,cex=0.9)
mtext("Cal yr BP",side=1,outer=T,line=1,  adj=0.525,cex=1)




##################################################################
## Calculates the Kendall Tau correlation for moving window
###################################################################

Tau<-function(ts.matrix,ages,window.size,step.size,type){
  library(Kendall)
  bin.start<-c(ages[length(ages)],(ages[length(ages)]-window.size))
  max.bins<-(ceiling(ages[length(ages)]-ages[1])/step.size)-2
  tau<-c()
  pval<-c()
  cor<-c()
  cor.p<-c()
  cor.p.ad<-c()
  
  midpoints<-c()
  nsamp<-c()
  if (type=="time") {
    bins<-matrix(NA,nrow=max.bins,ncol=2)
    for (i in 1:max.bins){
      bin.i<-bin.start-(i*step.size)
      bin.rows<-which(ages<=bin.i[1] & ages>=bin.i[2])
      ts.win<-ts.matrix[bin.rows,]
      nsamp[i]<-nrow(ts.win)
      tau[i]<-MannKendall(ts.win)$tau
      pval[i]<-MannKendall(ts.win)$sl
      cor[i]<-cor.test(ts.win[,1],ts.win[,2],, method = "spearman")$estimate
      cor.p[i]<-cor.test(ts.win[,1],ts.win[,2],, method = "spearman")$p.value
      cor.p.ad[i] <- p.adjust(cor.p[i], "holm")
      bins[i,]<-bin.i
    }	 
    midpoints<-bins[,2]+((bins[,1]-bins[,2])/2)
    out<-list(Taus=tau,Cor=cor,bins=bins,midpoints=midpoints,nsamp=nsamp, pval=pval, cor.p=cor.p,cor.p.ad=cor.p.ad)
  } else if (type=="sample") {
    bins<-matrix(NA,nrow=(nrow(ts.matrix)-window.size),ncol=2)
    for (i in (1:(nrow(ts.matrix)-window.size))){
      ts.win<-ts.matrix[c(i:(i+(window.size-step.size))),]
      nsamp[i]<-nrow(ts.win)
      bins[i,]<-ages[c(i,(i+(window.size-step.size)))]
      tau[i]<-MannKendall(ts.win)$tau
      pval[i]<-MannKendall(ts.win)$sl
      cor[i]<-cor.test(ts.win[,1],ts.win[,2], method = "spearman")$estimate
      cor.p[i]<-cor.test(ts.win[,1],ts.win[,2], method = "spearman")$p.value
      cor.p.ad[i] <- p.adjust(cor.p[i], "holm")
      
    }	
  } else {
    print("invalid type")
  }
  midpoints<-bins[,2]+((bins[,1]-bins[,2])/2)
  out<-list(Taus=tau,Cor=cor,bins=bins,midpoints=midpoints,nsamp=nsamp,pval=pval,cor.p=cor.p,cor.p.ad=cor.p.ad)
  return(out)
}


##read BC timeseries 30-year binned datasets
ts_BCdismilarity_north_lakes_bins <- read.csv("ts_BCdismilarity_north_lakes_bins.csv", row.names = 1)
ts_BCdismilarity_central_lakes_bins <- read.csv("ts_BCdismilarity_central_lakes_bins.csv", row.names = 1)
tsref_BCdismilarity_north_lakes_bins <- read.csv("tsref_BCdismilarity_north_lakes_bins.csv", row.names = 1)
tsref_BCdismilarity_central_lakes_bins <- read.csv("tsref_BCdismilarity_central_lakes_bins.csv", row.names = 1)

## non-binned data
ts_BCdismilarity_north_lakes <- read.csv("ts_BCdismilarity_north_lakes.csv", row.names = 1)
tsref_BCdismilarity_north_lakes <-read.csv("tsref_BCdismilarity_north_lakes.csv", row.names = 1)

## 



umayo <- ts_BCdismilarity_central_lakes_bins %>% 
  filter(id=="umayo") %>%
  filter(age<2010)

#test if BC magnitude is an effect of time
dif <- c(NA, diff(llav$age)) #NA
bc <- llav$timeseriesBC         

cor <- cor.test(dif, bc, method="spearman")
cor

# organise data - depending on lake sequences
bcTS <- plyr::join_all(list(llav, yah, pin, fondo), by='age', type='left') 
colnames(bcTS) <- c("bc_llaviucu", "age", "id", "bc_yahuarcocha", "id", 
                    "bc_pinan", "id", "bc_fondo", "id")
bcTS <- bcTS[, -c(3,5,7,9)]


bcTS <- plyr::join_all(list(llav, yah, pin, fondo, umayo), by='age', type='left') 
colnames(bcTS) <- c("bc_llaviucu", "age", "id", "bc_yahuarcocha", "id", 
                    "bc_pinan", "id", "bc_fondo", "id", "bc_umayo", "id")
bcTS <- bcTS[, -c(3,5,7,9, 11)]


bcTS <- plyr::join_all(list(umayo, llav, yah, pin, fondo), by='age', type='left') 
colnames(bcTS) <- c("bc_umayo", "age", "id", "bc_llaviucu", "id", "bc_yahuarcocha", "id", 
                    "bc_pinan", "id", "bc_fondo", "id")
bcTS <- bcTS[, -c(3,5,7,9, 11)]



#remove Fondococha for now
bcTS2 <- bcTS[, "bc_fondo"] 

#Check NAs and remove rows
row.has.na <- apply(bcTS, 1, function(x){any(is.na(x))})
sum(row.has.na)
bcTS <- bcTS[!row.has.na,]

#Check NAs and remove rows
row.has.na <- apply(df, 1, function(x){any(is.na(x))})
sum(row.has.na)
df <- df[!row.has.na,]


bcTS <- bcTS[,!names(bcTS) %in% c("age")]

#prepare datasets
#this is BC timeseries
x <- llav
y <- yah

seq1 <- merge(x, y, by="age") %>%
  select(timeseriesBC.x, timeseriesBC.y)

#this is BC timeseries from reference
seq1 <- merge(x, y, by="age") %>%
  select(timeseriesBC_fromreference.x, timeseriesBC_fromreference.y)

# 
seq1 <- as.matrix(seq1)
ages <- merge(x, y, by="age") %>% select(age)
ages <- as.matrix(ages)

#modspace <- as.matrix((as.data.frame(modspace)))

# Calculate tau correlation for moving window
window=5 #number of samples of amount of time encompased by the moving window
step=1 # how far to slide the window in each iteration
tsTau<-Tau(seq1,ages,window,step,"sample")



#tsTau<-Tau(modspace,ages,window,step,"sample")

#dev.new(width=5,height=6)
par(mfrow=c(4,1),mar=c(1,2,1,1),oma=c(4,4,2,1))

plot(tsTau$bins[,1],tsTau$Taus,xlim=c(max(ages),min(ages)),
     pch=16, type="b", ylab="Mann-Kendall coeff", xlab="cal yr BP", ylim=c(-1,1))

plot(tsTau$bins[,1],tsTau$Cor,xlim=c(max(ages),min(ages)),
     pch=16, type="b", ylab="Mann-Kendall coeff", xlab="cal yr BP")

# plot datapoints increasing and decreasing tends with tau pval < 0.05 
points(tsTau$bins[,1],tsTau$Taus,xlim=c(max(ages),min(ages)),ylim=c(min(tsTau$Taus),max(tsTau$Taus)),
       pch=16, col=ifelse(tsTau$pval<0.05 & tsTau$Taus>0, "red", 
                          ifelse(tsTau$pval<0.05 & tsTau$Taus<0 , "blue","black")))

points(tsTau$bins[,1],tsTau$Cor,xlim=c(max(ages),min(ages)),ylim=c(min(tsTau$Cor),max(tsTau$Cor)),
       pch=16, col=ifelse(tsTau$cor.p<0.05 & tsTau$Cor>0, "red", 
                          ifelse(tsTau$cor.p<0.05 & tsTau$Cor<0 , "blue","black")))


#points(ages,rep(max(tsTau$Taus)*1.2,length(ages)),pch="|",col="red")
text(max(ages),max(tsTau$Taus)*0.95,paste("windows=",window,"steps=",step),pos=4)
text(max(ages),-0.10,"Pinan-Fondococha",pos=4)
text(max(ages),max(tsTau$Taus)*0.25,"p < 0.05; tau>0",pos=4, col="red")
text(max(ages),max(tsTau$Taus)*0.05,"p < 0.05; tau<0",pos=4, col="blue")


# extract results Tau
altN <- data.frame(cbind(tsTau$bins[,1], tsTau$Taus, tsTau$pval))
colnames(altN) <- c("age", "tau", "pval")
altN$seq <- "Lla-Yah"

altS <- data.frame(cbind(tsTau$bins[,1], tsTau$Taus, tsTau$pval))
colnames(altS) <- c("age", "tau", "pval")
altS$seq <- "Lla-Fondo"

latS <- data.frame(cbind(tsTau$bins[,1], tsTau$Taus, tsTau$pval))
colnames(latS) <- c("age", "tau", "pval")
latS$seq <- "Lla-Pin"

latN <- data.frame(cbind(tsTau$bins[,1], tsTau$Taus, tsTau$pval))
colnames(latN) <- c("age", "tau", "pval")
latN$seq <- "Yah-Fondo"

pair5 <- data.frame(cbind(tsTau$bins[,1], tsTau$Taus, tsTau$pval))
colnames(pair5) <- c("age", "tau", "pval")
pair5$seq <- "Yah-Pin"

pair6 <- data.frame(cbind(tsTau$bins[,1], tsTau$Taus, tsTau$pval))
colnames(pair6) <- c("age", "tau", "pval")
pair6$seq <- "Fon-Pin"

# extract results Cor pearson
altN <- data.frame(cbind(tsTau$bins[,1], tsTau$Cor, tsTau$cor.p.ad))
colnames(altN) <- c("age", "tau", "pval")
altN$seq <- "Lla-Yah"

altS <- data.frame(cbind(tsTau$bins[,1], tsTau$Cor, tsTau$cor.p.ad))
colnames(altS) <- c("age", "tau", "pval")
altS$seq <- "Lla-Fondo"

latS <- data.frame(cbind(tsTau$bins[,1], tsTau$Cor, tsTau$cor.p.ad))
colnames(latS) <- c("age", "tau", "pval")
latS$seq <- "Lla-Pin"

latN <- data.frame(cbind(tsTau$bins[,1], tsTau$Cor, tsTau$cor.p.ad))
colnames(latN) <- c("age", "tau", "pval")
latN$seq <- "Yah-Fondo"

pair5 <- data.frame(cbind(tsTau$bins[,1], tsTau$Cor, tsTau$cor.p.ad))
colnames(pair5) <- c("age", "tau", "pval")
pair5$seq <- "Yah-Pin"

pair6 <- data.frame(cbind(tsTau$bins[,1], tsTau$Cor, tsTau$cor.p.ad))
colnames(pair6) <- c("age", "tau", "pval")
pair6$seq <- "Fon-Pin"


# make tidy
alt <- rbind(altN, altS, latN, latS, pair5, pair6)
write.csv(alt, "pairwise_spearman_BC_ref.csv")



##read pairwise Spearman rank correlations
pairwise_spearman_BC <- read.csv("pairwise_spearman_BC.csv", row.names = 1)


