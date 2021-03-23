#----------------------------------------------------------------------------------------------------
# Script: Diatom community dissimilarity time series: cumulative and successive bray-curtis distances
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#----------------------------------------------------------------------------------------------------

##loading libraries for functions used
library(analogue) 
library(tidyverse) 
library(ggplot2)
library(viridis)
library(zoo)

##Read in core and trainingset diatom abundances (absolute counts)
df <- readRDS("data/coresList.rds")

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  # comment the next line when predicting core trajectories in the timetrack analysis
  core <- core[, colSums(core) > 0] #select only present species
  return(core)
}

#Wrap-up the function over the list
cores <- lapply(seq_along(df), remove, cores=df)

#name list elements
names(cores) <- c("Fondococha", "Lagunillas", "Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

#drop trainingset from the core list
cores$trainingset <- NULL

#Calculate BC dissimilarity from successive samples
doBCsuccessive <- function(i, cores,..) {
  core <- cores[[i]]
  D <- analogue::distance(core/100, method="bray")
  # create an indicator for all diagonals in the matrix
  d <- row(D) - col(D)
  # use split to group on these values
  diagonal<-split(D, d)
  timeseriesBC<-unlist(diagonal["1"]) # select relevant one (diag = 1)
  upper_age <- df[[i]]$upper_age
  lower_age <- df[[i]]$lower_age
  cbind.data.frame(timeseriesBC, upper_age[-1], lower_age[-1]) 
}

#Wrap-up function over the cores list and name the cores
coresBCsucc <- lapply(seq_along(cores), doBCsuccessive, cores=cores)
names(coresBCsucc) <- names(cores)

#Extract core dataframes from the list
ts_BCdissimilarity <- plyr::ldply(coresBCsucc, data.frame)
ts_BCdissimilarity[,c(4,5)] <- NULL
colnames(ts_BCdissimilarity) <- c("lake", "timeseriesBC", "age")

# Select lakes of study
select.lakes <- paste(c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha"), collapse = '|')
ts_BCdissimilarity <- ts_BCdissimilarity %>% filter(str_detect(lake, select.lakes))

write.csv(ts_BCdissimilarity, "data/ts_BCdissimilarity.csv")

#Calculate BC dissimilarity from the reference assemblage (i.e. cumulative)
doBCcumulative <- function(i, cores,..) {
  core <- cores[[i]]
  D <- analogue::distance(core/100, method="bray")
  BCcum <-D[,ncol(D)] #take the oldest sample as reference to see the direction of change
  upper_age <- df[[i]]$upper_age
  lower_age <- df[[i]]$lower_age
  cbind.data.frame(BCcum, upper_age, lower_age) 
}

#Wrap-up function over the cores list and name the cores
coresBCcum <- lapply(seq_along(cores), doBCcumulative, cores=cores)
names(coresBCcum) <- names(cores)

#Calculate 5th percentile cut-off value of the cumulative BC dissimilarity time series for each lake
fifthpercentile <- list()
for (i in 1:length(coresBCcum)) {
  fifthpercentile[[i]] <- quantile(coresBCcum[[i]]$BCcum, probs = c(0.75, 0.95))
}
names(fifthpercentile) <- names(cores)

# Select lakes of study and extract from the list
lakes <- c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha")
fifthpercentile <- fifthpercentile[lakes]
coresBCcum <- coresBCcum[lakes]

##plot BC dissimilarity from the reference conditions
par(mfrow = c(2, 2))
par(mar = c(2.5, 3.5, 1, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(0, 0, 3, 0))

for (i in 1:length(coresBCcum)) {
  plot.ts(coresBCcum[[i]]$BCcum, cex = 0.6, cex.axis = 0.8,
          las = 1, pch = 19, main=names(coresBCcum[i]), ylim=c(0,1),
          ylab="BC Dissimilarity")
  abline(h=fifthpercentile[[i]][2], col="blue")
}

#Extract core dataframes from the list
ts_BCdissimilarity_ref <- plyr::ldply(coresBCcum, data.frame)
ts_BCdissimilarity_ref[,4] <- NULL #remove lower age variable
colnames(ts_BCdissimilarity_ref) <- c("lake", "timeseriesBC_fromreference", "age")

write.csv(ts_BCdissimilarity_ref, "data/ts_BCdissimilarity_ref.csv")

##BC dissimilarity time series
llav <- ts_BCdissimilarity %>% filter(lake=="Llaviucu")
yah <- ts_BCdissimilarity %>% filter(lake=="Yahuarcocha")
pin <- ts_BCdissimilarity %>% filter(lake=="Pinan")
fondo <- ts_BCdissimilarity %>% filter(lake=="Fondococha")

llav_ref <- ts_BCdissimilarity_ref %>% filter(lake=="Llaviucu")
yah_ref <- ts_BCdissimilarity_ref %>% filter(lake=="Yahuarcocha")
pin_ref <- ts_BCdissimilarity_ref %>% filter(lake=="Pinan")
fondo_ref <- ts_BCdissimilarity_ref %>% filter(lake=="Fondococha")

#read in fossil distance to modern analogues
dist_to_analogues_all_lakes<- read.csv("data/dist_to_analogues_all_lakes.csv",row.names = 1) 

## Set parameters for plotting (Figure 3)
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
abline(h=fifthpercentile$Pinan[2], col="black", lty=2, lwd=1) #5th cutoff value
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
abline(h=goodpoorbad$Yahuarcocha[2], col="black", lty=2,lwd=1) #5th cutoff value
lines(yah_ref$age, yah_ref$timeseriesBC_fromreference, lwd=1, col=seq_palette[2])
for(i in 1:4) abline(v=zones_yahuarcocha[i], col=i+1, lty=2, lwd=2)
dist_to_analogues_yah <- dist_to_analogues_all_lakes %>% 
  filter(lake=="Yahuarcocha") %>%
  filter(quality=="bad")
points(dist_to_analogues_yah$ages,rep(max(dist_to_analogues_yah$dist_to_analogues)*1.10,
                                      length(dist_to_analogues_yah$ages)),pch=25,bg="darkorange", col="black", cex=1.5)
legend("topright",legend = c("Yahuarcocha"), bg="transparent", bty="n",inset=0.015)

plot(fondo$age[1:79], fondo$timeseriesBC[1:79], type="n", ylab="BC\ndissimilarity", xlab=NA,xaxt='n',
     cex.lab=label.V,cex.axis=axis.V, ylim = c(0,1), xlim=c(-60,2000))
abline(v=c(1996:2063), col="grey", lty=1, lwd=3)
abline(h=goodpoorbad$Fondococha[2], col="black", lty=2,lwd=1) #5th cutoff value
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
abline(h=goodpoorbad$Llaviucu[2], col="black", lty=2,lwd=1) #5th cutoff value
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
long_matrixHumanHist <- read.csv("data/HYDE/long_matrixHumanHist.csv", row.names=1) %>%
  mutate(age=1950-years)
var <- "cropland"

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

#read pages2k temperature data
long_SAtempdata <- read.csv("data/paleoclimate/long_SAtempdata.csv", row.names = 1) %>%
  filter(method=="LNA_anomalies")

plot(long_SAtempdata$age_calBP, long_SAtempdata$Temp, type="l", lwd=2, ylab="Temperature\n\u00B0C",
     cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col="black", xlab=NA,xaxt='n')
legend("topright",legend = c("PAGES2K Temperature anomalies reconstruction"), 
       bg="transparent", bty="n",inset=0.010,cex=0.9)

### read South America hydroclimate records
pumacocha <- read.csv("data/paleoclimate/SouthAmerica_hydroclimate_records.csv", row.names = 1) %>%
  filter(core %in% c("pumacocha")) %>% filter(!d13C==-9999.00)

plot(pumacocha$age_calBP, pumacocha$d18O, type="l", lwd=2, ylab=expression(paste (delta^18, "O \u2030")),
     cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col="black", xlab=NA,xaxt='n')
legend("topright",legend = c("Pumacocha lake"), bg="transparent", bty="n",inset=0.010,cex=0.9)

#read archaecological time series
ecuadorspds <- read.csv("data/archaeologicalsites/ecuadorspd_ts.csv", row.names = 1)

plot(ecuadorspds$calBP, ecuadorspds$PrDens, type="h", lwd=2, ylab="SPD",
     cex.lab=label.V,cex.axis=axis.V, xlim=c(-60,2000), col="grey")
legend("topright", legend=c("Ecuador radiocarbon dates"), bg="transparent", bty="n",inset=0.010,cex=0.9)

mtext("Years Before Present",side=1,outer=T,line=1,  adj=0.525,cex=1)
dev.off()