##BC dissimilarity time series

## BC dissimilarity from consecutive and cumulative samples
#successive samples
## non-binned data
ts_BCdismilarity_north_lakes <- read.csv("data/ts_BCdismilarity_north_lakes.csv", row.names = 1)
tsref_BCdismilarity_north_lakes <- read.csv("data/tsref_BCdismilarity_north_lakes.csv", row.names = 1)


llav <- ts_BCdismilarity_north_lakes %>% filter(id=="llaviucu")
yah <- ts_BCdismilarity_north_lakes %>% filter(id=="yahuarcocha")
pin <- ts_BCdismilarity_north_lakes %>% filter(id=="pinan")
fondo <- ts_BCdismilarity_north_lakes %>% filter(id=="fondococha")

llav_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="llaviucu")
yah_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="yahuarcocha")
pin_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="pinan")
fondo_ref <- tsref_BCdismilarity_north_lakes %>% filter(id=="fondococha")

#read in fossil distance to modern analogues
dist_to_analogues_all_lakes<- read.csv("data/dist_to_analogues_all_lakes.csv",row.names = 1) 

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
long_matrixHumanHist <- read.csv("data/long_matrixHumanHist.csv", row.names=1) %>%
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