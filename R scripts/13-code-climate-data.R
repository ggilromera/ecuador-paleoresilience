setwd("/Volumes/xbenitogranell-data/0_project/data/climate")
#setwd("/nfs/xbenitogranell-data/0_project/data/climate")

library(dplyr)
library(tidyverse)
library(cowplot)
library(gtools)
library(sp)
library(raster)
library(ggvegan)
library(usdm)

#read geographical coordinates of hydrological records (from Andes metadatabase)
andes_clim_records <- read.csv("andes_climatic_records.csv") %>%
  filter(record %in% c("Ubaque", "Huagapo", "Pumacocha", "Palestina cave"))

# Read climate data
ubaque <- read.csv("bird-13C-ubaque.csv")
huagapo <- read.csv("data-huagapo.csv")
palestina <- read.csv("data-palestina.csv")
tigre <- read.csv("data-tigre.csv") %>% rename(d18O=tigre_perdido_d18O.PDB)
junco <- read.csv("data-junco.csv") %>% rename(age_AD=Age..cal.yrs.BP.)
pumacocha <- read.csv("data-pumacocha.csv")
quelccaya <- read.csv("data-quelccaya.csv")
umayo <- read.csv("data-umayo.csv")
triumfo <- read.csv("data-triumfo.csv")
cariaco <- read.csv("data-cariaco.csv") %>% rename(years_AD = age_AD) 
yanacocha <- read.csv("data-yanacocha.csv")
aculeo <- read.csv("data-aculeo.csv") %>% rename(years_AD = age_AD) 
lutacocha <- read.csv("data-lutacocha.csv") 
jahuacocha <- read.csv("data-jahuacocha.csv") 
SAproxydata <- read.csv("pages2k-proxydata.csv") %>%
  mutate(age_calBP=(1950-Year.C.E.)) %>%
  dplyr::select(-Year.C.E.)
SAtempdata <- read.csv("pages2k-tempdata.csv") %>% 
  mutate(age_calBP=(1950-Year.CE)) %>%
  dplyr::select(-Year.CE)


# Transform data for that all have agecalyrBP
ubaque$years_AD <- round((ubaque$bchron_age_model*(-1)+1950), digits=0) # transform so that time moves in the right direction, and 0 is year 0 of the common era
ubaque$age_calBP <- (1950-(ubaque$years_AD)) # transform so that time moves in the right direction, and 0 is year 0 of the common era

umayo$years_AD <- round((umayo$upper_age*(-1)+1950), digits=0) # transform so that time moves in the right direction, and 0 is year 0 of the common era
umayo$age_calBP <- round(umayo$upper_age, digits = 0)

palestina$age_calBP <- (1950-(palestina$years_AD))
tigre$age_calBP <- round((tigre$age_AD*(-1)+1950), digits=0)
junco$age_calBP <- round((junco$age_AD*(-1)+1950), digits=0)
quelccaya$age_calBP <- (1950-(quelccaya$years_AD))
pumacocha$age_calBP <- round(pumacocha$age_calBP, digits=0)
triumfo$years_AD <- round((triumfo$age_calBP*(-1)+1950), digits=0)
cariaco$age_calBP <- (1950-(cariaco$years_AD))
yanacocha$years_AD <- round((yanacocha$age_calBP*(-1)+1950), digits=0)
lutacocha$years_AD <- round((lutacocha$age_calBP*(-1)+1950), digits=0)
jahuacocha$years_AD <- round((jahuacocha$age_calBP*(-1)+1950), digits=0)
aculeo$age_calBP <- (1950-(aculeo$years_AD))


#create  a new variable containing core name
ubaque$core <- c(rep("ubaque", nrow(ubaque)))
ubaque <- ubaque %>% dplyr::select(years_AD, age_calBP, d13C, core)
ubaque$long <- -73.935
ubaque$lat <- 4.499

huagapo$core <- c(rep("huagapo", nrow(huagapo)))
huagapo <- huagapo %>% dplyr::select(years_AD, age_calBP, d18O, d13C, core)

palestina$core <- c(rep("palestina", nrow(palestina)))
palestina <- palestina %>% dplyr::select(years_AD, age_calBP, d18O, d13C, core)

tigre$core <- c(rep("tigreperdido", nrow(tigre)))

junco$core <- c(rep("junco", nrow(junco)))
junco <- junco %>% dplyr::select(age_calBP, Sand, Clay, Silt, core) 

pumacocha$core <- c(rep("pumacocha", nrow(pumacocha)))
pumacocha <- pumacocha %>% dplyr::select(years_AD, age_calBP, d18O, d13C, core)

quelccaya$core <- c(rep("quelccaya", nrow(quelccaya)))
quelccaya <- quelccaya %>% dplyr::select(years_AD, age_calBP, d18O, core)

umayo$core <- c(rep("umayo", nrow(umayo)))
umayo <- umayo %>% dplyr::select(years_AD, age_calBP, d18O, core)

triumfo$core <- c(rep("triumfo", nrow(triumfo)))
triumfo <- triumfo %>% dplyr::select(years_AD, age_calBP, d13C, core)

cariaco$core <- c(rep("cariaco", nrow(cariaco)))

yanacocha$core <- c(rep("yanacocha", nrow(yanacocha))) 
yanacocha <- yanacocha %>% mutate(CaTi = Ca..x/Ti..x) %>% dplyr::select(years_AD, age_calBP,CaTi,core)

aculeo$core <- c(rep("aculeo", nrow(aculeo)))

lutacocha$core <- c(rep("lutacocha", nrow(lutacocha)))
lutacocha <- lutacocha %>% mutate(CaTi = Ca..x/Ti..x) %>% dplyr::select(years_AD, age_calBP,CaTi,core)

jahuacocha$core <- c(rep("jahuacocha", nrow(jahuacocha)))
jahuacocha <- jahuacocha %>% mutate(CaTi = Ca.cps/Ti.cps) %>% dplyr::select(years_AD, age_calBP,CaTi,core)


#Merge dataframes
df <- bind_rows(huagapo, palestina, pumacocha, quelccaya, ubaque, umayo, triumfo, cariaco, yanacocha, aculeo,
                lutacocha, jahuacocha, tigre, junco)

#write.csv(df, "SouthAmerica_hydroclimate_records.csv")

# Create list
climateList <- split(df, df$core)


#Plot individual climate time series
ubaque <- df %>% filter(core=="ubaque") %>% filter(!d13C==-9999.00)
plt_ubaq <- ggplot(ubaque) + geom_line(aes(age_calBP, d13C)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.line.x = element_blank(), 
         axis.title.y = element_text(size=9),
        plot.margin = unit(c(0,0,0,1.5),"cm")) +
  theme(legend.position = "none") + 
  ylab(expression(paste (delta^13, "C \u2030")))+
  ggtitle("Ubaque") 
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  # #annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate("text", x = 1000, y = -20, label = "italic(MCA)", parse=TRUE) +
  # annotate("text", x = 400, y = -20, label = "italic(LIA)",parse=TRUE)
  

palestina <- df %>% filter(core=="palestina")
plt_pales <- ggplot(palestina) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.y = element_text(size=9)) +
  theme(legend.position = "none") + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=9),
        plot.margin = unit(c(0,0,0,1.5),"cm")) +
  ylab(expression(paste (delta^18, "O \u2030")))+
  ggtitle("Palestina") 
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") 
  #annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")
  # annotate("text", x = 1000, y = -5.6, label = "italic(MCA)", parse=TRUE) +
  # annotate("text", x = 400, y = -5.6, label = "italic(LIA)",parse=TRUE)

huagapo <- df %>% filter(core=="huagapo")
plt_huag <- ggplot(huagapo) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=9),
        plot.margin = unit(c(0.5,0.5,0.5,1.5),"cm")) +
  theme(legend.position = "none") + 
  ylab(expression(paste (delta^18, "O \u2030"))) +
  ggtitle("Huagapo")
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff")
  #annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")+
  # annotate("text", x = 1000, y = -12.1, label = "italic(MCA)", parse=TRUE) +
  # annotate("text", x = 400, y = -12.1, label = "italic(LIA)",parse=TRUE)


pumacocha <- df %>% filter(core=="pumacocha")
plt_pumac <- ggplot(pumacocha) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  # theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.line.x = element_blank(),
  theme(axis.title.y = element_text(size=9),
        plot.margin = unit(c(0,0,0,1.5),"cm")) +
  theme(legend.position = "none") +
  ggtitle("Pumacocha") +
  xlab("Cal years BP") +
  ylab(expression(paste (delta^18, "O \u2030")))
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") 
  #annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")+
  # annotate("text", x = 1000, y = -10.3, label = "italic(MCA)", parse=TRUE) +
  # annotate("text", x = 400, y = -10.3, label = "italic(LIA)",parse=TRUE)

tigre <- df %>% filter(core=="tigreperdido")
plt_tigre <- ggplot(tigre) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  # theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.line.x = element_blank(),
  theme(axis.title.y = element_text(size=9),
        plot.margin = unit(c(0,0,0,1.5),"cm")) +
  theme(legend.position = "none") +
  ggtitle("Tigre Perdido") +
  xlab("Cal years BP") +
  ylab(expression(paste (delta^18, "O \u2030")))+
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") 
#annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")+
# annotate("text", x = 1000, y = -10.3, label = "italic(MCA)", parse=TRUE) +
# annotate("text", x = 400, y = -10.3, label = "italic(LIA)",parse=TRUE)

junco <- df %>% filter(core=="junco")
plt_junco <- ggplot(junco) + geom_line(aes(age_calBP, Sand)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  # theme(axis.title.x=element_blank(),
  #        axis.text.x=element_blank(),
  #        axis.line.x = element_blank(),
  theme(axis.title.y = element_text(size=9),
        plot.margin = unit(c(0,0,0,1.5),"cm")) +
  theme(legend.position = "none") +
  ggtitle("Junco") +
  xlab("Cal years BP") +
  ylab("Sand (%)")
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") 
#annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")+
# annotate("text", x = 1000, y = -10.3, label = "italic(MCA)", parse=TRUE) +
# annotate("text", x = 400, y = -10.3, label = "italic(LIA)",parse=TRUE)


quelccaya <- df %>% filter(core=="quelccaya")
plt_quelcc <- ggplot(quelccaya) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=7.5),
        plot.margin = unit(c(0.5,0.5,0.5,1.5),"cm")) +
  theme(legend.position = "none") + 
  ylab(expression(paste (delta^18, "O \u2030")))
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  # annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")


umayo <- df %>% filter(core=="umayo") 
plt_umay <- ggplot(umayo) + geom_line(aes(age_calBP, d18O)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=7.5)) +
  theme(legend.position = "none") + 
  ylab(expression(paste (delta^18, "O \u2030")))+
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")



triumfo <- df %>% filter(core=="triumfo") %>% filter(!d13C==0)
plt_triumfo <- ggplot(triumfo) + geom_line(aes(age_calBP, d13C)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=7.5),
        plot.margin = unit(c(0,0,0,1.5),"cm")) +
  theme(legend.position = "none") + 
  ylab(expression(paste (delta^13, "C \u2030")))+
  ggtitle("Triunfo") 
  # annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  # annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  # annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")


cariaco <- df %>% filter(core=="cariaco") 
plt_cari <- ggplot(cariaco) + geom_line(aes(age_calBP, temperature)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,0)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=7.5),
        plot.margin = unit(c(0.5,0.5,0.5,1.5),"cm")) +
  theme(legend.position = "none") + 
  ylab("Mg/Ca") +
  ggtitle("Cariaco Basin") +
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")


yanacocha <- df %>% filter(core=="yanacocha")
plt_yanac <- ggplot(yanacocha) + geom_line(aes(age_calBP, CaTi)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,-40)) +
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.line.x = element_blank(), 
  #       axis.title.y = element_text(size=7.5)) +
  theme(legend.position = "none") + 
  ylab("Ca/Ti")+ xlab("Cal yr BP") +
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")



lutacocha <- df %>% filter(core=="lutacocha")
plt_lutac <- ggplot(lutacocha) + geom_line(aes(age_calBP, CaTi)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,-60)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=7.5)) +
  theme(legend.position = "none") + 
  ylab("Ca/Ti") +
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")



jahuacocha <- df %>% filter(core=="jahuacocha")
plt_jahua <- ggplot(jahuacocha) + geom_line(aes(age_calBP, CaTi)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,-60)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.line.x = element_blank(), 
        axis.title.y = element_text(size=7.5)) +
  theme(legend.position = "none") + 
  ylab("Ca/Ti") +
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")


aculeo <- df %>% filter(core=="aculeo")
plt_acul <- ggplot(aculeo) + geom_line(aes(age_calBP, temperature)) + 
  #scale_color_manual(values =c("#88101A", "#DB291A")) +
  theme_classic() + xlim(c(3000,-60)) +
  theme(axis.title.y = element_text(size=7.5)) +
  theme(legend.position = "none") + 
  ylab("Temperature (deg C)") + xlab("Cal yr BP") +
  annotate(geom="rect", xmin=1100, xmax=900, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D") +
  annotate(geom="rect", xmin=600, xmax=180, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#2b44ff") +
  annotate(geom="rect", xmin=150, xmax=0, ymin=-Inf, ymax=Inf, alpha=0.2, fill="#F8766D")


#plot records arranged by latitude
plot_hydroclimate_proxies <- plot_grid(plt_ubaq, plt_triumfo, plt_pales, plt_pumac, ncol = 1, align = "hv")
plot_hydroclimate_proxies

ggsave("plotTS_proxy_2.png", plot_hydroclimate_proxies, height = 8, width = 10)

# Plot pages2k temperature data
long_SAproxydata <- gather(SAproxydata, key=record, value=Temp, -age_calBP)
long_SAtempdata <- gather(SAtempdata, key=method, value=Temp, -age_calBP)

ggplot(long_SAproxydata, aes(x=age_calBP,y=Temp)) + 
  geom_line() +
  facet_wrap(~record, scales="free_y")

ggplot(long_SAtempdata, aes(x=age_calBP,y=Temp, colour=method, group=method)) + 
  geom_line() +
  ylab("Temperature anomalies degC")

write.csv(long_SAproxydata, "long_SAproxydata.csv")
write.csv(long_SAtempdata, "long_SAtempdata.csv")


# Plot whole climate time series
plot1 <- df %>%
  dplyr::select(years_AD, age_calBP, d18O, d13C, core) %>%
  filter(years_AD >= "0") %>%
  ggplot() +
  geom_line(aes(x = age_calBP, y = d18O)) +
  facet_wrap(~ core, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

plot1
ggsave("proxydataTS.png", plot1, height = 8, width = 10)



## PAGES2k v2.1.0 database
library(lipdR)
D <- readLipd("/Volumes/xbenitogranell-data/0_project/data/PAGES2K-southamerica")
#D <- readLipd("/nfs/xbenitogranell-data/0_project/data//PAGES2K-southamerica")

TS <- extractTs(D)
mts <- filterTs(TS, "paleoData_useInGlobalTemperatureAnalysis == TRUE")
tidyData <- tidyTs(mts)

plot.df <- tidyData %>% 
  filter(between(year,1000,2000)) %>% #only years from 1600 to 2000
  filter(interpretation1_variable == "T") %>% #only variables sensitive temperature
  #filter(archiveType=="lake sediment") #only lake sediment
  group_by(paleoData_TSid) %>% #group by column
  arrange(archiveType) #and sort by archiveType

plotTimeseriesStack(plot.df)

plotTimeseriesStack(plot.df,colorVar =  "archiveType")

plotTimeseriesStack(plot.df,labSize = 2,colorRamp = c("coral","black","blue","magenta","dark green"))




## Simulated paleoclimatic grids

#Temperature Seasonality 3000BP-present (1950)
setwd("/nfs/xbenitogranell-data/0_project/data/climate/TempSeasonality_3000BP")
setwd("/nfs/xbenitogranell-data/0_project/data/climate/PrecipSeasonality_3000BP/25Nov2019_0359PM_52.051s")
setwd("/nfs/xbenitogranell-data/0_project/data/climate/MAP_3000BP")
setwd("/nfs/xbenitogranell-data/0_project/data/climate/MAT_3000BP")


#
fs1 <- list.files(pattern = ".asc") # Locate files

os1 <- mixedsort(fs1) # Sort chronologically
os1 <- rev(os1) 
rs1 <- lapply(os1, raster) # Convert to raster
ss1 <- stack(rs1) # Stack rasters
proj4string(ss1) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")


# # ss1 <- crop(ss1, extent(-90, -30, -60, 20))
# # 
# r.maxs1 = cellStats(ss1,"max") # Cell statistics
# r.mins1 = cellStats(ss1, "min")
# #
# rs.scale1 <- ((ss1 - r.mins1) / (r.maxs1 - r.mins1) - 0.5 ) * 2 #Rescale from -1 to 1
# 
# vs1 <- calc(rs.scale1, fun=var) #Calculate variance in stack
# crs(vs1) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # Project
# #vs1 <- projectRaster(vs1, crs="+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")
# vas1 <- disaggregate(vs1, fact = 20) # Increase resolution
# #gs1 <- as(vas1, 'SpatialGridDataFrame')


# read lake sites
coorSites <- read.csv("/nfs/xbenitogranell-data/0_project/data/historic/coorSitesNAndes.csv", row.names = 1)
colnames(coorSites) <- c("sites", "Longitude", "Latitude")
#add Pinan
coorSites <- add_row(coorSites, sites="EpNGEO-J_Piñan1", Longitude = -78.444, Latitude = 0.506)
coorSites_df <- coorSites

coordinates(coorSites) <- c('Longitude', 'Latitude')

#put on the same geo-referencial
proj4string(coorSites) <- CRS("+proj=longlat +datum=WGS84 +units=m +ellps=WGS84 +towgs84=0,0,0")

# extract data
sitesClimate <- data.frame(coorSites@data[["sites"]],raster::extract(ss1, coorSites)) 

# Save results
  sitesTempSeasonality <- sitesClimate
  #colnames(sitesTempSeasonality) <- gsub("grid_data", "TempSeason", colnames(sitesTempSeasonality))
  
  colnames(sitesTempSeasonality) <- gsub("grid_data_", "", colnames(sitesTempSeasonality))
  colnames(sitesTempSeasonality) <- gsub("BP", "", colnames(sitesTempSeasonality))
  
  times <- as.character(seq(from=0, to=3000, by=30)) # Define all needed time steps
  sitesTempSeasonalityT <- as.data.frame(t(sitesTempSeasonality))[-1,]
  rownames(sitesTempSeasonalityT) <- times
  colnames(sitesTempSeasonalityT) <- as.character(sitesTempSeasonality[,1])
  #sitesTempSeasonalityT$variable <- "TempSeasonality"
  #sitesTempSeasonalityT$variable <- NULL
  indx <- sapply(sitesTempSeasonalityT, is.factor)
  sitesTempSeasonalityT[indx] <- lapply(sitesTempSeasonalityT[indx], function(x) as.numeric(as.character(x)))
  

##
  sitesPrecipSeasonality <- sitesClimate
  #colnames(sitesPrecipSeasonality) <- gsub("grid_data", "PrecipSeason", colnames(sitesPrecipSeasonality))
  colnames(sitesPrecipSeasonality) <- gsub("grid_data_", "", colnames(sitesPrecipSeasonality))
  colnames(sitesPrecipSeasonality) <- gsub("BP", "", colnames(sitesPrecipSeasonality))
  
  times <- as.character(seq(from=0, to=3000, by=30)) # Define all needed time steps
  sitesPrecipSeasonalityT <- as.data.frame(t(sitesPrecipSeasonality))[-1,]
  rownames(sitesPrecipSeasonalityT) <- times
  colnames(sitesPrecipSeasonalityT) <- as.character(sitesPrecipSeasonality[,1])
  #sitesPrecipSeasonalityT$variable <- "PrecipSeasonality"
  #sitesPrecipSeasonalityT$variable <- NULL
  indx <- sapply(sitesPrecipSeasonalityT, is.factor)
  sitesPrecipSeasonalityT[indx] <- lapply(sitesPrecipSeasonalityT[indx], function(x) as.numeric(as.character(x)))
  
  
##
  sitesMAP <- sitesClimate
  #colnames(sitesMAP) <- gsub("grid_data", "MAP", colnames(sitesMAP))
  colnames(sitesMAP) <- gsub("grid_data_", "", colnames(sitesMAP))
  colnames(sitesMAP) <- gsub("BP", "", colnames(sitesMAP))
  
  times <- as.character(seq(from=0, to=3000, by=30)) # Define all needed time steps
  sitesMAPT <- as.data.frame(t(sitesMAP))[-1,]
  rownames(sitesMAPT) <- times
  colnames(sitesMAPT) <- as.character(sitesMAP[,1])
  #sitesMAPT$variable <- "MAP"
  #sitesMAPT$variable <- NULL
  indx <- sapply(sitesMAPT, is.factor)
  sitesMAPT[indx] <- lapply(sitesMAPT[indx], function(x) as.numeric(as.character(x)))
  
  
##  
  sitesMAT <- sitesClimate
  #colnames(sitesMAT) <- gsub("grid_data", "MAT", colnames(sitesMAT))
  colnames(sitesMAT) <- gsub("grid_data_", "", colnames(sitesMAT))
  colnames(sitesMAT) <- gsub("BP", "", colnames(sitesMAT))
  
  times <- as.character(seq(from=0, to=3000, by=30)) # Define all needed time steps
  sitesMATT <- as.data.frame(t(sitesMAT))[-1,]
  rownames(sitesMATT) <- times
  colnames(sitesMATT) <- as.character(sitesMAT[,1])
  sitesMATT$variable <- "MAT"
  sitesMATT$variable <- NULL
  indx <- sapply(sitesMATT, is.factor)
  sitesMATT[indx] <- lapply(sitesMATT[indx], function(x) as.numeric(as.character(x)))
  
 
df <- analogue::join(sitesTempSeasonalityT, sitesPrecipSeasonalityT, sitesMAPT, sitesMATT)

## Bin climatic data
binProxyData <- function(i, cores, ...) {
  core <- cores[[i]]
  years <- rownames(core) #this is to get ages_AD
  
  # Function for binning samples in each core (Alistair Seddon)
  source("/nfs/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
  
  #make the age categories
  diatomBin1 <- binFunc(as.data.frame(core), as.numeric(years), 100, 0, 1700) ##
  diatomBin2 <- binFunc(as.data.frame(core), as.numeric(years), 30, 1700, 2000) ##
  
  #merge the two binning dataframes for North lakes
  rbind.data.frame(diatomBin1, diatomBin2)
}

#wrap up the function
binnedSimulated <- lapply(seq_along(df), binProxyData, cores = df)
names(binnedSimulated) <- names(df)

nams <- names(binnedSimulated)
for (i in seq_along(binnedSimulated)) {
  assign(paste0("", nams[i]), binnedSimulated[[i]])
}

select.lakes <- paste(c("Piñan1", "Yahurcch", "Fondocch1", "Llaviucu"), collapse = '|')
lakes_PrecipSeason <- sitesPrecipSeasonalityT[,str_detect(colnames(sitesPrecipSeasonalityT),
                                            select.lakes)]
lakes_PrecipSeason$variable <- "PrecipSeason"
  
lakes_TempSeason <- sitesTempSeasonalityT[,str_detect(colnames(sitesTempSeasonalityT),
                                                      select.lakes)]
lakes_TempSeason$variable <- "TempSeason"

lakes_MAP <- sitesMAPT[,str_detect(colnames(sitesMAPT),select.lakes)]
lakes_MAP$variable <- "MAP"

lakes_MAT <- sitesMATT[,str_detect(colnames(sitesMATT),select.lakes)]
lakes_MAT$variable <- "MAT"

#rbind data
climateDataSites <- rbind(lakes_PrecipSeason, lakes_TempSeason, lakes_MAP, lakes_MAT)
indx <- sapply(climateDataSites, is.factor)
climateDataSites[indx] <- lapply(climateDataSites[indx], function(x) as.numeric(as.character(x)))
climateDataSitesBinned <- split(climateDataSites, climateDataSites$variable)

climateDataSitesBinned$MAP$variable <- NULL
climateDataSitesBinned$MAT$variable <- NULL
climateDataSitesBinned$PrecipSeason$variable <- NULL
climateDataSitesBinned$TempSeason$variable <- NULL

#do interpolation between adjacent samples
# library(zoo)
# doInterpol <- function(i, cores, ...) {
#   core <- cores[[i]]
#   core <- as.data.frame(na.approx(core, na.rm = TRUE)) #do interpolation between adjacent samples
#   row.names(core) <- as.numeric(row.names(climateDataSitesBinned$PrecipSeason))
#   core
# }
# 
# climateDataSitesBinned_interPol <- lapply(seq_along(climateDataSitesBinned), doInterpol, cores = climateDataSitesBinned)
# names(climateDataSitesBinned_interPol) <- names(climateDataSitesBinned)


#save binnedSimulated Climate list
saveRDS(file="/nfs/xbenitogranell-data/0_project/data/historic/CommonEra/binnedClimateSimulated.rds", climateDataSitesBinned)

#Read in binnedSimulated Climate list
binnedClimateSimulated <- readRDS("/Volumes/xbenitogranell-data/0_project/data/historic/CommonEra/binnedClimateSimulated.rds")

matrixClimateHist <- plyr::ldply(binnedClimateSimulated, data.frame)

#this is to extract binned ages; only needs to do it once because binned ages are the same for all the cores
ages <- as.data.frame(as.numeric(row.names(binnedClimateSimulated$PrecipSeason)))
ages <- plyr::rename(ages,c("as.numeric(row.names(binnedClimateSimulated$PrecipSeason))"="age"))

matrixClimHist <- cbind(matrixClimateHist, ages)

long_matrixClimHist<-gather(matrixClimHist, key=lake, value=value, -age, -.id)
head(long_matrixClimHist)

select.lakes <- paste(c("Yahurcch", "Fondocch1"), collapse = '|')
long_matrixClimHist <- long_matrixClimHist %>%
  filter(str_detect(lake, select.lakes))


## plot the data
plt_clim_ts <- ggplot(long_matrixClimHist, aes(x = age, y = value, colour=lake)) +
  geom_line()+
  facet_wrap(~ .id, ncol = 1, scales = "free") +
  theme_bw()
plt_clim_ts


### HYDE human density
setwd("/Volumes/xbenitogranell-data/0_project/data/historic/CommonEra")
humanPopdf <- read.csv("HumanDensity_baseline.csv", row.names = 1)
colnames(humanPopdf) <- gsub("X", "", colnames(humanPopdf))
colnames(humanPopdf) <- gsub("_AD", "", colnames(humanPopdf))

select.lakes <- paste(c("Piñan1", "Yahurcch", "Fondocch1", "Llaviucu"), collapse = '|')
test <- humanPopdf %>%
  filter(str_detect(sites, select.lakes))
humanPopLakes <- as.data.frame(t(test))[-1,]
colnames(humanPopLakes) <- c("Yahuarcocha", "Fondococha", "Llaviucu", "Pinan")
humanPopLakes$variable <- "HumanDensity"


## HYDE cropland
cropland <- read.csv("Cropland.csv", row.names=1)
cropland_baseline <- read.csv("Cropland_baseline.csv", row.names=1)

colnames(cropland) <- gsub("X", "", colnames(cropland))
test <- cropland %>%
  filter(str_detect(site, select.lakes))
croplandLakes <- as.data.frame(t(test))[-1,]
colnames(croplandLakes) <- c("Yahuarcocha", "Fondococha", "Llaviucu", "Pinan")
croplandLakes$variable <- "cropland"

#rbind data
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

write.csv(long_matrixHumanHist, "long_matrixHumanHist.csv")

## read in Ecuador SPDs
ecuadorspds <- read.csv("/Volumes/xbenitogranell-data/0_project/data/archaeology/ecuadorspd_ts.csv", row.names = 1)

ecuadorspds$years_AD <- round((ecuadorspds$calBP*(-1)+1950), digits=0) # transform so that time moves in the right direction, and 0 is year 0 of the common era

# Bin SPDs for latter use in CCA
spd_bins1 <- binFunc(as.data.frame(ecuadorspds), as.numeric(ecuadorspds$years_AD), 100, 0, 1700) ##
spd_bins2 <- binFunc(as.data.frame(ecuadorspds), as.numeric(ecuadorspds$years_AD), 30, 1700, 2000) ##

spd_binned <- rbind.data.frame(spd_bins1, spd_bins2)


## bin Human predictors
binhumanData <- function(i, cores, ...) {
  core <- cores[[i]]
  years <- rownames(humanList$HumanDensity) #this is to get ages_AD
  core <- core[ , -which(colnames(core) %in% c("variable"))] 
  
  # Function for binning samples in each core (Alistair Seddon)
  #source("/nfs/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
  source("/volumes/xbenitogranell-data/0_project/R codes/useful scripts/binFunc.R")
  
  #make the age categories
  diatomBin1 <- binFunc(as.data.frame(core), as.numeric(years), 100, 0, 1700) ##
  diatomBin2 <- binFunc(as.data.frame(core), as.numeric(years), 30, 1700, 2000) ##
  
  #diatomBin = binFunc(diatoms, as.numeric(rownames(diatoms)), 40, 20, 2660)
  
  #merge the two binning dataframes for North lakes
  rbind.data.frame(diatomBin1, diatomBin2)
}

#wrap up the function
binnedHuman <- lapply(seq_along(humanList), binhumanData, cores = humanList)
names(binnedHuman) <- names(humanList)

#do Interpolation to fix 1000 yr BP
# HumanDataSitesBinned_interPol <- lapply(seq_along(binnedHuman), doInterpol, cores = binnedHuman)
# names(HumanDataSitesBinned_interPol) <- names(binnedHuman)

#save list
saveRDS(file="/nfs/xbenitogranell-data/0_project/data/historic/CommonEra/binnedHumanFootprint_baseline.rds", binnedHuman)


##
binnedHumanFootprint <- readRDS("/Volumes/xbenitogranell-data/0_project/data/historic/CommonEra/binnedHumanFootprint.rds")
binnedHumanFootprint <- readRDS("/Volumes/xbenitogranell-data/0_project/data/historic/CommonEra/binnedHumanFootprint_baseline.rds")

matrixHumanHist <- plyr::ldply(binnedHumanFootprint, data.frame)

#this is to extract binned ages; only needs to do it once because binned ages are the same for all the cores
ages <- as.data.frame(as.numeric(row.names(binnedHumanFootprint$cropland)))
ages <- plyr::rename(ages,c("as.numeric(row.names(binnedHumanFootprint$cropland))"="age"))

matrixHumanHist <- cbind(matrixHumanHist, ages)

long_matrixHumanHist<-gather(matrixHumanHist, key=lake, value=value, -age, -years, -.id)
#something odd with Llaviucu's cropland in 1890
long_matrixHumanHist[132,]$value <- 2.16
head(long_matrixHumanHist)

# long_matrixHumanHist <- long_matrixHumanHist %>% filter(.id=="cropland") %>%
#   mutate(logvar=sqrt(value))

## plot the data
plt_human_ts <- ggplot(long_matrixHumanHist, aes(x = age, y = value, colour=lake)) +
  geom_line()+
  scale_y_continuous(trans='log10') +
  facet_wrap(~ .id, ncol = 1, scales = "free") +
  theme_bw()
plt_human_ts


### core proxy binned data (ages AD)
load("/Volumes/xbenitogranell-data/0_project/data/training set/binnedProxies.RData") #this is proxy and diatom data for each lake

binnedProxies$umayo <- NULL

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  #core[is.na(core)] <- 0
  core <- core[, colSums(core != 0, na.rm = TRUE) > 0]
  return(core)
}
cores_binned <- lapply(seq_along(binnedProxies), remove, cores=binnedProxies)
#name list elements
names(cores_binned) <- c("Llaviucu", "Yahuarcocha", "Pinan", "Fondococha")

#extract 
nams <- names(cores_binned)
for (i in seq_along(cores_binned)) {
  assign(paste0("", nams[i]), cores_binned[[i]])
}

proxies <- c("d18O", "d13C", "C_N", "Si", "Ca", "K", "Fe_Mn", "Si_Ti", "Ti")

#diatoms
Yahuarcocha_diat <- Yahuarcocha[, -which(colnames(Yahuarcocha) %in% proxies)]
Llaviucu_diat <- Llaviucu[,-which(colnames(Llaviucu) %in% proxies)]
Pinan_diat <- Pinan[,-which(colnames(Pinan) %in% proxies)]
Fondococha_diat <- Fondococha[,-which(colnames(Fondococha) %in% proxies)]


#proxies
  #Llaviucu
  Llaviucu_proxy <- Llaviucu[, which(colnames(Llaviucu) %in% proxies)]
  Llaviucu_proxy <- Llaviucu_proxy[,-which(colnames(Llaviucu_proxy) %in% c("K", "Ca", "Ti", "Si"))]
  Llaviucu_pred <- cbind(Llaviucu_proxy, binnedHumanFootprint$cropland$Llaviucu, binnedHumanFootprint$HumanDensity$Llaviucu, binnedClimateSimulated$MAP$`EpNGEO-J_Llaviucu`,
                            binnedClimateSimulated$MAT$`EpNGEO-J_Llaviucu`, binnedClimateSimulated$TempSeason$`EpNGEO-J_Llaviucu`, binnedClimateSimulated$PrecipSeason$`EpNGEO-J_Llaviucu`,
                         spd_binned$PrDens)
  colnames(Llaviucu_pred)[3:9] <- c("Cropland", "HumanDensity", "MAP", "MAT", "TempSeason", "PrecipSeason", "SPD")
  
  Llaviucu_pred <- transform(Llaviucu_pred, Cropland=log10(Cropland+0.25))
  vifstep(Llaviucu_pred, th=5)
  
  #remove Cimatic data
  Llaviucu_pred$MAT <- NULL
  Llaviucu_pred$MAP <- NULL
  Llaviucu_pred$TempSeason <- NULL
  Llaviucu_pred$PrecipSeason <- NULL
  
  pairs(Llaviucu_pred, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
  
  write.csv(Llaviucu_pred, "Llaviucu_pred_historic.csv")
  
  #Yahuarcocha
  Yahuarcocha_proxy <- Yahuarcocha[, which(colnames(Yahuarcocha) %in% proxies)]
  Yahuarcocha_pred <- cbind(Yahuarcocha_proxy, binnedHumanFootprint$cropland$Yahuarcocha, binnedHumanFootprint$HumanDensity$Yahuarcocha, binnedClimateSimulated$MAP$`EpNGEO-F_Yahurcch`,
                            binnedClimateSimulated$MAT$`EpNGEO-F_Yahurcch`, binnedClimateSimulated$TempSeason$`EpNGEO-F_Yahurcch`, binnedClimateSimulated$PrecipSeason$`EpNGEO-F_Yahurcch`,
                            spd_binned$PrDens)
  colnames(Yahuarcocha_pred)[4:10] <- c("Cropland", "HumanDensity", "MAP", "MAT", "TempSeason", "PrecipSeason", "SPD")
  
  vifstep(Yahuarcocha_pred, th=5)
  
  Yahuarcocha_pred$TempSeason <- NULL
  Yahuarcocha_pred$MAT <- NULL
  Yahuarcocha_pred$MAP <- NULL
  Yahuarcocha_pred$PrecipSeason <- NULL
  
  pairs(Yahuarcocha_pred, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
  Yahuarcocha_pred <- transform(Yahuarcocha_pred, Cropland=log10(Cropland+0.25),
                             HumanDensity=log10(HumanDensity+0.25))
  
  write.csv(Yahuarcocha_pred, "Yahuarcocha_pred_historic.csv")
  
  
  #Pinan
  Pinan_proxy <- Pinan[,which(colnames(Pinan) %in% proxies)]
  Pinan_pred <- cbind(Pinan_proxy, binnedHumanFootprint$cropland$Pinan, binnedHumanFootprint$HumanDensity$Pinan,
                      binnedClimateSimulated$MAP$`EpNGEO-J_Piñan1`, binnedClimateSimulated$MAT$`EpNGEO-J_Piñan1`, binnedClimateSimulated$TempSeason$`EpNGEO-J_Piñan1`, 
                      binnedClimateSimulated$PrecipSeason$`EpNGEO-J_Piñan1`,
                      spd_binned$PrDens)
  colnames(Pinan_pred)[3:9] <- c("Cropland", "HumanDensity", "MAP", "MAT", "TempSeason", "PrecipSeason", "SPD")
  
  vifstep(Pinan_pred, th=5)
  
  Pinan_pred$TempSeason <- NULL
  Pinan_pred$PrecipSeason <- NULL
  Pinan_pred$MAT <- NULL
  Pinan_pred$MAP <- NULL
  
  pairs(Pinan_pred, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
  
  Pinan_pred <- transform(Pinan_pred, C_N=log10(C_N+0.25),
                                HumanDensity=log10(HumanDensity+0.25))
                                
  
  write.csv(Pinan_pred, "Pinan_pred_historic.csv")
  
  
  #Fondococha
  Fondococha_proxy <- Fondococha[,which(colnames(Fondococha) %in% proxies)]
  Fondococha_pred <- cbind(Fondococha_proxy, binnedHumanFootprint$cropland$Fondococha, binnedHumanFootprint$HumanDensity$Fondococha,
                           binnedClimateSimulated$MAP$`EpNGEO-J_Fondocch1`, binnedClimateSimulated$MAT$`EpNGEO-J_Fondocch1`, 
                           binnedClimateSimulated$TempSeason$`EpNGEO-J_Fondocch1`, binnedClimateSimulated$PrecipSeason$`EpNGEO-J_Fondocch1`,
                           spd_binned$PrDens)
  colnames(Fondococha_pred)[3:9] <- c("Cropland", "HumanDensity", "MAP", "MAT", "TempSeason", "PrecipSeason", "SPD")
  
  vifstep(Fondococha_pred, th=5)
  
  Fondococha_pred$TempSeason <- NULL
  Fondococha_pred$MAT <- NULL
  Fondococha_pred$MAP <- NULL
  Fondococha_pred$PrecipSeason <- NULL
  
  pairs(Fondococha_pred, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 
  
  Fondococha_pred <- transform(Fondococha_pred, 
                          HumanDensity=log10(HumanDensity+0.25),
                          Cropland=log10(Cropland+0.25))
  
  write.csv(Fondococha_pred, "Fondococha_pred_historic.csv")
  
  
# Perform CCA for each lake
axis.expl <- function(mod, axes = 1:2) {
    if(is.null(mod$CCA)) {
      sapply(axes, function(i) {
        100*mod$CA$eig[i]/mod$tot.chi
      })
    } else {
      sapply(axes, function(i) {
        100*mod$CCA$eig[i]/mod$tot.chi
      })
    }
  }
  
mod1 <- cca(decostand(Llaviucu_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Llaviucu_pred,
              subset = complete.cases(Llaviucu_diat), scale=TRUE)
(labs_llav <- axis.expl(mod1))
anova.cca(mod1, by="terms")

mod2 <- cca(decostand(Yahuarcocha_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Yahuarcocha_pred,
           subset = complete.cases(Yahuarcocha_diat), scale=TRUE)
(labs_yah <- axis.expl(mod2))
anova.cca(mod2, by="terms")

mod3 <- cca(decostand(Pinan_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Pinan_pred,
            subset = complete.cases(Pinan_diat), scale=TRUE)
(labs_pin <- axis.expl(mod3))
anova.cca(mod3, by="terms")

mod4 <- cca(decostand(Fondococha_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Fondococha_pred,
            subset = complete.cases(Fondococha_diat), scale=TRUE)
(labs_fondo <- axis.expl(mod4))
anova.cca(mod4, by="terms")


# plot CCAs
par(mfrow=c(2,2))
par(mar=c(3,2,2,3))

plot(mod4, display=c('bp', 'sites'), scaling=3, main="Fondococha")
plot(mod4, display=c('bp', 'species'), scaling=3, main="Fondococha")

plot(mod1, display=c('bp', 'sites'), scaling=3, main="Llaviucu")
plot(mod1, display=c('bp', 'species'), scaling=3, main="Llaviucu")

plot(mod3, display=c('bp', 'sites'), scaling=3, main="Piñan")
plot(mod3, display=c('bp', 'species'), scaling=3, main="Piñan")

plot(mod2, display=c('bp', 'sites'), scaling=3, main="Yahuarcocha")
plot(mod2, display=c('bp', 'species'), scaling=3, main="Yahuarcocha")

## Fortify the ordinations for ggploting--for each model
ford <- fortify(mod3, axes = 1:2)  # fortify the ordination
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores
species <- subset(ford, Score == 'species')  # take species scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

sign_llav<- paste(c("Cropland", "HumanDensity"), collapse = '|')
sign_yahu <- paste(c("Cropland", "d13C", "C_N", "d18O", "SPD"), collapse = '|')
sign_pin <- paste(c("C_N", "HumanDensity", "d13C"), collapse = '|')
sign_fondo <- paste(c("Cropland", "Fe_Mn"), collapse = '|')

#plot
cca_llaviucu <- ggplot() +
  theme_bw()+
  # geom_point(data = subset(ford, Score == 'sites'),
  #            mapping = aes(x = CCA1, y = CCA2)) + 
  geom_text(data = subset(ford, Score == 'sites'),
             mapping = aes(label=Label, x = CCA1, y = CCA2),size=3) + 
  geom_segment(data=arrows %>% filter(str_detect(Label, sign_llav)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="blue")+
  geom_segment(data=arrows %>% filter(!str_detect(Label, sign_llav)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="grey")+
  geom_text(data = arrows %>% filter(str_detect(Label, sign_llav)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="blue") +
  geom_text(data = arrows %>% filter(!str_detect(Label, sign_llav)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="grey") +
  xlab(paste0(names(labs_llav[1]), " (", sprintf("%.1f", labs_llav[1]), "%)"))+
  ylab(paste0(names(labs_llav[2]), " (", sprintf("%.1f", labs_llav[2]), "%)"))+
  coord_fixed()

cca_yahuarcocha <- ggplot() +
  theme_bw()+
  # geom_point(data = subset(ford, Score == 'sites'),
  #            mapping = aes(x = CCA1, y = CCA2)) + 
  geom_text(data = subset(ford, Score == 'sites'),
            mapping = aes(label=Label, x = CCA1, y = CCA2),size=3) + 
  geom_segment(data=arrows %>% filter(str_detect(Label, sign_yahu)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="blue")+
  geom_segment(data=arrows %>% filter(!str_detect(Label, sign_yahu)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="grey")+
  geom_text(data = arrows %>% filter(str_detect(Label, sign_yahu)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="blue") +
  geom_text(data = arrows %>% filter(!str_detect(Label, sign_yahu)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="grey") +
  xlab(paste0(names(labs_yah[1]), " (", sprintf("%.1f", labs_yah[1]), "%)"))+
  ylab(paste0(names(labs_yah[2]), " (", sprintf("%.1f", labs_yah[2]), "%)"))+
  coord_fixed()

cca_pinan <- ggplot() +
  theme_bw()+
  # geom_point(data = subset(ford, Score == 'sites'),
  #            mapping = aes(x = CCA1, y = CCA2)) + 
  geom_text(data = subset(ford, Score == 'sites'),
            mapping = aes(label=Label, x = CCA1, y = CCA2),size=3) + 
  geom_segment(data=arrows %>% filter(str_detect(Label, sign_pin)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="blue")+
  geom_segment(data=arrows %>% filter(!str_detect(Label, sign_pin)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="grey")+
  geom_text(data = arrows %>% filter(str_detect(Label, sign_pin)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="blue") +
  geom_text(data = arrows %>% filter(!str_detect(Label, sign_pin)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="grey") +
  xlab(paste0(names(labs_pin[1]), " (", sprintf("%.1f", labs_pin[1]), "%)"))+
  ylab(paste0(names(labs_pin[2]), " (", sprintf("%.1f", labs_pin[2]), "%)"))+
  coord_fixed()

cca_fondo <- ggplot() +
  theme_bw()+
  # geom_point(data = subset(ford, Score == 'sites'),
  #            mapping = aes(x = CCA1, y = CCA2)) + 
  geom_text(data = subset(ford, Score == 'sites'),
            mapping = aes(label=Label, x = CCA1, y = CCA2),size=3) + 
  geom_segment(data=arrows %>% filter(str_detect(Label, sign_fondo)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="blue")+
  geom_segment(data=arrows %>% filter(!str_detect(Label, sign_fondo)),
               mapping = aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
               arrow = arrow(length = unit(0.01, "npc")), colour="grey")+
  geom_text(data = arrows %>% filter(str_detect(Label, sign_fondo)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="blue") +
  geom_text(data = arrows %>% filter(!str_detect(Label, sign_fondo)), # crudely push labels away arrow heads
            mapping = aes(label = Label, x = CCA1 * 1.1, y = CCA2 * 1.1), colour="grey") +
  xlab(paste0(names(labs_fondo[1]), " (", sprintf("%.1f", labs_fondo[1]), "%)"))+
  ylab(paste0(names(labs_fondo[2]), " (", sprintf("%.1f", labs_fondo[2]), "%)"))+
  coord_fixed()

library(cowplot)
plt <- plot_grid(cca_pinan, cca_yahuarcocha, cca_fondo,
                 cca_llaviucu, align="h", ncol = 2, rel_widths = c(2,2),
                 labels = "auto")
plt

ggsave("cca_historical_new.png", plt, height = 8, width = 10)


# extract eigenvalues
var_eigenv_mod4 <- as.data.frame(mod4[["CCA"]][["biplot"]])
rowSums(var_eigenv_mod4[,1:6]) #first 6 CCA axes

var_eigenv_mod1 <- as.data.frame(mod1[["CCA"]][["biplot"]])
rowSums(var_eigenv_mod1[,1:6]) #first 6 CCA axes

var_eigenv_mod3 <- as.data.frame(mod3[["CCA"]][["biplot"]])
rowSums(var_eigenv_mod3[,1:6]) #first 6 CCA axes

var_eigenv_mod2 <- as.data.frame(mod2[["CCA"]][["biplot"]])
rowSums(var_eigenv_mod2[,1:6]) #first 6 CCA axes

# Extract statistically significant CCA axes ()
anova(mod1, by="axis")
ev <- as.vector(eigenvals(mod1, model = "constrained")) #extract eigenvalues for then broken stick

evplot <- function(ev) {
  # Broken stick model (MacArthur 1957) Author: Francois Gillet, 25 August 2012
  n <- length(ev)
  bsm <- data.frame(j=seq(1:n), p=0)
  bsm$p[1] <- 1/n
  for (i in 2:n) bsm$p[i] <- bsm$p[i-1] + (1/(n + 1 - i))
  bsm$p <- 100*bsm$p/n
  # Plot eigenvalues and % of variation for each axis
  op <- par(mfrow=c(2,1))
  barplot(ev, main="Eigenvalues", col="bisque", las=2)
  abline(h=mean(ev), col="red")
  legend("topright", "Average eigenvalue", lwd=1, col=2, bty="n")
  barplot(t(cbind(100*ev/sum(ev), bsm$p[n:1])), beside=TRUE, 
          main="% variation", col=c("bisque",2), las=2)
  legend("topright", c("% eigenvalue", "Broken stick model"), 
         pch=15, col=c("bisque",2), bty="n")
  par(op)
}

br <- evplot(ev)


llaviucu.axes<-scores(mod1,display=c("sites"), choice=c(1,2,3,4,5, 6), scaling=3)
yahuarcocha.axes <- scores(mod2, display=c("sites"), choice=c(1,2,3,4,5), scaling=3)

#then do Procrustes to test association between ordinations and 
# extract residuals of individual observations (years) (timetrack R scripts)

p.prot<-function(ord1,ord2){
  res1<-residuals(ord1)
  res2<-residuals(ord2)
  pprot<-protest(res1,res2,perm=10000)
  yonx<-residuals(pprot)
  yx<-as.matrix(yonx)
  return(yonx)
}

pro <- procrustes(llaviucu.axes, eig_training6, symmetric = FALSE)
