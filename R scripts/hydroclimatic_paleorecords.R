#---------------------------------------------------------------------------
# Script: Hydroclimatic paleorecords (supplementary figure S7) 
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------

# Load libraries for the functions used
library(tidyverse) #manipulate dataframes
library(cowplot) #plot composite ggplots
library(gtools)
library(sp)
library(usdm)

# Read in PAGES2k temperature anomalies (https://www.nature.com/articles/sdata201788) 
SAtempdata <- read.csv("data/paleoclimate/pages2k-tempdata.csv") %>% 
  mutate(age_calBP=(1950-Year.CE)) %>%
  dplyr::select(-Year.CE)

long_SAtempdata <- gather(SAtempdata, key=method, value=Temp, -age_calBP)

ggplot(long_SAtempdata, aes(x=age_calBP,y=Temp, colour=method, group=method)) + 
  geom_line() +
  ylab("Temperature anomalies degC")

#save South American temperature reconstructions in long format
write.csv(long_SAtempdata, "data/paleoclimate/long_SAtempdata.csv")

# Read in paleoclimate records
ubaque <- read.csv("data/paleoclimate/bird-13C-ubaque.csv") #DOI:10.1177/0959683617721324
palestina <- read.csv("data/paleoclimate/data-palestina.csv") #DOI:10.5194/cp-10-1967-2014
pumacocha <- read.csv("data/paleoclimate/data-pumacocha.csv") #DOI: 10.1073/pnas.1003719108
triumfo <- read.csv("data/paleoclimate/data-triumfo.csv") 

# Transform data for that all have agecalyrBP
ubaque$years_AD <- round((ubaque$bchron_age_model*(-1)+1950), digits=0) # transform so that time moves in the right direction, and 0 is year 0 of the common era
ubaque$age_calBP <- (1950-(ubaque$years_AD)) # transform so that time moves in the right direction, and 0 is year 0 of the common era
palestina$age_calBP <- (1950-(palestina$years_AD))
pumacocha$age_calBP <- round(pumacocha$age_calBP, digits=0)
triumfo$years_AD <- round((triumfo$age_calBP*(-1)+1950), digits=0)

#create  a new variable containing record name
ubaque$core <- c(rep("ubaque", nrow(ubaque)))
ubaque <- ubaque %>% dplyr::select(years_AD, age_calBP, d13C, core)
ubaque$long <- -73.935
ubaque$lat <- 4.499

palestina$core <- c(rep("palestina", nrow(palestina)))
palestina <- palestina %>% dplyr::select(years_AD, age_calBP, d18O, d13C, core)

pumacocha$core <- c(rep("pumacocha", nrow(pumacocha)))
pumacocha <- pumacocha %>% dplyr::select(years_AD, age_calBP, d18O, d13C, core)

triumfo$core <- c(rep("triumfo", nrow(triumfo)))
triumfo <- triumfo %>% dplyr::select(years_AD, age_calBP, d13C, core)

#Merge dataframes
df <- bind_rows(ubaque, palestina, pumacocha, triumfo)

#save dataset
write.csv(df, "data/paleoclimate/SouthAmerica_hydroclimate_records.csv")

# Create list
climateList <- split(df, df$core)

#Plot individual hydroclimate record time series
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


#plot records arranged by latitude
plot_hydroclimate_proxies <- plot_grid(plt_ubaq, plt_triumfo, plt_pales, plt_pumac, ncol = 1, align = "hv")
plot_hydroclimate_proxies

# Save plot (Supplementary Fig S7)
ggsave("figures/hydroclimate_timeseries.png", plot_hydroclimate_proxies, height = 8, width = 10)

