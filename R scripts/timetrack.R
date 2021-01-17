#clear workspace
rm(list=ls(all=TRUE))
dev.off()
#unload all loaded packages
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

##loading libraries for functions used
library(analogue) 
library(tidyverse) 
library(rioja) 
library(mgcv)
library(cluster)
library(ggplot2)
library(viridis)
source("R scripts/functions/functions.R")

## This chunk is to process absolute counts core data
#read diatom core datasets
mergedCores <- read.csv("data/mergedCores_diatomcounts.csv")[,-1] 
agedepth <- mergedCores[, names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat <- mergedCores[, !names(mergedCores) %in% c("depth", "upper_age", "lower_age", "lake")]
diat[is.na(diat)] <- 0

diatoms_save <- cbind(agedepth, diat)
coresList <- split(diatoms_save, diatoms_save$lake)

# this is function to calculate relative abundance from counts data
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
diat_red <- diat[, n.occur >2 & abund>3] #more than 3% of RA and present in >2 samples

#check rows with NA
row.has.na <- apply(diat_red, 1, function(x){any(is.na(x))})
sum(row.has.na)

diatoms_save <- cbind(agedepth, diat_red)

#read taxonomic harmonisation diatom names list
changes <- read.csv("data/old_new_nms_cores_counts.csv", stringsAsFactors = FALSE)
#column_new1: ecological groups
#column_new2: harmonized taxonomic names

#transform dataframe to tidy format
new <- diatoms_save %>%
  gather(key = taxa, value = count, -depth, -upper_age, -lower_age, -lake) %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes$old, to = changes$new_2)) %>%
  dplyr::group_by(depth, taxa, lake, upper_age, lower_age) %>%
  summarise(count = sum(count)) %>%
  filter(!count == 0) %>% #this is to remove empty samples (rows)
  filter(!upper_age == 0) %>% #this is to remove ages == 0 (triumfo and fondodocha record)
  spread(key = taxa, value = count) %>%
  as.data.frame()

cores_merged <- new 

## split cores by lakes and reassemble
coresList <- split(cores_merged, cores_merged$lake)

###############################################
## This chunk is to process diatom training set
##Read in environmental lake data
environmental_data_lakes <- read.csv("data/environmental_data_lakes.csv") %>%
  mutate(lake_depth_ratio=Lake_area/Depth_avg) %>%
  mutate(lake_catch_ratio=Lake_area/Wshd_area) %>%
  mutate(catch_vol_ratio=Wshd_area/Vol_total)

rownames(environmental_data_lakes) <- environmental_data_lakes$code

##Read in diatom sediment surface data
training <- read.csv("data/diatoms_trainingset.csv", row.names = 1) #with updated diatom taxonomy and selected spp (>3% of RA and present in >2 samples) plus Miriam's Llaviucu slides

#Regions
lake_regions <- read.csv("data/lakeregions.csv", row.names = 1)

##Merge training set and regions datasets
modern_lakes <- merge(training, lake_regions, by="row.names")

#transform dataframe to tidy format
df_thin <- modern_lakes %>%
  gather(key = taxa, value = count, -Row.names, -region)#don't gather region

#import dataframe with old and new names to group
changes_training <- read.csv("data/old_new_nms_trainingset.csv", stringsAsFactors = FALSE)

#spread
new <- df_thin %>%
  mutate(taxa = plyr::mapvalues(taxa, from = changes_training$old, to = changes_training$new_1)) %>%
  group_by(region, Row.names, taxa) %>%
  summarise(count = sum(count)) %>%
  spread(key = taxa, value = count)

levels(new$region)

#this is to reduce training set to northern Andean lakes only
select.regions <- paste(c("Ecuador", "Colombia", "Junin", "Cusco", "eastern"), collapse = '|')

##Do some filtering in the training set (i.e., remove datapoints that have spp with way too much abundance)
new <- new %>%
  #filter(!str_detect(region, remove.regions)) %>%
  #filter(!str_detect(Row.names, remove.lakes)) %>%
  filter(str_detect(region, select.regions)) %>%
  filter(Fragilaria.crotonensis < 40) %>%
  filter(Cyclostephanos.tholiformis < 40) %>%
  filter(Cyclostephanos.andinus < 40) %>%
  as.data.frame()

training <- new[, -which(names(new) %in% c("Row.names", "region"))]

#For surface plotting ordination: Merge diatom training set and environmental data of lakes
row.names(training) <- new[, which(names(new) %in% c("Row.names"))]
env_surf <- merge(training,environmental_data_lakes, by="row.names")

#For extracting spp from trainingset
training2 <- env_surf[,2:235]

#For extracting environmental variables from diatom training set
env_data_lakes <- env_surf[,236:ncol(env_surf)]
row.names(env_data_lakes) <- env_data_lakes$code

# For merging environmental dataset with lake regions
env_data_lakes <- merge(env_data_lakes,lake_regions, by="row.names")
row.names(env_data_lakes) <- env_data_lakes$code
env_data_lakes <- env_data_lakes[,-c(1,2)]

# Write environmnental lake dataset
write.csv(env_data_lakes, "data/lake_env.csv")

#combine cores and training set (join's analogue package function)
df <- analogue::join(coresList[[1]], coresList[[2]], coresList[[3]], coresList[[4]],
                     coresList[[5]], coresList[[6]], coresList[[7]], coresList[[8]], training2, verbose = TRUE)

#check NA in the list and name the list
listnans <- lapply(df, function(x) sum(is.na(x)))
names(df) <- c("Fondococha", "Lagunillas", "Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

#save the core list
saveRDS(df, "data/coresList.rds")

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  # comment the next line when predicting core trajectories in the timetrack analysis
  #core <- core[, colSums(core) > 0] #select only present species
  #core <- tran(core, "hellinger")
  return(core)
}

cores <- lapply(seq_along(df), remove, cores=df)

#name list elements
names(cores) <- c("Fondococha", "Lagunillas", "Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

#extract training set
training <- cores$trainingset

#drop trainingset from the core list
cores$trainingset <- NULL

#check NA in the list
listnans <- lapply(cores, function(x) sum(is.na(x)))

### Canonical Correspondence Analysis
# Subset candidate environmental variables
env_data <- env_data_lakes
variables <- c("pH", "Cond", "Water.T", "TP", "Depth_avg", "Ca", "Mg", "K", "Elevation",
               "MAT", "P.season", "MAP", "T.season", 
               "lake_depth_ratio", "lake_catch_ratio", "catch_vol_ratio",
               "HFP2009")
env_data <- env_data[,variables]

#transform variables to meet assumptions of homogeneity of variances
env_data <- transform(env_data, Water.T=log10(Water.T+0.25), Elevation=sqrt(Elevation),Cond=log10(Cond+0.25), Ca=log10(Ca+0.25), Mg=log10(Mg+0.25), K=log10(K+0.25), TP=log10(TP+0.25),  
                      Depth_avg=log10(Depth_avg+0.25), MAT=log10(MAT+0.25), P.season=log10(P.season+0.25), MAP=log10(MAP+0.25), T.season=log10(T.season+0.25),
                      lake_depth_ratio=log10(lake_depth_ratio+0.25), lake_catch_ratio=log10(lake_catch_ratio+0.25), catch_vol_ratio=log10(catch_vol_ratio+0.25),
                      HFP2009=log10(HFP2009+0.25))

#Run individual CCAs to select significant predictors on diatom data
training <- training[,colSums(training) > 0] 

ccaResult <- list()
for (i in 1:length(env_data)) {
  mod <- cca(training~env_data[,i], na=na.omit, subset = complete.cases(training), scale=TRUE)
  ccaResult$mod[[i]] <- mod
  #plot(ccaResult$mod[[i]], main=colnames(env_data[i]))
  print(anova(ccaResult$mod[[i]]))
}

names(ccaResult$mod) <- colnames(env_data)

#Subset significant variable for multivariate cca
select <- c("pH", "Cond", "Water.T", "TP", "Ca", "Mg", "K", "Elevation",
            "MAT", "P.season", "MAP", "T.season", 
            "HFP2009")

env_data_red <- env_data[,select]

##boxplots environmental data selected for CCA
par(mfrow = c(4, 4))
par(mar = c(2.5, 3.5, 1, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(0, 0, 3, 0))

for (i in 1:length(env_data_red)) {
  boxplot(env_data_red[,i], cex = 0.6, cex.axis = 0.8,
          las = 1, pch = 19, main=colnames(env_data_red[i])) 
}

#Impute missing values with mice package
# CCA env dataset
library(mice)
data <- env_data_red
tempData <- mice(data,m=5,maxit=50,meth='pmm',seed=500)
summary(tempData)
completedData <- complete(tempData,1)
env_data_red <- completedData

# Check collinearity (set 5 as threshold)
library(usdm)
vifstep(env_data_red, th=5)

##2nd CCA constrained to environmental data with core samples added passively
select <- c("pH", "Cond", "Water.T", "TP", "Ca",
            "MAT", "P.season", "MAP", "T.season", 
            "HFP2009")

env_data_red2 <- env_data_red[,select]

#redo CCA without collinear variables
training <- training[,colSums(training) > 0] 

#duplicate training dataframe before hellinger transformation for subset the more abundant species
training_plt <- training

training <- tran(training, method="hellinger") #give better results transforming
mod_training <- cca(training~., data=env_data_red2, scale=TRUE)
plot(mod_training, scaling = 3)

#Plot eigenvalues and percentages of variation of an ordination object
anova(mod_training, by="axis")
ev <- as.vector(eigenvals(mod_training, model = "constrained")) #extract eigenvalues for then broken stick

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

# Function to determine % of variation explained by each CCA axis
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
(labs <- axis.expl(mod_training))
# labs[1] <- c(18.1)
# labs[2] <- c(16.2)


#Make core sample predictions to add passively into CCA
lake <- "Llaviucu"
lakedepth <- "llaviucu"
pred1<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred1 <- as.data.frame(pred1[,1:2])
pred1$lake <- lake
pred1$years <- coresList[[lakedepth]]$upper_age

lake <- "Yahuarcocha"
lakedepth <- "yahuarcocha"
pred2<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred2 <- as.data.frame(pred2[,1:2])
pred2$lake <- lake
pred2$years <- coresList[[lakedepth]]$upper_age

lake <- "Pinan"
lakedepth <- "pinan"
pred3<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred3 <- as.data.frame(pred3[,1:2])
pred3$lake <- lake
pred3$years <- coresList[[lakedepth]]$upper_age

lake <- "Fondococha"
lakedepth <- "fondococha"
pred4<-predict(mod_training, newdata=decostand(cores[[lake]], "hellinger"), type="wa")
pred4 <- as.data.frame(pred4[,1:2])
pred4$lake <- lake
pred4$years <- coresList[[lakedepth]]$upper_age

pred_cores <- rbind(pred1, pred2, pred3, pred4)
pred_cores <- split(pred_cores, pred_cores$lake)

##### Plot for paper (Figure 4)
# set palette colours
seq_palette <- viridis(4)

png("figures/CCA_timetrack.png", width=10, height=8, units="in", res=300)
layout(matrix(1:2, ncol = 2))

# Extract species scores
scrs <- scores(mod_training, display = "species", scaling = 3)
take <- colnames(training_plt) %in% rownames(scrs)
TAXA <- which(colSums(training_plt[,take] > 0) > 15 & (apply(training_plt[,take]^2, 2, max) > 10)) 
plot(mod_training, display="species", scaling=3,type="n", xlab="", ylab="")
title("Species", adj = 0.3, line = 0.2, cex.main=1, font.main=2, xpd=NA)
title(xlab = paste0(names(labs[1]), " (", sprintf("%.1f", labs[1]), "%)"))
title(ylab = paste0(names(labs[2]), " (", sprintf("%.1f", labs[2]), "%)"))

#customize species scores by ecological groups; must be done before abbreviate
#from TAXA
take1 <- unique(changes_training$new_1)
take2 <- take1 %in% names(TAXA)
nmsdiat <- changes_training$new_2[take2] #FAIL most likely due to repeated names; change manually

nmsdiat <- nmsdiat[-43]
nmsdiat[3] <- "tycoplanktonic"
nmsdiat[4] <- "tycoplanktonic"
nmsdiat[5] <- "tycoplanktonic"
nmsdiat[6] <- "epiphytic"
nmsdiat[7] <- "freshwater planktic"
nmsdiat[10] <- "freshwater planktic"
nmsdiat[11] <- "epiphytic"
nmsdiat[12] <- "benthic"
nmsdiat[24] <- "tycoplanktonic"
nmsdiat[26] <- "freshwater planktic"
nmsdiat[27] <- "tycoplanktonic"
nmsdiat[32] <- "epiphytic"
nmsdiat[40] <- "tycoplanktonic"
nmsdiat[41] <- "tycoplanktonic"
nmsdiat[42] <- "tycoplanktonic"

#points(mod_training, display="species", cex = 0.5, scaling = 3, pch=20, select = TAXA)
points(mod_training, display="species", cex = 0.8, scaling = 3, pch=20, select = TAXA,
       col=as.factor(nmsdiat))

# add abbreviated species names
spnames_abb <- abbreviate(names(training_plt), minlength=10) 
colnames(training_plt) <- spnames_abb
rownames(scrs) <- colnames(training_plt)
take <- colnames(training_plt) %in% rownames(scrs)
TAXA <- which(colSums(training_plt[,take] > 0) > 15 & (apply(training_plt[,take]^2, 2, max) > 10)) 

  # #Extract most abundant spp scores 
  # take_scrs <- rownames(scrs) %in% names(TAXA)
  # scrsdf <- data.frame(scrs)
  # scrs_spp_CCA1 <- as.data.frame(scrsdf[,1][take_scrs])
  # scrs_spp_CCA1_2 <- as.data.frame(cbind(scrs_spp_CCA1, scrsdf[,2][take_scrs]))
  # colnames(scrs_spp_CCA1_2) <-c("CCA1", "CCA2")
  # rownames(scrs_spp_CCA1_2) <- names(TAXA)
  # scrs_spp_CCA1_2$guild <- nmsdiat
  # 
  # #Plot ssp sccores
  # png("CCAaxes_speciesscore.png", width=10, height=8, units="in", res=300)
  # par(mfrow = c(1, 2))
  # par(mar = c(5, 4, 2, 2))
  # datCCA1 <- scrs_spp_CCA1_2[order(scrs_spp_CCA1_2$CCA1),] 
  # barplot(datCCA1$CCA1, col=as.factor(datCCA1$guild), xlab="CCA axis 1 species score",
  #         horiz=TRUE)
  # datCCA2 <- scrs_spp_CCA1_2[order(scrs_spp_CCA1_2$CCA2),]
  # barplot(datCCA2$CCA2, col=as.factor(datCCA2$guild), xlab="CCA axis 2 species score",
  #         horiz=TRUE)
  # legend(-1.5, 40, legend = as.character(unique(datCCA2$guild)), bty = "n",
  #        col = as.factor(unique(datCCA2$guild)), pch = 15,
  #        cex=0.8)
  # dev.off()

training_plt <- tran(training_plt, method="hellinger") #give better results transforming
mod_training_plt <- cca(training_plt~., data=env_data_red2, scale=TRUE)
ordipointlabel(mod_training_plt, display = "species", scaling = 3,
               select = TAXA, cex = 0.7, add = TRUE)

legend("bottomleft", legend = as.character(unique(nmsdiat)), bty = "n",
                      col = as.factor(unique(nmsdiat)), pch = 21, pt.bg = as.factor(unique(nmsdiat)),
                      cex=0.7)

## Plot site (lake) scores
#plot(mod_training, display=c('bp', 'sites'), xlab="", ylab="")
plot(mod_training, type="n", xlab="", ylab="")

###
#Calculate species turnover modern lakes
A <- analogue::distance(training, method="bray")

# create an indicator for all diagonals in the matrix
d <- row(A) - col(A)
# use split to group on these values
diagonal<-split(A, d)
spturn<-unlist(diagonal["1"]) # select relevant one (diag = 1)

###
# CCA-species turnover ordisurf
surf <- ordisurf(mod_training$CCA$u[-1,], spturn,
                 method = "REML", select = FALSE, add = TRUE,
                 col = "red", cex=1, lwd.cl = 1)
#summary(surf)
legend(-1.5, -3.1, col="red", lty=1, cex=0.8,
       legend= c("Species turnover"), bty='n')

#calibrate () calls predict.gam
ordResult <- list()
for(i in seq_along(pred_cores)) {
  newdat <- pred_cores[[i]][,1:2]
  pred_spturn <- calibrate(surf, newdata = newdat, se.fit=TRUE)
  ordResult$pred_spturn[[i]] <- as.numeric(pred_spturn$fit)
  ordResult$pred_spturn_error[[i]] <- as.numeric(pred_spturn$se.fit)
  ordResult$fitted.values[[i]] <- fitted(surf)
  RMSE <- function(error) { sqrt(mean(error^2)) }
  ordResult$RMSE[[i]] <-  RMSE(surf$residuals)
}
names(ordResult$pred_spturn) <- names(pred_cores)
names(ordResult$pred_spturn_error) <- names(pred_cores)

#quantile(fitted(surf))

#extract mean and sd species turnover for each core and add into histogram
data_mean_sd <- list()
for(i in 1:length(ordResult[[1]])) {
  data_mean_sd$mean[[i]]<-mean(as.numeric(lapply(ordResult$pred_spturn,"[",n=i)))
  data_mean_sd$sd[[i]]<-sd(as.numeric(lapply(ordResult$pred_spturn,"[",n=i)))
}

hist(fitted(surf), main="Modern species turnover")
for(i in 1:4) abline(v=data_mean_sd$mean[i], col=i+1, lty=2, lwd=2)
cols <- c("red", "mediumblue", "cadetblue1", "green")
legend("topright", col=cols, lty=1, cex=0.8,
       legend= c("Fondococha", "Piñan", "Yahuarcocha", "Llaviucu"), bty='n')

### Add CCA site scores
points(mod_training, display = 'sites', pch=16, cex=0.7, col="lightgrey")
text(mod_training, display = "bp", cex=0.7, col="black")
title("Sites", adj = 0.3, line = 0.2, cex.main=1, font.main=2, xpd=NA)

#Add core trajectories
points(pred_cores$Pinan[,1:2], type="l", col=seq_palette[1])
points(pred_cores$Yahuarcocha[,1:2], type="l", col=seq_palette[2])
points(pred_cores$Fondococha[,1:2], type="l", col=seq_palette[3])
points(pred_cores$Llaviucu[,1:2], type="l", col=seq_palette[4])

#Add top-bottom core samples
points(pred_cores$Pinan[1, ], pch = 24, cex = 1.6, bg=seq_palette[1],
       col = "grey")
points(pred_cores$Yahuarcocha[1, ], pch = 24, cex = 1.6, bg=seq_palette[2],
       col = "grey")
points(pred_cores$Fondococha[1, ], pch = 24, cex = 1.6, bg=seq_palette[3],
       col = "grey")
points(pred_cores$Llaviucu[1, ], pch = 24, cex = 1.6, bg=seq_palette[4],
       col = "grey")

points(pred_cores$Pinan[nrow(pred_cores$Pinan),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[1])
points(pred_cores$Yahuarcocha[nrow(pred_cores$Yahuarcocha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[2])
points(pred_cores$Fondococha[nrow(pred_cores$Fondococha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[3])
points(pred_cores$Llaviucu[nrow(pred_cores$Llaviucu),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[4])

legend(-1.5, -2.6, pch = c(24, 22), bg=c("forestgreen", "forestgreen"),
       col=c("black","black"),  cex = 0.8,
       legend = c("Core top","Core bottom"), bty = "n")


#par(fig = c(0.02, 0.22, 0.72, 0.97), new = TRUE) #top left
par(fig = c(0.42,0.62, 0.72, 0.97), new = TRUE) #mid left for 1:2 species sires CCA plot

par(mar = c(0, 0, 0, 0))
plot(pred_cores$Pinan[,1:2], type = "n", axes = FALSE, main="")
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Pinan[,1:2], type="l", col=seq_palette[1])
points(pred_cores$Pinan[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[1])
points(pred_cores$Pinan[nrow(pred_cores$Pinan),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[1])
text(pred_cores$Pinan[1,], labels=pred_cores$Pinan$years[1], cex = 0.8, offset = 0.5, pos = 3, col="black", 
     bg="white")
text(pred_cores$Pinan[nrow(pred_cores$Pinan),], labels=pred_cores$Pinan$years[nrow(pred_cores$Pinan)], 
     cex = 0.8, offset = 0.5, pos = 3, col="black", bg="white")
title("Piñan", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

## Plot core insets
#par(fig = c(0.48,0.68, 0.72, 0.97), new = TRUE) #top right
par(fig = c(0.78,0.98, 0.72, 0.97), new = TRUE) #top right for 1:2 species sires CCA plot

par(mar = c(0, 0, 0, 0))
plot(pred_cores$Yahuarcocha[,1:2], type = "n", axes = FALSE)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Yahuarcocha[,1:2], type="l", col=seq_palette[2])
points(pred_cores$Yahuarcocha[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[2])
points(pred_cores$Yahuarcocha[nrow(pred_cores$Yahuarcocha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[2])
text(pred_cores$Yahuarcocha[1,], labels=pred_cores$Yahuarcocha$years[1], cex=0.8, offset=0.5, pos=4, col="black")
text(pred_cores$Yahuarcocha[nrow(pred_cores$Yahuarcocha),], labels=pred_cores$Yahuarcocha$years[nrow(pred_cores$Yahuarcocha)],
     cex = 0.8, offset = 0.5, pos = 4, col="black")
title("Yahuarcocha", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

#par(fig = c(0.02, 0.22, 0.02, 0.27), new = TRUE) #bottom left
par(fig = c(0.42,0.62, 0.02, 0.27), new = TRUE) #mid bottom for 1:2 species sires CCA plot
par(mar = c(0, 0, 0, 0))
plot(pred_cores$Fondococha[,1:2], type = "n", axes=FALSE)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Fondococha[,1:2], type="l", col=seq_palette[3])
points(pred_cores$Fondococha[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[3])
points(pred_cores$Fondococha[nrow(pred_cores$Fondococha),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[3])
text(pred_cores$Fondococha[1,], labels=pred_cores$Fondococha$years[1], cex=0.8, offset=0.5, pos=3, col="black")
text(pred_cores$Fondococha[nrow(pred_cores$Fondococha),], labels=pred_cores$Fondococha$years[nrow(pred_cores$Fondococha)],
     cex=0.8, offset = 0.5, pos = 4, col="black")
title("Fondococha", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

#par(fig = c(0.48,0.68, 0.02, 0.27), new = TRUE) #bottom right
par(fig = c(0.78,0.98, 0.02, 0.27), new = TRUE) #bottom right for 1:2 species sires CCA plot
par(mar = c(0, 0, 0, 0))
plot(pred_cores$Llaviucu[,1:2], type = "n", axes = FALSE)
u <- par("usr")
rect(u[1], u[3], u[2], u[4], col = "white")
points(pred_cores$Llaviucu[,1:2], type="l", col=seq_palette[4])
points(pred_cores$Llaviucu[1, ], pch = 24, cex = 1.6,
       col = "grey", bg=seq_palette[4])
points(pred_cores$Llaviucu[nrow(pred_cores$Llaviucu),], pch = 22, cex = 1.6,
       col = "grey", bg=seq_palette[4])
text(pred_cores$Llaviucu[1,], labels=pred_cores$Llaviucu$years[1], cex=0.8, offset=0.5, pos=2, col="black")
text(pred_cores$Llaviucu[nrow(pred_cores$Llaviucu),], labels=pred_cores$Llaviucu$years[nrow(pred_cores$Llaviucu)], 
     cex=0.8, offset=0.5, pos=3, col="black")
title("Llaviucu", adj = 0.5, line = 0.2, cex.main=0.9, font.main=2, xpd=NA)
box()

#save plot
dev.off()




