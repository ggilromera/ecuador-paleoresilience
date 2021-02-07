#---------------------------------------------------------------------------
# Script: Canonical Correspondence Analysis on diatom core abundances against geochemical variables
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------

# Load libraries for the functions used
library(dplyr)
library(tidyverse)
library(cowplot)
library(gtools)
library(vegan)
library(ggvegan)
library(usdm)

### Read in core proxy binned data (ages AD)
load("data/binnedProxies.RData") #this is proxy and diatom data for each lake
binnedHuman <- readRDS("data/binnedHuman.rds") #this is historical human footprint predictors binned

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

#extract dataframes
nams <- names(cores_binned)
for (i in seq_along(cores_binned)) {
  assign(paste0("", nams[i]), cores_binned[[i]])
}

# create a vector containing the proxy names to extract
proxies <- c("d18O", "d13C", "C_N", "Si", "Ca", "K", "Fe_Mn", "Si_Ti", "Ti")

# Extract diatoms
Yahuarcocha_diat <- Yahuarcocha[, -which(colnames(Yahuarcocha) %in% proxies)]
Llaviucu_diat <- Llaviucu[,-which(colnames(Llaviucu) %in% proxies)]
Pinan_diat <- Pinan[,-which(colnames(Pinan) %in% proxies)]
Fondococha_diat <- Fondococha[,-which(colnames(Fondococha) %in% proxies)]


#Extract proxies
#Llaviucu
Llaviucu_proxy <- Llaviucu[, which(colnames(Llaviucu) %in% proxies)]
Llaviucu_proxy <- Llaviucu_proxy[,-which(colnames(Llaviucu_proxy) %in% c("K", "Ca", "Ti", "Si"))]
Llaviucu_pred <- cbind(Llaviucu_proxy, binnedHuman$cropland$Llaviucu, binnedHuman$HumanDensity$Llaviucu)
                       
colnames(Llaviucu_pred)[3:4] <- c("Cropland", "HumanDensity")

Llaviucu_pred <- transform(Llaviucu_pred, Cropland=log10(Cropland+0.25))

#check collinearity among predictors
vifstep(Llaviucu_pred, th=5)

#write.csv(Llaviucu_pred, "Llaviucu_pred_historic.csv")

#Yahuarcocha
Yahuarcocha_proxy <- Yahuarcocha[, which(colnames(Yahuarcocha) %in% proxies)]
Yahuarcocha_pred <- cbind(Yahuarcocha_proxy, binnedHuman$cropland$Yahuarcocha, binnedHuman$HumanDensity$Yahuarcocha)
colnames(Yahuarcocha_pred)[4:5] <- c("Cropland", "HumanDensity")

#check collinearity among predictors
vifstep(Yahuarcocha_pred, th=5)

Yahuarcocha_pred <- transform(Yahuarcocha_pred, Cropland=log10(Cropland+0.25),
                              HumanDensity=log10(HumanDensity+0.25))
#write.csv(Yahuarcocha_pred, "Yahuarcocha_pred_historic.csv")

#Pinan
Pinan_proxy <- Pinan[,which(colnames(Pinan) %in% proxies)]
Pinan_pred <- cbind(Pinan_proxy, binnedHuman$cropland$Pinan, binnedHuman$HumanDensity$Pinan)
colnames(Pinan_pred)[3:4] <- c("Cropland", "HumanDensity")

#check collinearity among predictors
vifstep(Pinan_pred, th=5)

Pinan_pred <- transform(Pinan_pred, C_N=log10(C_N+0.25),
                        HumanDensity=log10(HumanDensity+0.25))
#write.csv(Pinan_pred, "Pinan_pred_historic.csv")


#Fondococha
Fondococha_proxy <- Fondococha[,which(colnames(Fondococha) %in% proxies)]
Fondococha_pred <- cbind(Fondococha_proxy, binnedHuman$cropland$Fondococha, binnedHuman$HumanDensity$Fondococha)
colnames(Fondococha_pred)[3:4] <- c("Cropland", "HumanDensity")

#check collinearity among predictors
vifstep(Fondococha_pred, th=5)

pairs(Fondococha_pred, diag.panel = panel.hist, upper.panel = panel.smooth, lower.panel = panel.cor, gap = 0, cex.labels = 1, cex=1.5, font.labels = 1) 

Fondococha_pred <- transform(Fondococha_pred, 
                             HumanDensity=log10(HumanDensity+0.25),
                             Cropland=log10(Cropland+0.25))

#write.csv(Fondococha_pred, "Fondococha_pred_historic.csv")


# Perform CCA for each lake core
# Function to extract % of explained variability for CCA axes 1 and 2
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

mod_Llaviucu <- cca(decostand(Llaviucu_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Llaviucu_pred,
            subset = complete.cases(Llaviucu_diat), scale=TRUE)
(labs_llav <- axis.expl(mod_Llaviucu))
# Test which proxies explain a significant portion of diatom variability
anova.cca(mod_Llaviucu, by="terms")

mod_Yahuarcocha <- cca(decostand(Yahuarcocha_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Yahuarcocha_pred,
            subset = complete.cases(Yahuarcocha_diat), scale=TRUE)
(labs_yah <- axis.expl(mod_Yahuarcocha))
anova.cca(mod_Yahuarcocha, by="terms")

mod_Pinan <- cca(decostand(Pinan_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Pinan_pred,
            subset = complete.cases(Pinan_diat), scale=TRUE)
(labs_pin <- axis.expl(mod_Pinan))
anova.cca(mod_Pinan, by="terms")

mod4_Fondococha <- cca(decostand(Fondococha_diat, "hell", na.rm=TRUE)~., na=na.omit, data=Fondococha_pred,
            subset = complete.cases(Fondococha_diat), scale=TRUE)
(labs_fondo <- axis.expl(mod_Pinan))
anova.cca(mod4_Fondococha, by="terms")


# plot CCAs (base R)
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

### Plot CCAs (ggplot)
# Fortify the ordinations for ggploting--for each model
ford <- fortify(mod_Llaviucu, axes = 1:2)  # fortify the ordination
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores
species <- subset(ford, Score == 'species')  # take species scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

# create a vector with significant variables
sign_llav<- paste(c("Cropland", "HumanDensity"), collapse = '|')

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


# Fortify the ordinations for ggploting--for each model
ford <- fortify(mod_Yahuarcocha, axes = 1:2)  # fortify the ordination
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores
species <- subset(ford, Score == 'species')  # take species scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

# create a vector with significant variables
sign_yahu <- paste(c("Cropland", "d13C", "C_N", "d18O"), collapse = '|')

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

# Fortify the ordinations for ggploting--for each model
ford <- fortify(mod_Pinan, axes = 1:2)  # fortify the ordination
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores
species <- subset(ford, Score == 'species')  # take species scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

# create a vector with significant variables
sign_pin <- paste(c("C_N", "HumanDensity"), collapse = '|')

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

# Fortify the ordinations for ggploting--for each model
ford <- fortify(mod4_Fondococha, axes = 1:2)  # fortify the ordination
take <- c('CCA1', 'CCA2')  # which columns contain the scores we want
arrows <- subset(ford, Score == 'biplot')  # take only biplot arrow scores
species <- subset(ford, Score == 'species')  # take species scores

## multiplier for arrows to scale them to the plot range
mul <- ggvegan:::arrowMul(arrows[, take],
                          subset(ford, select = take, Score == 'sites'))
arrows[, take] <- arrows[, take] * mul  # scale biplot arrows

# create a vector with significant variables
sign_fondo <- paste(c("Cropland", "HumanDensity"), collapse = '|')

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

# Plot composite CCA plot (Fig S5 of the manuscript)
plt_cca <- plot_grid(cca_pinan, cca_yahuarcocha, cca_fondo,
                 cca_llaviucu, align="h", ncol = 2, rel_widths = c(2,2),
                 labels = "auto")
plt_cca

ggsave("figures/cca_historical_new.png", plt, height = 8, width = 10)

