#---------------------------------------------------------------------------
# Script: Hierarchical cluster analysis
# Paper: Ecological resilience of tropical Andean lakes: a paleolimnological perspective
# Author: Benito, X.
# e-mail: xavier.benito.granell@gmail.com
#---------------------------------------------------------------------------

##loading libraries for functions used
library(analogue) #to join diatom datasets on their common spp
library(tidyverse)
library(cluster) 
library(rioja) #to perform constrained hierarchical clustering

## HRead in the diatom core list with absolute counts
df <- readRDS("data/coresList.rds")

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  core <- core[, colSums(core) > 0] #select only present species
  return(core)
}

#Wrap-up function over the list
cores <- lapply(seq_along(df), remove, cores=df)

#name list elements
names(cores) <- c("Fondococha", "Lagunillas", "Llaviucu", "Pinan", "Titicaca", "Triumfo", "Umayo", "Yahuarcocha", "trainingset")

#drop trainingset from the core list
cores$trainingset <- NULL

lakes <- c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha")
cores <- cores[lakes]

## Perform MULTIPLE constrained cluster analyses
doCluster <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, colSums(core) > 0] #select only present species
  core <- decostand(core, method="hellinger") #Hellinger transform relative abundance data
  diss <- vegan::vegdist(core, method="bray")
  clust <- chclust(diss, method="coniss")
  #bstick(clust)
  clust
}

coreCluster <- lapply(seq_along(cores), doCluster, cores = cores)
names(coreCluster) <- lakes

## plot Clusters
par(mfrow=c(2,2))
par(mar=c(4,2,2,2))

for(i in seq_along(coreCluster)) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  plot(coreCluster[[i]], hang=-1)
  title(names(cores[i]))
  par(op)
}

# Estimate broken stick model for significant zones (visual inspection)
for (i in seq_along(coreCluster)) {
  bstick(coreCluster[[i]])
  title(names(cores[i]))
}

#extract core clusters
nams <- names(coreCluster)
for (i in seq_along(coreCluster)) {
  assign(paste0("cCluster_", nams[i]), coreCluster[[i]])
}

#Extract significant diatom zones
# a splits = a+1 zones (k), so split the groups of the ccluster analysis accordingly
sigClustLlav <- cutree(coreCluster$Llaviucu, k=4)
locate <- cumsum(rle(sigClustLlav)$lengths)+1
zones_llaviucu <- coresList$llaviucu$upper_age[locate]

sigClustPin <- cutree(coreCluster$Pinan, k=4)
locate <- cumsum(rle(sigClustPin)$lengths)+1
zones_pinan <- coresList$pinan$upper_age[locate]

sigClustYah <- cutree(coreCluster$Yahuarcocha, k=4)
locate <- cumsum(rle(sigClustYah)$lengths)+1
zones_yahuarcocha <- coresList$yahuarcocha$upper_age[locate]

sigClustFondo <- cutree(coreCluster$Fondococha, k=4)
locate <- cumsum(rle(sigClustFondo)$lengths)+1
zones_fondococha <- coresList$fondococha$upper_age[locate]
