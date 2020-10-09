#Read core list
readRDS("data/")

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
names(coreCluster) <- names(cores)

## plot Clusters
par(mfrow=c(3,3))
par(mar=c(4,2,2,2))

for(i in seq_along(coreCluster)) {
  op <- par(mar = c(5,4,1,1) + 0.1)
  plot(coreCluster[[i]], hang=-1)
  title(names(coreCluster[i]))
  par(op)
}

# Estimate broken stick model for significant zones (visual inspection)
for (i in seq_along(coreCluster)) {
  core <- coreCluster[[i]]
  n <- ncol(cores[[i]])
  k <- seq(1,n)
  sumFrac <- 1/k
  bstick <- rep(NA,n)
  for(j in 1:n) bstick[j] = 1/n*sum(sumFrac[j:n])
  
  # Compare variance explained by each split
  clustVarEx <- rev(diff(core$height)/max(core$height))
  plot(k[1:10], clustVarEx[1:10], type="o", col="black", xlab="Number of Splits", ylab="% Var Expl")
  title(names(coreCluster[i]))
  points(k[1:10], bstick[1:10], type="o", col="red")
  
}

# Estimate broken stick model for significant zones (visual inspection)
for (i in seq_along(coreCluster)) {
  bstick(coreCluster[[i]])
  title(names(coreCluster[i]))
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

sigClustPin <- cutree(cCluster_Pinan, k=4)
locate <- cumsum(rle(sigClustPin)$lengths)+1
zones_pinan <- coresList$pinan$upper_age[locate]

#sigClustTriumfo <- cutree(cCluster_Triumfo, k=3)
#sigClustUmy <- cutree(cCluster_Umayo, k=3)
sigClustYah <- cutree(cCluster_Yahuarcocha, k=4)
locate <- cumsum(rle(sigClustYah)$lengths)+1
zones_yahuarcocha <- coresList$yahuarcocha$upper_age[locate]

sigClustFondo <- cutree(cCluster_Fondococha, k=4)
locate <- cumsum(rle(sigClustFondo)$lengths)+1
zones_fondococha <- coresList$fondococha$upper_age[locate]
