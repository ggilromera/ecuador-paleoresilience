
## THIS IS NO LONGER NEEEDED --> BrayCurtis_timeseries

## Here to load the diatom core list with absolute counts
df <- readRDS("data/coresList.rds")

#Remove empty spp resulting from merging dataframes (and drop year & depths vars)
remove <- function(i, cores, ...) {
  core <- cores[[i]]
  core <- core[, -which(names(core) %in% c("depth", "upper_age", "lower_age", "lake", "AgeCE"))] # drop year & depths vars
  # comment the next line when predicting core trajectories in the timetrack analysis
  core <- core[, colSums(core) > 0] #select only present species
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

#Calculate BC dissimilarity from the reference assemblage
doBCcumulative <- function(i, cores,..) {
  core <- cores[[i]]
  D <- analogue::distance(core/100, method="bray")
  SQD<-D[,ncol(D)] #take the oldest sample as reference to see the direction of change
  return(SQD)
}

lakes <- c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha")

coresBCcum <- lapply(seq_along(cores[lakes]), doBCcumulative, cores=cores[lakes])
names(coresSQD) <- c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha")

goodpoorbad <- list()
for (i in 1:length(coresBCcum)) {
  goodpoorbad[[i]] <- quantile(coresBCcum[[i]], probs = c(0.75, 0.95))
}
names(goodpoorbad) <- c("Fondococha", "Llaviucu", "Pinan", "Yahuarcocha")


##plot BC dissimilarity from the reference conditions
par(mfrow = c(2, 2))
par(mar = c(2.5, 3.5, 1, 0.5))
par(mgp = c(1.5, 0.5, 0))
par(oma = c(0, 0, 3, 0))

for (i in 1:length(coresBCcum)) {
  plot.ts(coresSQD[[i]], cex = 0.6, cex.axis = 0.8,
          las = 1, pch = 19, main=names(coresBCcum[i]), ylim=c(0,1),
          ylab="BC Dissimilarity")
  abline(h=goodpoorbad[[i]][2], col="blue")
}


