## ADD A DESCRIPTION HERE

binFunc = function(xDF, Ages, binWidth, minBin, maxBin){
  
  coreMat = xDF
  
  coreBin = seq(minBin, maxBin, binWidth)
  sampleBin = vector("list", length(coreBin)-1)
  
  for(i in 1:length(sampleBin)){
    toBin = which(Ages>=coreBin[i] & Ages<coreBin[i+1])
    sampleBin[[i]] = coreMat[toBin,]
  }
  
  coreRes = matrix(NA, ncol=ncol(xDF), nrow=length(sampleBin))
  if(ncol(xDF) > 1) {
    for(i in 1:length(sampleBin)) coreRes[i,] = apply(sampleBin[[i]], 2, function(x)mean(x, na.rm=T))
  } else { 
    for(i in 1:length(sampleBin)) coreRes[i,] = mean(sampleBin[[i]], na.rm=T)
  }
  colnames(coreRes) = colnames(xDF)
  rownames(coreRes) = coreBin[-1]
  
  coreRes <- as.data.frame(coreRes)
  coreRes
  
}


ttplot <- function(obj, choices = 1:2, pch = c(21, 2), cex = 1,
                   col = c("#888a85", "#e8a600"),
                   lwd = 3, surf = NULL, ccol = "#3465a4",
                   clwd = 2, labcex = 1.2, indicate = TRUE,
                   tpch = 21, bpch = 22, pcex = 1.6,
                   pcol = "black", pbg = "forestgreen",
                   legend = TRUE,
                   ...) {
    scrs <- scores(obj$ord, choices = choices, scaling = obj$scaling,
                   display = "sites")
    fits <- fitted(obj)[, choices]
    xlim <- range(scrs[, 1], fits[, 1])
    ylim <- range(scrs[, 2], fits[, 2])
    plt <- plot(obj$ord, choices = choices, scaling = obj$scaling,
                type = "n", display = "sites", ...,
                ylim = ylim, xlim = xlim)
    points(obj$ord, choices = choices, scaling = obj$scaling,
           display = "sites", pch = pch[1], col = col[1],
           cex = cex, bg = col[1])
    if(!is.null(surf)) {
        op <- par(lwd = clwd)
        on.exit(par(op))
        plot(surf, add = TRUE, col = ccol, labcex = labcex)
        par(op)
    }
    lines(fits, pch = pch[2], col = col[2], lwd = lwd)
    if(indicate) {
        points(fits[1, , drop = FALSE], pch = tpch, cex = pcex,
               col = pcol, bg = pbg)
        points(fits[nrow(fits), , drop = FALSE], pch = bpch, cex = pcex,
               col = pcol, bg = pbg)
    }
    if(legend) {
        legend("topright", col = pcol, pt.bg = pbg, pt.cex = pcex,
               pch = c(tpch, bpch), cex = 1.2,
               legend = c("Top","Bottom"), bty = "n")
    }
}

# this is for using in for loop
multittplot <- function(obj, choices = 1:2, pch = c(21, 2), cex = 1,
                   col = c("#888a85", "#e8a600"),
                   lwd = 3, surf = NULL, ccol = "#3465a4",
                   clwd = 2, labcex = 1.2, indicate = TRUE,
                   tpch = 21, bpch = 22, pcex = 1.6,
                   pcol = "black", pbg = "forestgreen",
                   legend = TRUE,
                   ...) {
  scrs <- scores(obj$ord, choices = choices, scaling = obj$scaling,
                 display = "sites")
  fits <- fitted(obj)[, choices]
  xlim <- range(scrs[, 1], fits[, 1])
  ylim <- range(scrs[, 2], fits[, 2])
  plt <- plot(obj$ord, choices = choices, scaling = obj$scaling,
              type = "n", display = "sites", ...,
              ylim = ylim, xlim = xlim)
  points(obj$ord, choices = choices, scaling = obj$scaling,
         display = "sites", pch = pch[1], col = col[1],
         cex = cex, bg = col[1])
  if(!is.null(surf)) {
    op <- par(lwd = clwd)
    on.exit(par(op))
    plot(surf, add = TRUE, col = ccol, labcex = labcex)
    par(op)
  }
  #lines(fits, pch = pch[2], col = col[2], lwd = lwd)
  if(indicate) {
    points(fits[1, , drop = FALSE], pch = tpch, cex = pcex,
           col = pcol, bg = pbg)
    points(fits[nrow(fits), , drop = FALSE], pch = bpch, cex = pcex,
           col = pcol, bg = pbg)
  }
  if(legend) {
    # legend("bottomright", col = pcol, pt.bg = pbg, pt.cex = pcex,
    #        pch = c(tpch, bpch), cex = 1.2,
    #        legend = c("Top","Bottom"), bty = "n")
    # 
    legend("bottomleft", col = c(pcol, pcol, acol), pt.bg = c(pbg, pbg, abg, NA), pt.cex = pcex,
           pch = c(tpch, bpch, apch), lty = c(NA, NA, NA, 1), cex = 0.8,
           legend = c("Top","Bottom", "Breakpoint (cal yr BP)"), bty = "n")
    
    
  }
}


## PCA plot
PCAplot <- function(obj, choices = 1:2, pch = c(21, 2), cex = 1,
                    col = c("grey", "grey"), surf=NULL,
                    ccol = "#3465a4",
                    clwd = 2, labcex = 1.2, indicate = TRUE,
                    tpch = 21, bpch = 22, pcex = 1.6,
                    pcol = "black", pbg = "forestgreen",
                    legend = TRUE,
                    ...) {
  scrs <- scores(obj, choices = choices, scaling = 3,
                 display = "sites")
  xlim <- range(scrs[, 1], scrs[, 1])
  ylim <- range(scrs[, 2], scrs[, 2])
  plt <- plot(obj, choices = choices, scaling = 3,
              type = "n", display = "sites", ...,
              ylim = ylim, xlim = xlim)
  points(obj, choices = choices, scaling = 3,
         display = "sites", pch = pch[1], col = col[1],
         cex = cex, bg = col[1])
  
  lines(scrs, pch = pch[2], col = col[2], lwd = lwd)
  
  if(!is.null(surf)) {
    op <- par(lwd = clwd)
    on.exit(par(op))
    plot(surf, add = TRUE, col = ccol, labcex = labcex)
    par(op)
  }
  
  if(indicate) {
    points(scrs[1, , drop = FALSE], pch = tpch, cex = pcex,
           col = pcol, bg = pbg)
    points(scrs[nrow(scrs), , drop = FALSE], pch = bpch, cex = pcex,
           col = pcol, bg = pbg)
  }
  if(legend) {
    legend("topright", col = pcol, pt.bg = pbg, pt.cex = pcex,
           pch = c(tpch, bpch), cex = 1.2,
           legend = c("Top","Bottom"), bty = "n")
  }
}

