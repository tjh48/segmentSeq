#plotLocusDistribution <- function(loci, weights, adjustLength = TRUE, main) {
#  CHGpos <- (end(loci@coordinates) + start(loci@coordinates)) / 2 + c(0, cumsum(chrlens))[match(as.character(seqnames(loci@coordinates)), seqlevels(loci@coordinates))]

#  if(adjustLength)
#    weights <- weights * width(loci@coordinates)
  
#  weights <- weights / sum(weights, na.rm = TRUE)
  
#  plot(density(CHGpos[!is.na(weights)], weights = weights[!is.na(weights)], adjust = 1e-3, cut = 0), main = main)
#  abline(v = cumsum(chrlens)[-length(chrlens)], col = "red", lty = 2, lwd = 3)
#  abline(v = c(centromeres[,1] + c(0, cumsum(chrlens)[-length(chrlens)]), centromeres[,2] + c(0, cumsum(chrlens)[-length(chrlens)])), lty = 2, lwd = 2, col = "blue")
#}

#plotMethDist <- function(mD, subtract, main ="", ...) {
#  reps <- unique(gsub("\\..*", "", levels(mD@replicates)))
#  if(missing(subtract)) subtract = 0
#  listdist <- lapply(1:length(reps), function(ii) {
#    rep <- reps[ii]
#    z <- do.call("cbind", lapply(grep(rep, mD@replicates), function(jj) as.numeric(mD@Cs[,jj]) / as.numeric(mD@Cs[,jj] + mD@Ts[,jj])))
#    wts <- rowSums(z, na.rm = TRUE) / rowSums(!is.na(z))
#    plotMethDistribution(mD, weights = wts, main = "", col = rainbow(length(reps))[ii], add = (ii != 1), lty = ii, lwd = 3, subtract, ...)
#  })
#  legend(col = rainbow(length(reps)), legend = gsub("_", "/", reps), lty = 1:length(reps), lwd = 3, x = "topright")
#  listdist
#}

plotMethDistribution <- function(meth, samples, subtract, centromeres, main = "", add = FALSE, ...) {
  weights = rowSums(meth@Cs[,samples,drop = FALSE]) / rowSums(meth@Cs[,samples, drop = FALSE] + meth@Ts[,samples, drop = FALSE])
  chrlens <- seqlengths(meth@alignments)
  CHGpos <- (end(meth@alignments) + start(meth@alignments)) / 2 + c(0, cumsum(chrlens))[match(as.character(seqnames(meth@alignments)), seqlevels(meth@alignments))]
  
  if(missing(weights)) weights <- 1
  
  weights <- weights
  
  dens <- density(CHGpos[!is.na(weights)], weights = weights[!is.na(weights)] / sum(weights, na.rm = TRUE), adjust = 1e-3, cut = 0)
  dens$y <- dens$y * diff(dens$x)[1]
  if(!missing(subtract))
    dens$y <- dens$y - subtract

  ytext <- min(dens$y) - diff(range(dens$y)) * 0.1  
  if(!add) {
    plot(dens, axes = FALSE,
         ylim = c(min(dens$y) - diff(range(dens$y)) * 0.15, max(dens$y) + diff(range(dens$y)) * 0.1), main = main, ...)
    axis(2)
    abline(v = cumsum(chrlens)[-length(chrlens)], col = "red", lty = 2, lwd = 3)    
    if(!missing(centromeres)) abline(v = c(centromeres[,1] + c(0, cumsum(chrlens)[-length(chrlens)]), centromeres[,2] + c(0, cumsum(chrlens)[-length(chrlens)])), lty = 2, lwd = 2, col = "blue")
    text(seqlevels(meth@alignments), srt = 30, adj = 1, y = ytext, x = cumsum(c(0, chrlens[-length(chrlens)])) + chrlens / 2)
  } else {
    lines(dens, ...)
  }
  dens
}
                         
