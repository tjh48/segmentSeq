
summariseLoci <- function(cD, perReplicate = TRUE)
  {
    if(perReplicate) {
      colSums(exp(cD@locLikelihoods))
    } else {
      sum(1 - exp(rowSums(log(1-exp(cD@locLikelihoods)))))
    }
  }

.controlFDR <- function(likes, FDR) {
  selnum <- max(which(cumsum((1 - sort(likes, decreasing = TRUE)) / 1:length(likes)) < FDR))
  if(selnum > 0) sellikes <- sort(order(likes, decreasing = TRUE)[1:selnum]) else sellikes <- integer()
  sellikes
}

.controlFWER <- function(likes, FWER) {
  llsum <- likes
  selnum <- max(which(1 - cumprod(sort(llsum, decreasing = TRUE)) < FWER))
  if(selnum > 0) sellikes <- sort(order(llsum, decreasing = TRUE)[1:selnum]) else sellikes <- integer()
  sellikes
}

selectLoci <- function(cD, likelihood, FDR, FWER, perReplicate = TRUE) {
  if(!missing(likelihood)) {
    selLoc <- which(rowSums(cD@locLikelihoods > log(likelihood)) > 0)
  } else {
    if(!missing(FDR)) {
      controlFunction <- .controlFDR
      controlCrit <- FDR
    } else if(!missing(FWER)) {
      controlFunction <- .controlFWER
      controlCrit <- FWER
    } else stop ("No criterion for locus selection given.")    
    if(perReplicate) {
      selLoc <- sort(do.call("union", lapply(1:ncol(cD@locLikelihoods), function(jj) controlFunction(exp(cD@locLikelihoods[,jj]), controlCrit))))
    } else {
      selLoc <- controlFunction(1 - exp(rowSums(log(1 - exp(cD@locLikelihoods)))), controlCrit)
    }
  }
  if(length(selLoc) == 0) stop("No loci found for the given selection criterion.")
  cD[selLoc,]
}
