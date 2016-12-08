mergeMethSegs <- function(segs, aD, gap, cl) {

  strandMerge <- function(strandSegs) {
    if(nrow(strandSegs) > 0) {
      #strandSegs@locLikelihoods[is.na(strandSegs@locLikelihoods)] <- Inf
      whsim <- which(rowSums(strandSegs@locLikelihoods[-1,,drop = FALSE] == strandSegs@locLikelihoods[-nrow(strandSegs),,drop = FALSE], na.rm = TRUE) == rowSums(!is.na(strandSegs@locLikelihoods[-1,, drop = FALSE]) & !is.na(strandSegs@locLikelihoods[-nrow(strandSegs),,drop =FALSE])) &
                     seqnames(strandSegs@coordinates)[-1] == seqnames(strandSegs@coordinates)[-nrow(strandSegs)] &
                     start(strandSegs@coordinates)[-1] - end(strandSegs@coordinates)[-nrow(strandSegs)] + 1 < gap)
      if(length(whsim)) {    
        strandSegs@locLikelihoods[strandSegs@locLikelihoods == Inf] <- NA
        endgroup <- which(diff(whsim) != 1)
        merges <- cbind(whsim[c(1, endgroup + 1)], whsim[c(endgroup, length(whsim))] + 1)  
        end(strandSegs@coordinates)[merges[,1]] <- end(strandSegs@coordinates)[merges[,2]]
        
        countMerge <- getCounts(strandSegs@coordinates[merges[,1]], aD, cl = cl)
        strandSegs@data[merges[,1],,] <- array(c(countMerge$Cs, countMerge$Ts), c(dim(countMerge$Cs), 2))        
        strandSegs <- strandSegs[-unlist(lapply(1:nrow(merges), function(ii) (merges[ii,1] + 1):merges[ii,2])),]
      }
      strandSegs
    }
    return(strandSegs)
  }

  if(any(strand(segs@coordinates) == "+")) 
      plusSegs <- strandMerge(segs[which(strand(segs@coordinates) == "+"),]) else plusSegs <- segs[0,]

  if(any(strand(segs@coordinates) == "-")) 
      minusSegs <- strandMerge(segs[which(strand(segs@coordinates) == "-"),]) else minusSegs <- segs[0,]
  if(any(strand(segs@coordinates) == "*")) 
      nsSegs <- strandMerge(segs[which(strand(segs@coordinates) == "*"),]) else nsSegs <- segs[0,]

  mergesegs <- .mergeListLoci(list(plusSegs, minusSegs, nsSegs))
  
  mergesegs <- mergesegs[order(as.factor(seqnames(mergesegs@coordinates)), start(mergesegs@coordinates), end(mergesegs@coordinates)),]
  mergesegs
}  
