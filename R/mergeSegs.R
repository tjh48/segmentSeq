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
        strandSegs@data[merges[,1],] <- countMerge$Cs
        strandSegs@pairData[merges[,1],] <- countMerge$Ts
        
        strandSegs <- strandSegs[-unlist(lapply(1:nrow(merges), function(ii) (merges[ii,1] + 1):merges[ii,2])),]
      }
      strandSegs
    }
    return(strandSegs)
  }

  plusSegs <- strandMerge(segs[which(strand(segs@coordinates) == "+"),])
  minusSegs <- strandMerge(segs[which(strand(segs@coordinates) == "-"),])
  nsSegs <- strandMerge(segs[which(strand(segs@coordinates) == "*"),])

  mergesegs <- new("methData",
                   data = rbind(plusSegs@data, minusSegs@data, nsSegs@data),
                   pairData = rbind(plusSegs@pairData, minusSegs@pairData, nsSegs@pairData),
                   replicates = plusSegs@replicates,
                   coordinates = c(plusSegs@coordinates, minusSegs@coordinates, nsSegs@coordinates),
                   seglens = rbind(plusSegs@seglens, minusSegs@seglens, nsSegs@seglens),
                   libsizes = segs@libsizes,
                   pairLibsizes = segs@pairLibsizes,
                   locLikelihoods = rbind(plusSegs@locLikelihoods, minusSegs@locLikelihoods, nsSegs@locLikelihoods))
  
  mergesegs <- mergesegs[order(as.factor(seqnames(mergesegs@coordinates)), start(mergesegs@coordinates), end(mergesegs@coordinates)),]
  mergesegs
}  
