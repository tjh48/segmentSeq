.processPosteriors <- function(lociPD, nullPD, aD, lociCutoff = 0.5, nullCutoff = 0.9, getLikes = TRUE, cl)
  {
    selLoci <- lociPD[which(rowSums(lociPD@posteriors >= log(lociCutoff), na.rm = TRUE) > 0),]
    selNull <- nullPD[which(rowSums(nullPD@posteriors >= log(nullCutoff), na.rm = TRUE) > 0),]

    message("Checking overlaps...", appendLF = FALSE)
    locAccept <- sapply(unique(selLoci@replicates), function(rep)
                        {
                          accepts <- rep(FALSE, nrow(selLoci))
                          repLoci <- which(selLoci@posteriors[,rep] >= log(lociCutoff))
                          accepts[repLoci] <- !getOverlaps(selLoci@annotation[repLoci,], selNull@annotation[selNull@posteriors[,rep] >= log(nullCutoff),], overlapType = "contains", whichOverlaps = FALSE, cl = cl)
                          message(".", appendLF = FALSE)
                          return(accepts)
                        })
    message("done.", appendLF = TRUE)
    
    
    lociSegs <- selLoci[rowSums(locAccept) > 0,]
    segAccept <- locAccept[rowSums(locAccept) > 0,,drop = FALSE]
    
    ordLSD <- NULL
    ordLSD[order(rowSums(segAccept), lociSegs@annotation$end - lociSegs@annotation$start, decreasing = TRUE)] <- 1:nrow(lociSegs)
                 
    filLSD <- filterSegments(lociSegs@annotation, orderOn = ordLSD, decreasing = FALSE)
    filSegs <- lociSegs[filLSD, ]
    filSegs <- filSegs[with(filSegs@annotation, order(as.factor(chr), start, end)), ]
    
    potDiscards <- unique(unlist(lapply(unique(nullPD@annotation$chr), function(chrom)
                                        which(nullPD@annotation$chr == chrom)[with(nullPD@annotation[nullPD@annotation$chr == chrom,],
                                             which(end - start + 1 < 300 & nullClass == "empty" &
                                                   (end + 1) %in% c(filSegs@annotation$start[filSegs@annotation$chr == chrom], aD@chrlens[aD@chrs == chrom]) &
                                                   (start - 1) %in% c(filSegs@annotation$end[filSegs@annotation$chr == chrom], 0)))]
                                        )))

    discardNulls <- nullPD[potDiscards,]
    
    leftExtendSeg <- unlist(lapply(unique(filSegs@annotation$chr), function(chrom)
                                which(filSegs@annotation$chr == chrom)[match(discardNulls@annotation$start[discardNulls@annotation$chr == chrom] - 1, filSegs@annotation$end[filSegs@annotation$chr == chrom])]))
    rightExtendSeg <- unlist(lapply(unique(filSegs@annotation$chr), function(chrom)
                                    which(filSegs@annotation$chr == chrom)[match(discardNulls@annotation$end[discardNulls@annotation$chr == chrom] + 1, filSegs@annotation$start[filSegs@annotation$chr == chrom])]))

    extensions <- apply(cbind(potDiscards, leftExtendSeg, rightExtendSeg), 1, function(discards) {
      leftExt <- !is.na(discards[2]) & !any(which(nullPD@posteriors[discards[1],] >= log(nullCutoff)) %in% which(filSegs@posteriors[discards[2],] >= log(lociCutoff)))
      rightExt <- !is.na(discards[3]) & !any(which(nullPD@posteriors[discards[1],] >= log(nullCutoff)) %in% which(filSegs@posteriors[discards[3],] >= log(lociCutoff)))
      if(!leftExt & !rightExt) return(c(0, 0))
      if(leftExt & !rightExt) return(c(nullPD@seglens[discards[1],1], 0))
      if(!leftExt & rightExt) return(c(0, nullPD@seglens[discards[1],1]))
      if(leftExt & rightExt) {
        lenweight <- filSegs@seglens[discards[2],1] + filSegs@seglens[discards[3],1]
        extLeft <- round(nullPD@seglens[discards[1],1] * filSegs@seglens[discards[2],1] / lenweight + runif(1, -1e-3, 1e-3) * as.numeric(filSegs@seglens[discards[2],1] == filSegs@seglens[discards[3],1]))
        extRight <- nullPD@seglens[discards[1],1] - extLeft
        return(c(extLeft, extRight))
      }
    })

    extSegs <- filSegs
    extSegs@annotation$end[leftExtendSeg[!is.na(leftExtendSeg)]] <- extSegs@annotation$end[leftExtendSeg[!is.na(leftExtendSeg)]] + extensions[1,!is.na(leftExtendSeg)]
    extSegs@annotation$start[rightExtendSeg[!is.na(rightExtendSeg)]] <- extSegs@annotation$start[rightExtendSeg[!is.na(rightExtendSeg)]] - extensions[2,!is.na(rightExtendSeg)]
    
    if (getLikes & nrow(extSegs) > 1) likeSegs <- lociLikelihoods(aD = aD, cD = extSegs, cl = cl) else likeSegs <- extSegs
        
    likeSegs
  }
