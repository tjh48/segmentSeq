.processPosteriors <- function(lociPD, nullPD, aD, lociCutoff = 0.9, nullCutoff = 0.9, getLikes = TRUE, cl)
  {
    if(nrow(lociPD) > 0) selLoci <- lociPD[which(rowSums(lociPD@locLikelihoods >= log(lociCutoff), na.rm = TRUE) > 0),] else selLoci <- lociPD
    if(nrow(nullPD) > 0) selNull <- nullPD[which(rowSums(nullPD@locLikelihoods >= log(nullCutoff), na.rm = TRUE) > 0),] else selNull <- nullPD
      
    message("Checking overlaps...", appendLF = FALSE)
    locAccept <- matrix(sapply(levels(selLoci@replicates), function(rep)
                        {
                          repCol <- which(levels(selLoci@replicates) == rep)
                          accepts <- rep(FALSE, nrow(selLoci))
                          repLoci <- which(selLoci@locLikelihoods[,repCol] >= log(lociCutoff))
                          if(nrow(selNull) > 0 && length(repLoci) > 0) accepts[repLoci] <- !getOverlaps(selLoci@coordinates[repLoci,], selNull@coordinates[which(selNull@locLikelihoods[,repCol] >= log(nullCutoff)),], overlapType = "contains", whichOverlaps = FALSE, cl = cl) else accepts <- rep(TRUE, nrow(selLoci))
                          message(".", appendLF = FALSE)
                          return(accepts)
                        }), nrow = nrow(lociPD))
    message("done.", appendLF = TRUE)

    message("Selecting loci...", appendLF = FALSE)
    
    lociSegs <- selLoci[rowSums(locAccept) > 0,]
    segAccept <- locAccept[rowSums(locAccept) > 0,,drop = FALSE]

    gc()
    
    ordLSD <- NULL
    ordLSD[order(rowSums(segAccept), width(lociSegs@coordinates), decreasing = TRUE)] <- 1:nrow(lociSegs)
                 
    filLSD <- .filterSegments(lociSegs@coordinates, orderOn = ordLSD, decreasing = FALSE)
    filSegs <- lociSegs[filLSD, ]
    filSegs <- filSegs[order(as.character(seqnames(filSegs@coordinates)), start(filSegs@coordinates), end(filSegs@coordinates)),]

    message("done!")
    
    potDiscards <- unique(unlist(lapply(levels(seqnames(nullPD@coordinates)), function(chrom)
                                        {
                                          whChrom <- which(seqnames(nullPD@coordinates) == chrom)
                                          nullChr <- ranges(nullPD@coordinates)[whChrom,]
                                          
                                          whChrom[which(width(nullChr) < 300 & rowSums(nullPD@data[whChrom,]) == 0 &
                                                        (end(nullChr) + 1) %in% c(start(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom], seqlengths(nullPD@coordinates)[levels(seqnames(nullPD@coordinates)) == chrom]) &
                                                        (start(nullChr) - 1) %in% c(end(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom], 0))
                                                  ]
                                        })))
                                        
    if(length(potDiscards) > 0)
      {
        discardNulls <- nullPD[potDiscards,]
        
        
        message("Extending loci...", appendLF = FALSE)
        
        leftExtendSeg <- unlist(lapply(levels(seqnames(filSegs@coordinates)), function(chrom) {
          message(".", appendLF = FALSE)
          which(seqnames(filSegs@coordinates) == chrom)[match(start(discardNulls@coordinates)[as.character(seqnames(discardNulls@coordinates)) == chrom] - 1,
                          end(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom])]      
        }))
        
        rightExtendSeg <- unlist(lapply(levels(seqnames(filSegs@coordinates)), function(chrom) {
          message(".", appendLF = FALSE)
          which(seqnames(filSegs@coordinates) == chrom)[match(end(discardNulls@coordinates)[as.character(seqnames(discardNulls@coordinates)) == chrom] + 1,
                          start(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom])]      
        }))
        
        extensions <- apply(cbind(potDiscards, leftExtendSeg, rightExtendSeg), 1, function(discards) {
          leftExt <- !is.na(discards[2]) & !any(which(nullPD@locLikelihoods[discards[1],] >= log(nullCutoff)) %in% which(filSegs@locLikelihoods[discards[2],] >= log(lociCutoff)))
          rightExt <- !is.na(discards[3]) & !any(which(nullPD@locLikelihoods[discards[1],] >= log(nullCutoff)) %in% which(filSegs@locLikelihoods[discards[3],] >= log(lociCutoff)))
          if(!leftExt & !rightExt) return(c(0, 0))
          if(leftExt & !rightExt) return(c(nullPD@seglens[discards[1],1], 0))
          if(!leftExt & rightExt) return(c(0, nullPD@seglens[discards[1],1]))
          if(leftExt & rightExt) {
            lenweight <- filSegs@seglens[discards[2],1] + filSegs@seglens[discards[3],1]
            extLeft <- round(nullPD@seglens[discards[1],1] * filSegs@seglens[discards[2],1] / lenweight + runif(1, -1e-3, 1e-3) * as.numeric(filSegs@seglens[discards[2],1] == filSegs@seglens[discards[3],1]))
            extRight <- nullPD@seglens[discards[1],1] - extLeft
            message(".", appendLF = FALSE)
            return(c(extLeft, extRight))
            
          }
        })
        
        message("done!")
        
        
        extSegs <- filSegs
        
        if(length(extensions) > 0)
          {
            end(extSegs@coordinates)[leftExtendSeg[!is.na(leftExtendSeg)]] <- end(extSegs@coordinates)[leftExtendSeg[!is.na(leftExtendSeg)]] + extensions[1,!is.na(leftExtendSeg)]
            start(extSegs@coordinates)[rightExtendSeg[!is.na(rightExtendSeg)]] <- start(extSegs@coordinates)[rightExtendSeg[!is.na(rightExtendSeg)]] - extensions[2,!is.na(rightExtendSeg)]
          }
      } else extSegs <- filSegs    
    
    if (getLikes & nrow(extSegs) > 1) likeSegs <- lociLikelihoods(aD = aD, cD = extSegs, cl = cl) else likeSegs <- extSegs
    
    colnames(likeSegs@locLikelihoods) <- levels(likeSegs@replicates)
    
    likeSegs
  }
