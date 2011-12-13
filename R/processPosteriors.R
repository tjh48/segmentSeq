.processPosteriors <- function(lociPD, nullPD, emptyPD, aD, lociCutoff = 0.9, nullCutoff = 0.9, getLikes = TRUE, cl)
  {
    if(nrow(lociPD) > 0) selLoci <- lociPD[which(rowSums(lociPD@locLikelihoods >= log(lociCutoff), na.rm = TRUE) > 0),] else selLoci <- lociPD
    if(nrow(nullPD) > 0) selNull <- nullPD[which(rowSums(nullPD@locLikelihoods >= log(nullCutoff), na.rm = TRUE) > 0),] else selNull <- nullPD

    if(missing(emptyPD))
      emptyPD <- nullPD[rowSums(nullPD@data) == 0,]
      
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
    
    selLoci <- selLoci[rowSums(locAccept) > 0,]
    segAccept <- locAccept[rowSums(locAccept) > 0,,drop = FALSE]

    rm(selNull, locAccept)
    
    gc()
    
    ordLSD <- NULL
    ordLSD[order(rowSums(segAccept), width(selLoci@coordinates), decreasing = TRUE)] <- 1:nrow(selLoci)
                 
    filLSD <- .filterSegments(selLoci@coordinates, orderOn = ordLSD, decreasing = FALSE)
    filSegs <- selLoci[filLSD, ]
    filSegs <- filSegs[order(as.character(seqnames(filSegs@coordinates)), start(filSegs@coordinates), end(filSegs@coordinates)),]

    message("done!")
    
    potDiscards <- unique(unlist(lapply(seqlevels(emptyPD@coordinates), function(chrom)
                                        {
                                          whChrom <- which(seqnames(emptyPD@coordinates) == chrom)
                                          nullChr <- ranges(emptyPD@coordinates)[whChrom,]
                                          
                                          whChrom[which(
                                                        (end(nullChr) + 1) %in% c(start(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom], seqlengths(emptyPD@coordinates)[levels(seqnames(emptyPD@coordinates)) == chrom]) &
                                                        (start(nullChr) - 1) %in% c(end(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom], 0))
                                                  ]
                                        })))
                                        
    if(length(potDiscards) > 0)
      {
        discardNulls <- emptyPD[potDiscards,]
        
        
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
          leftExt <- !is.na(discards[2]) & !any(which(emptyPD@locLikelihoods[discards[1],] >= log(nullCutoff)) %in% which(filSegs@locLikelihoods[discards[2],] >= log(lociCutoff)))
          rightExt <- !is.na(discards[3]) & !any(which(emptyPD@locLikelihoods[discards[1],] >= log(nullCutoff)) %in% which(filSegs@locLikelihoods[discards[3],] >= log(lociCutoff)))
          if(!leftExt & !rightExt) return(c(0, 0))
          if(leftExt & !rightExt) return(c(width(emptyPD@coordinates)[discards[1]], 0))
          if(!leftExt & rightExt) return(c(0, width(emptyPD@coordinates)[discards[1]]))
          if(leftExt & rightExt) {
            lenweight <- width(filSegs@coordinates)[discards[2]] + width(filSegs@coordinates)[discards[3]]
            extLeft <- round(width(emptyPD@coordinates)[discards[1]] * width(filSegs@coordinates)[discards[2]] / lenweight + runif(1, -1e-3, 1e-3) * as.numeric(width(filSegs@coordinates)[discards[2]] == width(filSegs@coordinates)[discards[3]]))
            extRight <- width(emptyPD@coordinates)[discards[1]] - extLeft
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
