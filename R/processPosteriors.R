# modification on git from copied files
.extractStrand <- function(sD, strand) {
  sDPlus <- sD[,grep(strand, colnames(sD@Cs))]
  sDPlus@replicates <- as.factor(gsub(paste("\\.", strand, sep = ""), "", as.character(sDPlus@replicates[grep(strand, sDPlus@replicates)])))
  if(strand == "plus") strand(sDPlus@coordinates) = "+" else if(strand == "minus") strand(sDPlus@coordinates) = "-"
  sDPlus@locLikelihoods <- sDPlus@locLikelihoods[,grep(strand, colnames(sDPlus@locLikelihoods))]
  colnames(sDPlus@locLikelihoods) <- gsub(paste("\\.", strand, sep = ""), "", colnames(sDPlus@locLikelihoods))
  colnames(sDPlus@Cs) <- gsub(paste("\\.", strand, sep = ""), "", colnames(sDPlus@Cs))
  colnames(sDPlus@Ts) <- gsub(paste("\\.", strand, sep = ""), "", colnames(sDPlus@Ts))
  sDPlus
}

.processPosteriors <- function(lociPD, nullPD, emptyPD, aD, lociCutoff = 0.9, nullCutoff = 0.9, extendLoci, getLikes = FALSE, verbose = TRUE, cl)
  {
      if(missing(emptyPD)) emptyPD <- NULL
    
      strandSegs <- lapply(levels(strand(lociPD@coordinates)), function(ss)
          if(any(strand(lociPD@coordinates) == ss))
              {
                  if(!is.null(emptyPD)) strandEmpty <- emptyPD[which(strand(emptyPD@coordinates) %in% list("+", "-", c("+", "-", "*"))[[which(c("+", "-", "*") == ss)]]),] else strandEmpty <- NULL
                  if(verbose) message("Strand ", ss)
                  strandSegs <- .processPosts(
                      lociPD = lociPD[which(strand(lociPD@coordinates) == ss),],
                      nullPD = nullPD[which(strand(nullPD@coordinates) %in% list("+", "-", c("+", "-", "*"))[[which(c("+", "-", "*") == ss)]]),],
                      emptyPD = strandEmpty,
                      aD, lociCutoff, nullCutoff, extendLoci = extendLoci, getLikes = FALSE,
                      verbose = verbose, cl = cl)
                  strandSegs
              })
      strandSegs <- strandSegs[which(!sapply(strandSegs, is.null))]

    segs <- .mergeListLoci(strandSegs)
    
#    if(class(lociPD) == "segMeth") {
#      segs <- new("methData",
#                  data = do.call("abind", c(lapply(splitSeg, function(x) x@data), along = 1))
#                  data = do.call("abind", c(strandSegs, list(along = 1))),
#                  replicates = lociPD@replicates,
#                  coordinates = do.call("c", lapply(strandSegs, function(x) if(!is.null(x) && length(x@coordinates) > 0) ret#urn(x@coordinates) else return(GRanges(seqinfo = seqinfo(x@coordinates))))),                  
#                  locLikelihoods = do.call("rbind", lapply(strandSegs, function(x) if(!is.null(x)) x@locLikelihoods)),
#                  sampleObservables = list(nonconversion = lociPD@nonconversion))
#      } else {
#        segs <- new("lociData",
#                    data = do.call("rbind", lapply(strandSegs, function(x) if(!is.null(x)) x@data)),
#                    replicates = lociPD@replicates,
#                    coordinates = do.call("c", lapply(strandSegs, function(x) if(!is.null(x)) return(x@coordinates) else return(GRanges()))),
#                    seglens = do.call("c", lapply(strandSegs, function(x) if(!is.null(x)) x@rowObservables$seglens)),
#                    locLikelihoods = do.call("rbind", lapply(strandSegs, function(x) if(!is.null(x)) x@locLikelihoods)))
#      }
    
    
    segs <- segs[order(as.factor(seqnames(segs@coordinates)), start(segs@coordinates), end(segs@coordinates)),]      
    segs
  }

.processPosts <- function(lociPD, nullPD, emptyPD, aD, lociCutoff = 0.9, nullCutoff = 0.9, getLikes = FALSE, verbose = TRUE, extendLoci, cl)
  {
      if(!is.null(cl))
          clusterEvalQ(cl, rm(list = ls()))

      if(nrow(lociPD) > 0) selLoci <- lociPD[which(rowSums(lociPD@locLikelihoods >= log(lociCutoff), na.rm = TRUE) > 0),] else selLoci <- lociPD
      if(nrow(nullPD) > 0) selNull <- nullPD[which(rowSums(nullPD@locLikelihoods >= log(nullCutoff), na.rm = TRUE) > 0),] else selNull <- nullPD
      
#    if(missing(emptyPD) & class(lociPD) == "segData")
#      emptyPD <- nullPD[rowSums(sapply(1:ncol(nullPD), function(jj) as.integer(nullPD@data[,jj]))) == 0,]

    if(verbose) message("Checking overlaps...", appendLF = FALSE)
#      locAccept <- do.call("cbind", 
#                           lapply(levels(selLoci@replicates), function(rep)
      for(rep in levels(selLoci@replicates))
          {
              if(verbose) message(".", appendLF = FALSE)
              repCol <- which(levels(selLoci@replicates) == rep)
              accepts <- rep(FALSE, nrow(selLoci))
              repLoci <- which(selLoci@locLikelihoods[,repCol] >= log(lociCutoff))
              accepts[repLoci] <- TRUE
              if(nrow(selNull) > 0 && length(repLoci) > 0)
                  {
                      repNull <- sort(selNull@coordinates[which(selNull@locLikelihoods[,repCol] >= log(nullCutoff)),])
                      if(length(repNull) > 1) repNull <- repNull[c(1, which(diff(start(repNull)) != 0))+1,]
                      accepts[repLoci] <- !.getOverlaps(selLoci@coordinates[repLoci,], repNull, overlapType = "contains", whichOverlaps = FALSE, cl = NULL)
                  }

#              if(nrow(selNull) > 0 && length(repLoci) > 0) {
#                  repNull <- sort(selNull@coordinates[which(selNull@locLikelihoods[,repCol] >= log(nullCutoff)),])
#                  if(length(repNull) > 1) repNull <- repNull[c(1, which(diff(start(repNull)) != 0))+1,]
#                  
#                  repNull <- repNull[c(which(diff(end(repNull)) != 0), length(repNull)),]
#                  whichCN <- sort(unique(subjectHits(findOverlaps(repNull, selLoci@coordinates[repLoci,], type = "within"))))                                     
#                  accepts[repLoci[whichCN]] <- FALSE
#              }
              selLoci@locLikelihoods[,which(levels(replicates(selLoci)) == rep)] <- log(as.integer(accepts))
          }
      if(verbose) message("done.", appendLF = TRUE)
    
    if(verbose) message("Selecting loci...", appendLF = FALSE)

    locTrue <- unique(which(selLoci@locLikelihoods == 0, arr.ind = TRUE)[,1])
    selLoci <- selLoci[locTrue,]
      segAccept <- selLoci@locLikelihoods == 0
    
    rm(selNull, locTrue)
    
    gc()


    filterOnNumberLength <- function(segAccept, selLoci)
      {
        selTrue <- sort(unique(unlist(lapply(1:ncol(segAccept), function(rep) which(segAccept[,rep])))))
        if(length(selTrue) == 0) selTrue <- 0
        segAccept <- segAccept[selTrue,,drop = FALSE]
        selLoci <- selLoci[selTrue,,drop = FALSE]
        
        ordLSD <- NULL
        ordLSD[order(rowSums(do.call("cbind", lapply(1:ncol(segAccept), function(rep) as.integer(segAccept[,rep])))), width(selLoci@coordinates), decreasing = TRUE)] <- 1:nrow(selLoci)
        
        filLSD <- .filterSegments(selLoci@coordinates, orderOn = ordLSD, decreasing = FALSE)
        filSegs <- selLoci[filLSD,]
        filSegs <- filSegs[order(as.factor(seqnames(filSegs@coordinates)), start(filSegs@coordinates), end(filSegs@coordinates)),]
        filSegs
      }

    filSegs <- filterOnNumberLength(segAccept, selLoci)

    if(verbose) message("done!")

    if(extendLoci)
      {
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
            if(verbose) message("Extending loci...", appendLF = FALSE)
            
            leftExtendSeg <- unlist(lapply(levels(seqnames(filSegs@coordinates)), function(chrom) {
              if(verbose) message(".", appendLF = FALSE)
              which(seqnames(filSegs@coordinates) == chrom)[match(start(discardNulls@coordinates)[as.character(seqnames(discardNulls@coordinates)) == chrom] - 1,
                              end(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom])]      
            }))
            
            rightExtendSeg <- unlist(lapply(levels(seqnames(filSegs@coordinates)), function(chrom) {
              if(verbose) message(".", appendLF = FALSE)
              which(seqnames(filSegs@coordinates) == chrom)[match(end(discardNulls@coordinates)[as.character(seqnames(discardNulls@coordinates)) == chrom] + 1,
                              start(filSegs@coordinates)[as.character(seqnames(filSegs@coordinates)) == chrom])]      
            }))        
            
                                        #        return(list(rightExtendSeg, leftExtendSeg, potDiscards, emptyPD))
            
            leftLikes <- matrix(NA, nrow = length(leftExtendSeg), ncol = ncol(filSegs@locLikelihoods))
            leftLikes[!is.na(leftExtendSeg),] <- sapply(1:ncol(filSegs@locLikelihoods), function(ii) as.double(filSegs@locLikelihoods[leftExtendSeg[!is.na(leftExtendSeg)],ii]) >= log(lociCutoff))
            
            rightLikes <- matrix(NA, nrow = length(rightExtendSeg), ncol = ncol(filSegs@locLikelihoods))
            rightLikes[!is.na(rightExtendSeg),] <- sapply(1:ncol(filSegs@locLikelihoods), function(ii) as.double(filSegs@locLikelihoods[rightExtendSeg[!is.na(rightExtendSeg)],ii]) >= log(lociCutoff))
            
            emptyLoc <- sapply(1:ncol(emptyPD@locLikelihoods), function(ii) as.double(emptyPD@locLikelihoods[potDiscards,ii]) >= log(nullCutoff))
            
            leftExt <- (!is.na(leftExtendSeg) & rowSums(emptyLoc & leftLikes, na.rm = TRUE) == 0)
            rightExt <- (!is.na(rightExtendSeg) & rowSums(emptyLoc & rightLikes, na.rm = TRUE) == 0)
            
            extensions <- matrix(NA, ncol = 2, nrow = length(potDiscards))
            extensions[!leftExt & !rightExt,] <- c(0, 0)
            extensions[leftExt & !rightExt,] <- cbind(width(emptyPD@coordinates)[potDiscards[leftExt & !rightExt]], 0)
            extensions[!leftExt & rightExt,] <- cbind(0, width(emptyPD@coordinates)[potDiscards[!leftExt & rightExt]])
            
            lenweights <- width(filSegs@coordinates)[leftExtendSeg] + width(filSegs@coordinates)[rightExtendSeg]
            extLeft <- round(width(emptyPD@coordinates)[potDiscards] * width(filSegs@coordinates)[leftExtendSeg] / lenweights + runif(1, -1e-3, 1e-3) * as.numeric(width(filSegs@coordinates)[leftExtendSeg] == width(filSegs@coordinates)[rightExtendSeg]))
            extRight <- width(emptyPD@coordinates)[potDiscards] - extLeft
            extensions[leftExt & rightExt] <- cbind(extLeft, extRight)[leftExt & rightExt,]
            
                                        #            extensions <- t(apply(cbind(potDiscards, leftExtendSeg, rightExtendSeg), 1, function(discards) {
                                        #              
                                        #              leftExt <- !is.na(discards[2]) && !any(which(unlist(emptyPD@locLikelihoods[discards[1],,drop = TRUE]) >= log(nullCutoff)) %in% which(unlist(filSegs@locLikelihoods[discards[2],,drop = TRUE]) >= log(lociCutoff)))
                                        #              rightExt <- !is.na(discards[3]) && !any(which(unlist(emptyPD@locLikelihoods[discards[1],,drop = TRUE]) >= log(nullCutoff)) %in% which(unlist(filSegs@locLikelihoods[discards[3],,drop = TRUE]) >= log(lociCutoff)))
                                        #              if(!leftExt & !rightExt) return(c(0, 0))
                                        #              if(leftExt & !rightExt) return(c(width(emptyPD@coordinates)[discards[1]], 0))
                                        #              if(!leftExt & rightExt) return(c(0, width(emptyPD@coordinates)[discards[1]]))
                                        #              if(leftExt & rightExt) {
                                        #                lenweight <- width(filSegs@coordinates)[discards[2]] + width(filSegs@coordinates)[discards[3]]
                                        #                extLeft <- round(width(emptyPD@coordinates)[discards[1]] * width(filSegs@coordinates)[discards[2]] / lenweight + runif(1, -1e-3, 1e-3) * as.numeric(width(filSegs@coordinates)[discards[2]] == width(filSegs@coordinates)[discards[3]]))
                                        #                extRight <- width(emptyPD@coordinates)[discards[1]] - extLeft
                                        #                return(c(extLeft, extRight))                
                                        #              }
                                        #            }))
            
            if(verbose) message("done!")
            
            
            extSegs <- filSegs
            
            if(length(extensions) > 0)
              {
                end(extSegs@coordinates)[leftExtendSeg[!is.na(leftExtendSeg)]] <- end(extSegs@coordinates)[leftExtendSeg[!is.na(leftExtendSeg)]] + extensions[!is.na(leftExtendSeg),1]
                
                start(extSegs@coordinates)[rightExtendSeg[!is.na(rightExtendSeg)]] <- start(extSegs@coordinates)[rightExtendSeg[!is.na(rightExtendSeg)]] - extensions[!is.na(rightExtendSeg),2]
              }
          } else extSegs <- filSegs
      } else extSegs <- filSegs

    #extSegs <- .convertSegToLoci(extSegs)

    if (getLikes & nrow(extSegs) > 1) likeSegs <- lociLikelihoods(aD = aD, cD = extSegs, cl = cl) else likeSegs <- extSegs    
    if(ncol(likeSegs@locLikelihoods) > 0) colnames(likeSegs@locLikelihoods) <- levels(likeSegs@replicates)
    likeSegs
  }
