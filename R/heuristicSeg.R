.zeroInMeth <- function(aD)
  {
    zeroCs <- rowSums(sapply(1:ncol(aD), function(ii) as.integer(aD@Cs[,ii]))) == 0
    
    chrBreaks <- which(seqnames(aD@alignments)[-nrow(aD)] != seqnames(aD@alignments)[-1])
    chrBreaks <- cbind(c(1, chrBreaks + 1), c(chrBreaks, nrow(aD)))
    
    whichZero <- which(zeroCs)
    
    zeroBlocks <- do.call("rbind", lapply(1:nrow(chrBreaks), function(ii) {
      chrZeros <- whichZero[whichZero >= chrBreaks[ii,1] & whichZero <= chrBreaks[ii,2]]      
      adjZeros <- cbind(chrZeros[c(1, which(diff(chrZeros) > 1) + 1)], chrZeros[c(which(diff(chrZeros) > 1), length(chrZeros))])
    }))
    
    zeroCoords <- GRanges(seqnames(aD@alignments)[zeroBlocks[,1]], IRanges(start = start(aD@alignments)[zeroBlocks[,1]], end = end(aD@alignments)[zeroBlocks[,2]]))
    zeroAD <- getMethCounts(zeroCoords, aD = aD, preFiltered = FALSE, as.matrix = FALSE, cl = cl)
    
    zeroAD    
  }      


heuristicSeg <- function(sD, aD, RKPM = 1000, gap = 100, prop = 0.2, subRegion = NULL, largeness = 1e8, getLikes = TRUE, verbose = TRUE, cl = NULL)
  {
    if(!is.null(subRegion))
      {
        sD <- sD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(sD@coordinates)) == subRegion$chr[ii] & end(sD@coordinates) >= subRegion$start[ii] & start(sD@coordinates) <= subRegion$end[ii]))),]
        aD <- aD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(aD@alignments)) == subRegion$chr[ii] & end(aD@alignments) >= subRegion$start[ii] & start(aD@alignments) <= subRegion$end[ii]))),]
      }
            
    if(prod(dim(sD)) > largeness)
      {
        sDsplit <- .splitSD(sD, largeness)
        
        message("Segmentation split into ", length(sDsplit), " parts.")
        
        splitSeg <- lapply(1:length(sDsplit), function(ii) {
          message("Segmenting; Part ", ii, " of ", length(sDsplit))
          .partheuristicSeg(sDsplit[[ii]], aD = aD, bimodality = FALSE, verbose = TRUE, cl = NULL, RKPM = RKPM, gap = gap, prop = prop)      
        })
        
        lD <- .mergeListLoci(splitSeg)
      } else lD <- .partheuristicSeg(sD, aD = aD, bimodality = FALSE, verbose = TRUE, cl = NULL, RKPM = RKPM, gap = gap, prop = prop)
        
    if(getLikes & !missing(aD)) lD <- lociLikelihoods(lD, aD, cl = cl) else if(getLikes & missing(aD)) warning("I can't calculate locus likelihoods without an aD object. You can run lociLikelihoods on the output of this function for the same result.")

    lD
  }


.partheuristicSeg <- function(sD, aD, bimodality = FALSE, RKPM = 1000, gap = 100, prop = 0.2, subRegion = NULL, verbose = TRUE, locCutoff = 0.9, cl)
  {     
    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)

    sD <- sD[order(as.factor(seqnames(sD@coordinates)), as.integer(start(sD@coordinates)), as.integer(end(sD@coordinates))),]

    dupStarts <- which(fastUniques(cbind(as.character(seqnames(sD@coordinates)), as.integer(start(sD@coordinates)))))

    if(verbose) message("Evaluating potential loci...", appendLF = TRUE)

    replicates <- sD@replicates

    seglens <- matrix(width(ranges(sD@coordinates)), ncol = 1)

    densityFunction <- function(rep, sD, seglens, RKPM, gap, bimodality)
      {
        locDens <- rowSums(sapply(which(sD@replicates == rep), function(jj)
                       as.integer(sD@data[,jj]) / seglens[,1] / sD@libsizes[jj])) / sum(sD@replicates == rep)
        message(".", appendLF = FALSE)
        
        if(!bimodality) return(locDens * 1e9 > RKPM) else {
          logDens <- log10(locDens)          
          nzlocs <- which(logDens != -Inf)
          bimodalSep(logDens[nzlocs], bQ = c(0, 1))
          message("done!")
          return(logDens > locCutoff)
        }
      }

    methFunction <- function(methD, prop, locCutoff, nullP = FALSE)
      {
        Cs <- sapply(1:ncol(methD@Cs), function(ii) as.integer(methD@Cs[,ii]))
        Ts <- sapply(1:ncol(methD@Cs), function(ii) as.integer(methD@Ts[,ii]))

        combCs <- sapply(levels(methD@replicates), function(rep) rowSums(Cs[,methD@replicates == rep,drop = FALSE]))
        combTs <- sapply(levels(methD@replicates), function(rep) rowSums(Ts[,methD@replicates == rep,drop = FALSE]))

        # Beta-distribution here is derived from being the improper (B(0,0)) conjugate prior of a binomial
        
        if(!nullP) {
          pbetas <- suppressWarnings(log(pbeta(prop, combCs, combTs, log = FALSE, lower.tail = FALSE)))        
          pbetas[combCs > 0 & combTs == 0] <- 1
          pbetas[combCs == 0 & combTs > 0] <- -Inf
        } else {
          pbetas <- suppressWarnings(log(pbeta(prop, combCs, combTs, log = FALSE, lower.tail = TRUE)))
          pbetas[combCs > 0 & combTs == 0] <- -Inf
          pbetas[combCs == 0 & combTs > 0] <- 0
        }
        
        locM <- pbetas > log(locCutoff)
      }

    if(class(sD) == "methSegs") {
      locM <- methFunction(sD, prop = prop, locCutoff = locCutoff)
      sD@locLikelihoods <- DataFrame(log(locM))

      empties <- .zeroInMeth(aD = aD)

      
      potnullD <- .constructNulls(empties@alignments, sD, sD@coordinates[rowSums(locM == 1, na.rm = TRUE) > 0,], withinOnly = FALSE)

      potnullD@Cs <- getCounts(potnullD@coordinates, new("alignmentData", data = aD@Cs, alignments = aD@alignments, libnames = aD@libnames, libsizes = rep(1, ncol(aD)), replicates = aD@replicates), cl = cl)
      
      potnullD@Ts <- getCounts(potnullD@coordinates, new("alignmentData", data = aD@Ts, alignments = aD@alignments, libnames = aD@libnames, libsizes = rep(1, ncol(aD)), replicates = aD@replicates), cl = cl)
      potnullD <- potnullD[width(potnullD@coordinates) > gap,]

      nullM <- methFunction(potnullD, prop = prop, locCutoff = locCutoff, nullP = TRUE)
      potnullD@locLikelihoods <- DataFrame(log(nullM))
      
      message("...done!")
      seg <- .processPosteriors(lociPD = sD, nullPD = potnullD, lociCutoff = 1, nullCutoff = 1, getLikes = FALSE, cl = cl)
      return(seg)

      } else {
        locM <- do.call("cbind", lapply(levels(sD@replicates), densityFunction, sD = sD, seglens = seglens, RKPM = RKPM, gap = gap, bimodality = bimodality))                
    
        selLoc <- rowSums(locM, na.rm = TRUE) > 0
        potlociD <- sD[selLoc,]
    
        potlociD@replicates <- as.factor(sD@replicates)    
        potlociD@locLikelihoods <- DataFrame(log(locM[selLoc,,drop = FALSE]))
        
        if(class(sD) == "segData") {
          potnullD <- new("segData",
                          libsizes = sD@libsizes,
                          replicates = sD@replicates,
                          data = DataFrame(matrix(nrow = 0, ncol = ncol(sD))))
          emptyPD <- new("segData", libsizes = sD@libsizes, replicates = sD@replicates, data = DataFrame(matrix(nrow = 0, ncol = ncol(sD))))
        } else if(class(sD) == "methSegs") {
          potnullD <- new("methSegs",
                          replicates = sD@replicates,
                          Cs = DataFrame(matrix(nrow = 0, ncol = ncol(sD))),
                          Ts = DataFrame(matrix(nrow = 0, ncol = ncol(sD))))
          emptyPD <- new("methSegs", replicates = sD@replicates,
                         Cs = DataFrame(matrix(nrow = 0, ncol = ncol(sD))),
                         Ts = DataFrame(matrix(nrow = 0, ncol = ncol(sD))))
        }
        
        
        
                                        #    if(examineNulls) {
        emptyNulls <- gaps(sD@coordinates[dupStarts,])
        emptyNulls <- emptyNulls[strand(emptyNulls) == "*",]    
        
        if(bimodality) {
          emptyNulls <- emptyNulls[getOverlaps(emptyNulls, potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL),]
          gap <- 10^(bimodalSep(log10(end(emptyNulls) - start(emptyNulls) + 1)))
        }
        
        emptyPD@coordinates <- emptyNulls[width(emptyNulls) <= as.integer(floor(gap)),]
        emptyPD@locLikelihoods <- DataFrame(matrix(-Inf, ncol = length(levels(sD@replicates)), nrow = length(emptyPD@coordinates)))
        
        emptyPD@data <- DataFrame(matrix(0, ncol = length(sD@replicates), nrow = length(emptyPD@coordinates)))
        
        emptyNulls <- emptyNulls[width(emptyNulls) > as.integer(floor(gap)),]
        whover <- getOverlaps(emptyNulls, potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)
        emptyNulls <- emptyNulls[whover,]
        
        expandRanges <- suppressWarnings(do.call("c", lapply(seqlevels(sD@coordinates), function(ii) {
          subRange <- ranges(sD@coordinates[seqnames(sD@coordinates) == ii,])
          right <- c(start(subRange), 1 + seqlengths(sD@coordinates)[ii])[findInterval(end(subRange), c(start(subRange), 1 + seqlengths(sD@coordinates)[ii])) + 1] - 1
          left <- c(0, unique(end(subRange)))[findInterval(start(subRange), unique(end(subRange)) + 0.5) + 1] + 1      
          
          
          if(length(subRange) > 0) {
            chrRanges <- GRanges(seqnames = ii, IRanges(start = left, end = right))
            chlevels <- rep(NA, length(seqlevels(sD@coordinates)))
            chlevels[seqlevels(sD@coordinates) == ii] <- 1L
            seqinfo(chrRanges, new2old = chlevels) <- seqinfo(sD@coordinates)
          } else {
            chrRanges <- GRanges()
            seqinfo(chrRanges) <- seqinfo(sD@coordinates)
          }
          chrRanges
        })))
        
        whSDboth <- which(start(expandRanges) != start(sD@coordinates) & end(expandRanges) != end(sD@coordinates))
        whSDboth <- whSDboth[getOverlaps(expandRanges[whSDboth,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
        bothAnn <- expandRanges[whSDboth,]
        whSDboth <- whSDboth[width(bothAnn) > gap]
        bothAnn <- bothAnn[width(bothAnn) > gap,]
        
        whSDleft <- which(start(expandRanges) != start(sD@coordinates))
        whSDleft <- whSDleft[getOverlaps(sD@coordinates[whSDleft,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
        leftAnn <- expandRanges[whSDleft,]
        end(leftAnn) <- end(sD@coordinates[whSDleft])
        whSDleft <- whSDleft[width(leftAnn) > gap]
        leftAnn <- leftAnn[width(leftAnn) > gap,]
        
        whSDright <- which(end(expandRanges) != end(sD@coordinates))
        whSDright <- whSDright[getOverlaps(sD@coordinates[whSDright,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
        rightAnn <- expandRanges[whSDright,]
        start(rightAnn) <- start(sD@coordinates[whSDright])
        whSDright <- whSDright[width(rightAnn) > gap]
        rightAnn <- rightAnn[width(rightAnn) > gap,]
        
        gc()
        
        zeroData <- do.call("DataFrame", sapply(1:ncol(sD), function(x) Rle(0, length(emptyNulls))))
        colnames(zeroData) <- colnames(sD@data)
        
        nullData <- rbind(zeroData,
                          sD@data[whSDboth,,drop = FALSE],
                          sD@data[whSDleft,,drop = FALSE],
                          sD@data[whSDright,,drop = FALSE])
        nullAnnotation <- c(emptyNulls, bothAnn, leftAnn, rightAnn)
        
        rm(whSDboth, whSDleft, whSDright, bothAnn, leftAnn, rightAnn, emptyNulls, expandRanges, zeroData)
        gc()      
        
        if(nrow(nullData) > 0)
          {
            nulDens <- do.call("cbind", lapply(levels(replicates), function(rep)       
                                               rowSums(sapply(which(replicates == rep), function(jj)
                                                              as.integer(nullData[,jj]) / width(nullAnnotation) / sD@libsizes[jj])) / sum(replicates == rep)
                                               ))
            
            if(!bimodality) nulM <- (nulDens * 1e9 < RKPM) & width(nullAnnotation) > gap else {
              nulM <- (log10(nulDens) < locCutoff) & width(nullAnnotation) > gap
            }
          }
        
        whNull <- which(rowSums(nulM) > 0)
        nulM <- nulM[whNull,,drop = FALSE]
        nullData <- nullData[whNull,,drop = FALSE]
        nullAnnotation <- nullAnnotation[whNull,]
        gc()
        
        if (is.null(subRegion)) {
          nulSub <- 1:length(nullAnnotation)
        } else nulSub <- sort(unique(c(unlist(apply(subRegion, 1, 
                                                    function(sR) which(seqnames(nullAnnotation) == as.character(sR[1]) &
                                                                       start(nullAnnotation) >= as.numeric(sR[2]) &
                                                                       end(nullAnnotation) <= as.numeric(sR[3])))))))
        
        if(length(nullAnnotation) == 0) nulSub <- NULL
        potnullD@data <- nullData[nulSub,,drop = FALSE]
        potnullD@coordinates <- nullAnnotation[nulSub,]
        potnullD@locLikelihoods <- DataFrame(log(nulM[nulSub,, drop = FALSE]))
        
        rm(nullData, nullAnnotation, nulM)
        
        if (is.null(subRegion)) {
          locSub <- 1:nrow(potlociD)      
        } else locSub <- sort(unique(c(unlist(apply(subRegion, 1, 
                                                    function(sR) which(as.character(seqnames(potlociD@coordinates)) == as.character(sR[1]) &
                                                                       start(potlociD@coordinates) >= as.numeric(sR[2]) &
                                                                       end(potlociD@coordinates) <= as.numeric(sR[3])))))))
        
        if(nrow(potlociD) == 0) locSub <- NULL
        
        potlociD <- potlociD[locSub,]
        
        gc()
        
        seg <- .processPosteriors(lociPD = potlociD, nullPD = potnullD, emptyPD = emptyPD, aD = aD, lociCutoff = 1, nullCutoff = 1, getLikes = FALSE, cl = cl)
      }
    seg
  }
