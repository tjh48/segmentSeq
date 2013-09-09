heuristicSeg <- function(sD, aD, gap = 100, RKPM = 1000, prop = 0.2, locCutoff = 0.99, subRegion = NULL, largeness = 1e8, getLikes = TRUE, verbose = TRUE, cl = NULL)
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
          .partheuristicSeg(sD[sDsplit[[ii]],], aD = aD, bimodality = FALSE, verbose = TRUE, cl = NULL, RKPM = RKPM, gap = gap, prop = prop, locCutoff = locCutoff)      
        })
        
        lD <- .mergeListLoci(splitSeg)
      } else lD <- .partheuristicSeg(sD, aD = aD, bimodality = FALSE, verbose = TRUE, cl = NULL, RKPM = RKPM, gap = gap, prop = prop, locCutoff = locCutoff)
        
    if(getLikes & !missing(aD)) lD <- lociLikelihoods(lD, aD, cl = cl) else if(getLikes & missing(aD)) warning("I can't calculate locus likelihoods without an aD object. You can run lociLikelihoods on the output of this function for the same result.")

    lD
  }


.partheuristicSeg <- function(sDP, aD, bimodality = FALSE, RKPM = 1000, gap = 100, prop = 0.2, subRegion = NULL, verbose = TRUE, locCutoff, cl)
  {
    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)

    sDP <- sDP[order(as.factor(seqnames(sDP@coordinates)), as.integer(start(sDP@coordinates)), as.integer(end(sDP@coordinates))),]
    sDPSmall <- fastUniques(data.frame(chr = as.character(seqnames(sDP@coordinates)), start = as.numeric(start(sDP@coordinates))))
    
    dupStarts <- which(fastUniques(cbind(as.character(seqnames(sDP@coordinates)), as.integer(start(sDP@coordinates)))))

    if(verbose) message("Evaluating potential loci...", appendLF = FALSE)

    replicates <- sDP@replicates

    seglens <- matrix(width(ranges(sDP@coordinates)), ncol = 1)

    densityFunction <- function(rep, sDP, seglens, RKPM, gap, bimodality)
      {
        locDens <- rowSums(sapply(which(sDP@replicates == rep), function(jj)
                       as.integer(sDP@data[,jj]) / seglens[,1] / sDP@libsizes[jj])) / sum(sDP@replicates == rep)
        message(".", appendLF = FALSE)
        
        if(!bimodality) return(locDens * 1e9 > RKPM) else {
          logDens <- log10(locDens)          
          nzlocs <- which(logDens != -Inf)
          bimodalSep(logDens[nzlocs], bQ = c(0, 1))
          message("done!")
          return(logDens > locCutoff)
        }
      }

    if(class(sDP) == "segMeth") {

      sDP@locLikelihoods <- (.methFunction(sDP, prop = prop, locCutoff = locCutoff))
      
      internalNulls <- sDP[rowSums(!sDP@locLikelihoods, na.rm = TRUE) > 0 & width(sDP@coordinates) > gap,]

      message(".", appendLF = FALSE)
      empties <- .zeroInMeth(aD = aD, smallSegs = sDP@coordinates[sDPSmall])
      message(".", appendLF = FALSE)
      if(length(empties) > 0) {
        potnullD <- .constructMethNulls(empties, sDP, sDP@coordinates[rowSums(sDP@locLikelihoods, na.rm = TRUE) > 0,])
        message(".", appendLF = FALSE)
        potnullD <- potnullD[width(potnullD@coordinates) > gap,]

        counts <- getCounts(potnullD@coordinates, aD, cl = cl)
        
        potnullD@Cs <- counts$Cs
        potnullD@Ts <- counts$Ts

        colnames(potnullD@Cs) <- colnames(internalNulls@Cs)
        colnames(potnullD@Ts) <- colnames(internalNulls@Cs)
        
      } else potnullD <- new("segMeth")      

      message(".", appendLF = FALSE)
      potnulls <- new("segMeth",
                      Cs = rbind(potnullD@Cs, internalNulls@Cs),
                      Ts = rbind(potnullD@Ts, internalNulls@Ts),
                      replicates = internalNulls@replicates,
                      coordinates = c(potnullD@coordinates, internalNulls@coordinates), nonconversion = sDP@nonconversion)
      
      potnulls@locLikelihoods <- log(.methFunction(potnulls, prop = prop, locCutoff = locCutoff, nullP = TRUE))
      

#      save(sDP, potnulls, file = "temp/processPosteriorVariables.RData")

      sDP@locLikelihoods <- log(sDP@locLikelihoods)
      
      message("...done!")
      seg <- .processPosteriors(lociPD = sDP, nullPD = potnulls, lociCutoff = 1, nullCutoff = 1, getLikes = FALSE, cl = cl)
      gc()
      return(seg)

    } else {
      locM <- do.call("cbind", lapply(levels(sDP@replicates), densityFunction, sDP = sDP, seglens = seglens, RKPM = RKPM, gap = gap, bimodality = bimodality))                
      
      selLoc <- rowSums(locM, na.rm = TRUE) > 0
      potlociD <- sDP[selLoc,]
      
      potlociD@replicates <- as.factor(sDP@replicates)    
      potlociD@locLikelihoods <- (log(locM[selLoc,,drop = FALSE]))
      
      potnullD <- new("segData",
                      libsizes = sDP@libsizes,
                      replicates = sDP@replicates,
                      data = (matrix(nrow = 0, ncol = ncol(sDP))))
      emptyPD <- new("segData", libsizes = sDP@libsizes, replicates = sDP@replicates, data = (matrix(nrow = 0, ncol = ncol(sDP))))
        
                                        #    if(examineNulls) {
      emptyNulls <- gaps(sDP@coordinates[dupStarts,])
      emptyNulls <- emptyNulls[strand(emptyNulls) == "*",]    
      
      if(bimodality) {
        emptyNulls <- emptyNulls[getOverlaps(emptyNulls, potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL),]
        gap <- 10^(bimodalSep(log10(end(emptyNulls) - start(emptyNulls) + 1)))
      }
        
        emptyPD@coordinates <- emptyNulls[width(emptyNulls) <= as.integer(floor(gap)),]
        emptyPD@locLikelihoods <- (matrix(-Inf, ncol = length(levels(sDP@replicates)), nrow = length(emptyPD@coordinates)))
        
        emptyPD@data <- (matrix(0, ncol = length(sDP@replicates), nrow = length(emptyPD@coordinates)))
        
        emptyNulls <- emptyNulls[width(emptyNulls) > as.integer(floor(gap)),]
        whover <- getOverlaps(emptyNulls, potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)
        emptyNulls <- emptyNulls[whover,]
        
        expandRanges <- suppressWarnings(do.call("c", lapply(seqlevels(sDP@coordinates), function(ii) {
          subRange <- ranges(sDP@coordinates[seqnames(sDP@coordinates) == ii,])
          right <- c(start(subRange), 1 + seqlengths(sDP@coordinates)[ii])[findInterval(end(subRange), c(start(subRange), 1 + seqlengths(sDP@coordinates)[ii])) + 1] - 1
          left <- c(0, unique(end(subRange)))[findInterval(start(subRange), unique(end(subRange)) + 0.5) + 1] + 1      
          
          
          if(length(subRange) > 0) {
            chrRanges <- GRanges(seqnames = ii, IRanges(start = left, end = right))
            chlevels <- rep(NA, length(seqlevels(sDP@coordinates)))
            chlevels[seqlevels(sDP@coordinates) == ii] <- 1L
            seqinfo(chrRanges, new2old = chlevels) <- seqinfo(sDP@coordinates)
          } else {
            chrRanges <- GRanges()
            seqinfo(chrRanges) <- seqinfo(sDP@coordinates)
          }
          chrRanges
        })))
        
        whSDboth <- which(start(expandRanges) != start(sDP@coordinates) & end(expandRanges) != end(sDP@coordinates))
        whSDboth <- whSDboth[getOverlaps(expandRanges[whSDboth,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
        bothAnn <- expandRanges[whSDboth,]
        whSDboth <- whSDboth[width(bothAnn) > gap]
        bothAnn <- bothAnn[width(bothAnn) > gap,]
        
        whSDleft <- which(start(expandRanges) != start(sDP@coordinates))
        whSDleft <- whSDleft[getOverlaps(sDP@coordinates[whSDleft,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
        leftAnn <- expandRanges[whSDleft,]
        end(leftAnn) <- end(sDP@coordinates[whSDleft])
        whSDleft <- whSDleft[width(leftAnn) > gap]
        leftAnn <- leftAnn[width(leftAnn) > gap,]
        
        whSDright <- which(end(expandRanges) != end(sDP@coordinates))
        whSDright <- whSDright[getOverlaps(sDP@coordinates[whSDright,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
        rightAnn <- expandRanges[whSDright,]
        start(rightAnn) <- start(sDP@coordinates[whSDright])
        whSDright <- whSDright[width(rightAnn) > gap]
        rightAnn <- rightAnn[width(rightAnn) > gap,]
        
        gc()
        
        zeroData <- do.call("cbind", lapply(1:ncol(sDP), function(x) rep(0L, length(emptyNulls))))
        colnames(zeroData) <- colnames(sDP@data)
        
        nullData <- rbind(zeroData,
                          sDP@data[whSDboth,,drop = FALSE],
                          sDP@data[whSDleft,,drop = FALSE],
                          sDP@data[whSDright,,drop = FALSE])
        nullAnnotation <- c(emptyNulls, bothAnn, leftAnn, rightAnn)
        
        rm(whSDboth, whSDleft, whSDright, bothAnn, leftAnn, rightAnn, emptyNulls, expandRanges, zeroData)
        gc()      
        
        if(nrow(nullData) > 0)
          {
            nulDens <- do.call("cbind", lapply(levels(replicates), function(rep)       
                                               rowSums(sapply(which(replicates == rep), function(jj)
                                                              as.integer(nullData[,jj]) / width(nullAnnotation) / sDP@libsizes[jj])) / sum(replicates == rep)
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
        potnullD@locLikelihoods <- (log(nulM[nulSub,, drop = FALSE]))
        
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
