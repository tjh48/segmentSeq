heuristicSeg <- function(sD, aD, gap = 50, RKPM = 1000, prop, locCutoff = 0.9, nullCutoff = 0.9, subRegion = NULL, largeness = 1e8, getLikes = TRUE, verbose = TRUE, tempDir = NULL, cl = NULL, recoverFromTemp = FALSE, trimMeth = FALSE)
  {
    if(!is.null(tempDir)) dir.create(tempDir, showWarnings = FALSE)

    if(inherits(aD, "alignmentMeth") && (missing(prop) || !is.numeric(prop)))
        stop("'prop variable must be defined as a numeric variable; see 'thresholdFinder' function.")
        
    
    if(!is.null(subRegion))
      {
        sD <- sD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(sD@coordinates)) == subRegion$chr[ii] & end(sD@coordinates) >= subRegion$start[ii] & start(sD@coordinates) <= subRegion$end[ii]))),]
        aD <- aD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(aD@alignments)) == subRegion$chr[ii] & end(aD@alignments) >= subRegion$start[ii] & start(aD@alignments) <= subRegion$end[ii]))),]
      }
                              
    if(prod(dim(sD)) > largeness)
      {
        if(recoverFromTemp & !is.null(tempDir) & file.exists(paste(tempDir, "/sDsplit.RData", sep = ""))) {
          if(verbose) message("Recovering splitting from temporary file...")
          load(paste(tempDir, "/sDsplit.RData", sep = ""))
        } else {
            sDsplit <- .splitSD(sD, largeness)
            if(!is.null(tempDir))
                save(sDsplit, file = paste(tempDir, "/sDsplit.RData", sep = ""))
        }
        
        if(verbose) message("Segmentation split into ", length(sDsplit), " parts.")

        splitSeg <- lapply(1:length(sDsplit), function(ii) {
            if(verbose) message("Segmenting; Part ", ii, " of ", length(sDsplit))
            
            strandSeg <- lapply(unique(strand(sD@coordinates)), function(strand) {
                if(verbose) message("Strand: ", strand)
                existingSDP <- FALSE
                                        # better name recognition for temporary files from sD[sDsplit[[ii]]] ranges needed              
              if(recoverFromTemp & !is.null(tempDir))
                  {
                      existingTemp <- dir(tempDir, pattern = paste("seg", ii, "_\\", strand, "_", sep = ""), full.names = TRUE)
                      if(length(existingTemp) == 1) {
                          if(verbose) message("Recovering from temporary file...")
                          load(existingTemp)
                          existingSDP = TRUE
                      }
                  }
# sDP = sD[intersect(which(strand(sD@coordinates) == strand), sDsplit[[ii]]),]; aDP = aD[strand(aD@alignments) == strand,]; bimodality = FALSE; verbose = verbose; cl = NULL; RKPM = RKPM; gap = gap; prop = prop; locCutoff = locCutoff; largeness = largeness; tempDir = tempDir
                if(!existingSDP) {
                    sDP <- .partheuristicSeg(
                        sDP = sD[intersect(which(strand(sD@coordinates) == strand), sDsplit[[ii]]),],
                      aDP = aD,
                      bimodality = FALSE, verbose = verbose, cl = cl, RKPM = RKPM, gap = gap, prop = prop, locCutoff = locCutoff, nullCutoff = nullCutoff, largeness = largeness, tempDir = tempDir)
                  
                  if(!is.null(tempDir)) save(sDP,
                                             file = paste(tempDir, "/seg", ii, "_", strand, "_", "largeness_", largeness,
                                        #paste(apply(apply(as.data.frame(range(sDP@coordinates, ignore.strand = FALSE))[,1:3], 1, as.character), 2, function(x) paste(gsub(" ", "", x), collapse = "_")), collapse = "__"),
                                                 "_locLikes.RData", sep = "")
                                             )
              }
              return(sDP)
          })
            if(length(strandSeg) == 1) sDP <- strandSeg[[1]] else sDP <- .mergeListLoci(strandSeg)
          return(sDP)
      })        
        lD <- .mergeListLoci(splitSeg)
    } else {
        lD <- .partheuristicSeg(sDP = sD, aDP = aD, bimodality = FALSE, verbose = verbose, cl = cl, RKPM = RKPM, gap = gap, prop = prop, locCutoff = locCutoff, nullCutoff = nullCutoff, largeness = largeness, tempDir = tempDir)
    }
    if(class(aD) == "alignmentMeth")
        {
            if(trimMeth) lD <- .trimMeth(lD, aD, prop, cl)
            if(nrow(sD@data) == 0) {
                gcts <- getCounts(lD@coordinates, aD, cl = cl)
                lD@data <- array(c(gcts$Cs, gcts$Ts), dim = c(dim(gcts$Cs), 2))
            }
            lD@data <- round(lD@data)
        } else {
            if(nrow(sD@data) == 0)
                lD@data <- round(getCounts(lD@coordinates, aD, cl = cl))
        }    
    
    if(getLikes & !missing(aD)) lD <- lociLikelihoods(lD, aD, cl = cl) else if(getLikes & missing(aD)) warning("I can't calculate locus likelihoods without an aD object. You can run lociLikelihoods on the output of this function for the same result.")    
    
    lD
  }


.partheuristicSeg <- function(sDP, aDP, bimodality = FALSE, RKPM = 1000, gap = 100, prop = 0.2, subRegion = NULL, verbose = TRUE, locCutoff, nullCutoff, largeness = 1e8, tempDir = tempDir, cl)
  {
      fastUniques <- function(x)
          if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)
      
      if(nrow(sDP) == 0) {
                                        #sDP <- .convertSegToLoci(sDP)
          sDP@locLikelihoods <- matrix(nrow = 0, ncol = length(levels(sDP@replicates)))
          return(sDP)
      }
      sDP <- sDP[order(as.factor(seqnames(sDP@coordinates)), as.integer(start(sDP@coordinates)), as.integer(end(sDP@coordinates))),]
      if(verbose) message("Number of candidate loci: ", nrow(sDP), appendLF = TRUE)
      sDPSmall <- which(fastUniques(data.frame(chr = as.character(seqnames(sDP@coordinates)), start = as.numeric(start(sDP@coordinates)))))
      
      dupStarts <- which(fastUniques(cbind(as.character(seqnames(sDP@coordinates)), as.integer(start(sDP@coordinates)))))

      replicates <- sDP@replicates

      seglens <- matrix(width(ranges(sDP@coordinates)), ncol = 1)

      densityFunction <- function(rep, sDP, seglens, RKPM, gap, bimodality)
          {
              locDens <- rowSums(sapply(which(sDP@replicates == rep), function(jj)
                  as.integer(sDP@data[,jj]) / seglens[,1] / libsizes(sDP)[jj])) / sum(sDP@replicates == rep)
                                        #if(verbose) message(".", appendLF = FALSE)
              
              if(!bimodality) return(locDens * 1e9 > RKPM) else {
                  logDens <- log10(locDens)          
                  nzlocs <- which(logDens != -Inf)
                  bimodalSeparator(logDens[nzlocs], minperc = 10)
                  if(verbose) message("...done!")
                  return(Rle(logDens > locCutoff))
              }
          }
      
      if(class(aDP) == "alignmentMeth") {
                                        #aDP <- normaliseNC(aDP)
          breaks <- ceiling(prod(dim(sDP)) / largeness)
          if(breaks > 1) splitCalc <- split(1:nrow(sDP), cut(1:nrow(sDP), breaks = breaks, labels = FALSE)) else splitCalc <- list(1:nrow(sDP))        
          
          nullLikes <- NULL
          locLikes <- NULL
          
          if(verbose & length(splitCalc) > 1) message("")
          if(verbose & length(splitCalc) > 1) message("Splitting calculation into ", length(splitCalc), " steps...", appendLF = FALSE)
          
        for(ii in 1:length(splitCalc)) {
            if(verbose & length(splitCalc) > 1) message(ii, appendLF = FALSE)
            if(verbose & length(splitCalc) > 1) message(".", appendLF = FALSE)
            sDPx <- sDP[splitCalc[[ii]],]
            if(nrow(sDPx@data) == 0) {
                counts <- getCounts(sDPx@coordinates, aD = aDP, cl = cl)
                sDPx@data <- array(c(counts$Cs, counts$Ts), dim = c(dim(counts$Cs), 2))
            }
            lL <- .methFunction(sDPx, prop = prop, locCutoff = locCutoff)
            locLikes <- rbind(locLikes, lL)
            nL <- !lL
            nL[rowSums(nL) > 0,] <- .methFunction(sDPx[rowSums(nL) > 0,], prop = prop, locCutoff = nullCutoff, nullP = TRUE)
            nullLikes <- rbind(nullLikes, nL)
        }
                                        # sDP@locLikelihoods <- locLikes      
        
                                        #      internalNulls <- sDP[rowSums(!sDP@locLikelihoods, na.rm = TRUE) > 0 & width(sDP@coordinates) > gap,]
        
        empties <- .zeroInMeth(aD = aDP, smallSegs = sDP@coordinates[sDPSmall])
        
        if(length(empties) > 0) {
            potnullD <- .constructMethNulls(emptyNulls = empties, sDP = sDP, locDef = sDP@coordinates[which(rowSums(locLikes) > 0),], minlen = gap)
            
            if(verbose) message("Number of candidate nulls: ", nrow(potnullD), appendLF = FALSE)
            
            if(verbose) message(".", appendLF = FALSE)
                                        #        potnullD <- potnullD[width(potnullD@coordinates) > gap,]
            
            breaks <- ceiling(prod(dim(potnullD)) / largeness)
            if(breaks > 1) splitCalc <- split(1:nrow(potnullD), cut(1:nrow(potnullD), breaks = breaks, labels = FALSE)) else splitCalc <- list(1:nrow(potnullD))        

            if(nrow(potnullD) > 0) {
                potnullD@locLikelihoods <- do.call("rbind", lapply(1:length(splitCalc), function(ii) {
                    if(verbose) message(".", appendLF = FALSE)
                    potnullDx <- potnullD[splitCalc[[ii]],]
                    if(nrow(potnullDx@data) == 0) {
                        counts <- getCounts(potnullDx@coordinates, aD = aDP, cl = cl)
                        potnullDx@data <- array(c(counts$Cs, counts$Ts), dim = c(dim(counts$Cs), 2))
                    }          
                    
                    colnames(potnullDx@data) <- colnames(sDP@data)
                    potnullDx@sampleObservables$nonconversion <- aDP@nonconversion
                    ll <- .methFunction(potnullDx, prop = prop, locCutoff = nullCutoff, nullP = TRUE)
                }))
                colnames(potnullD@locLikelihoods) <- colnames(nullLikes)        
            }                                             
        } else potnullD <- new("lociData", locLikelihoods = matrix(NA, nrow = 0, ncol = nlevels(sDP@replicates)))
      
          if(verbose) message(".", appendLF = FALSE)
          potnullD <- new("lociData",
                          replicates = sDP@replicates,
                          coordinates = c(sDP@coordinates, potnullD@coordinates), sampleObservables = sDP@sampleObservables, locLikelihoods = log(rbind(nullLikes, potnullD@locLikelihoods)))

        #potnullD@locLikelihoods <- #do.call("cbind", lapply(as.list(potnullD@locLikelihoods), function(x) log(x)))
          sDP@locLikelihoods <- log(locLikes) #do.call("cbind", lapply(as.list(locLikes), function(x) log(x)))      
          
          if(!is.null(tempDir))
              save(sDP, potnullD,
                   file = paste(tempDir, "/processPosteriorVariables.RData", sep = ""))
          
          if(verbose) message("...done!")
          seg <- .processPosteriors(lociPD = sDP, nullPD = potnullD, lociCutoff = 1, nullCutoff = 1, getLikes = FALSE, extendLoci = FALSE, verbose = verbose, cl = cl)
          seg@sampleObservables <- sDP@sampleObservables
          gc()
          
          return(seg)
          
      } else {
          if(nrow(sDP@data) == 0)
              sDP@data <- getCounts(sDP@coordinates, aDP, cl = cl)

      if(verbose) message("Evaluating candidate loci...", appendLF = FALSE)
      
      locM <- do.call("DataFrame", lapply(levels(sDP@replicates), densityFunction, sDP = sDP, seglens = seglens, RKPM = RKPM, gap = gap, bimodality = bimodality))                
      if(verbose) message("done.")
      
      selLoc <- .rowSumDF(locM, na.rm = TRUE) > 0
      potlociD <- sDP[selLoc,]

      if(nrow(potlociD) == 0) return(potlociD)#return(.convertSegToLoci(potlociD))
      
      potlociD@replicates <- as.factor(sDP@replicates)    
      potlociD@locLikelihoods <- do.call("cbind", lapply(as.list(locM[selLoc,,drop = FALSE]), log))
      
      potnullD <- new("lociData",
                      libsizes = libsizes(sDP),
                      replicates = sDP@replicates,
                      data = (matrix(nrow = 0, ncol = ncol(sDP))))
      emptyPD <- new("lociData", libsizes = libsizes(sDP), replicates = sDP@replicates, data = (matrix(nrow = 0, ncol = ncol(sDP))))
        
                                        #    if(examineNulls) {
      emptyNulls <- gaps(sDP@coordinates[dupStarts,])
      emptyNulls <- emptyNulls[strand(emptyNulls) == "*",]    
      
      if(bimodality) {
                                        
          emptyNulls <- emptyNulls[!is.na(findOverlaps(emptyNulls, potlociD@coordinates, type = "within", select = "first")),]
                                        #getOverlaps(emptyNulls, potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL),]
        gap <- 10^(bimodalSeparator(log10(end(emptyNulls) - start(emptyNulls) + 1)))
      }
        
      emptyPD@coordinates <- emptyNulls[width(emptyNulls) <= as.integer(floor(gap)),]
      emptyPD@locLikelihoods <- do.call("cbind", lapply(1:length(levels(sDP@replicates)), function(jj) rep(-Inf, length(emptyPD@coordinates))))
      colnames(emptyPD@locLikelihoods) <- levels(sDP@replicates)
        
      emptyPD@data <- (matrix(0, ncol = length(sDP@replicates), nrow = length(emptyPD@coordinates)))
      
      emptyNulls <- emptyNulls[width(emptyNulls) > as.integer(floor(gap)),]
      whover <- !is.na(findOverlaps(emptyNulls, potlociD@coordinates, type = "within", select = "first"))
                                        #whover <- getOverlaps(emptyNulls, potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)
      emptyNulls <- emptyNulls[whover,]
      
      expandRanges <- suppressWarnings(do.call("c", lapply(seqlevels(sDP@coordinates), function(ii) {
          message(ii)
          subRange <- ranges(sDP@coordinates[seqnames(sDP@coordinates) == ii,])
          right <- c(start(subRange), 1 + max(end(subRange)))[findInterval(end(subRange), c(start(subRange), 1 + max(end(subRange)))) + 1] - 1
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
                                        #whSDboth <- whSDboth[getOverlaps(expandRanges[whSDboth,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
      whSDboth <- whSDboth[!is.na(findOverlaps(expandRanges[whSDboth,], potlociD@coordinates, type = "within", select = "first"))]
      bothAnn <- expandRanges[whSDboth,]
      whSDboth <- whSDboth[width(bothAnn) > gap]
      bothAnn <- bothAnn[width(bothAnn) > gap,]
        
      whSDleft <- which(start(expandRanges) != start(sDP@coordinates))
                                        #whSDleft <- whSDleft[getOverlaps(sDP@coordinates[whSDleft,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
      whSDleft <- whSDleft[!is.na(findOverlaps(sDP@coordinates[whSDleft,], potlociD@coordinates, type = "within", select = "first"))]
      leftAnn <- expandRanges[whSDleft,]
      end(leftAnn) <- end(sDP@coordinates[whSDleft])
      whSDleft <- whSDleft[width(leftAnn) > gap]
      leftAnn <- leftAnn[width(leftAnn) > gap,]
      
      whSDright <- which(end(expandRanges) != end(sDP@coordinates))
                                        #whSDright <- whSDright[getOverlaps(sDP@coordinates[whSDright,], potlociD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL)]
      whSDright <- whSDright[!is.na(findOverlaps(sDP@coordinates[whSDright,], potlociD@coordinates, type = "within", select = "first"))]
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
                  (rowSums(sapply(which(replicates == rep), function(jj)
                      as.integer(nullData[,jj]) / width(nullAnnotation) / libsizes(sDP)[jj])) / sum(replicates == rep))
                                                 ))
              
              if(!bimodality) nulM <- (nulDens * 1e9 < RKPM) & width(nullAnnotation) > gap else {
                  nulM <- (log10(nulDens) < locCutoff) & width(nullAnnotation) > gap
              }
              whNull <- which(rowSums(nulM) > 0)
          } else {
              whNull <- integer()
              nulM <- matrix(nrow = 0, ncol = length(levels(replicates)))
          }
      
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
                                        #print(colnames(nulM))
                                        #potnullD@locLikelihoods <- .matrix2Rle(log(nulM[nulSub,, drop = FALSE]))
      potnullD@locLikelihoods <- log(nulM[nulSub,, drop = FALSE])
      
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
      
      seg <- .processPosteriors(lociPD = potlociD, nullPD = potnullD, emptyPD = emptyPD, aD = aDP, lociCutoff = 1, nullCutoff = 1, extendLoci = TRUE, getLikes = FALSE, verbose = verbose, cl = cl)
  }
    seg
}


.trimMeth <- function(hSut, aM, prop, cl) {

    aMU <- aM@alignments
    values(aMU) <- NULL
    aMU <- sort(unique(aMU))

    repeat {
        startHS <- new("lociData",
                       coordinates = GRanges(seqnames(hSut@coordinates), IRanges(start(hSut@coordinates), start(hSut@coordinates)), strand = strand(hSut@coordinates)),
                       replicates = replicates(hSut))
        endHS <- new("lociData",
                     coordinates = GRanges(seqnames(hSut@coordinates), IRanges(end(hSut@coordinates), end(hSut@coordinates)), strand = strand(hSut@coordinates)),
                     replicates = replicates(hSut))
                
        startHS@sampleObservables$nonconversion = aM@nonconversion
        startCount <- getCounts(startHS@coordinates, aM, cl = cl)
        startHS@data <- array(c(startCount$Cs, startCount$Ts), dim = c(dim(startCount$Cs), 2))
                              
        
        endHS@sampleObservables$nonconversion = aM@nonconversion
        endCount <- getCounts(endHS@coordinates, aM, cl = cl)
        endHS@data <- array(c(endCount$Cs, endCount$Ts), dim = c(dim(endCount$Cs), 2))        
        
        mFS <- .methFunction(startHS, prop = prop, locCutoff = 0.5, nullP = TRUE)
        mFE <- .methFunction(endHS, prop = prop, locCutoff = 0.5, nullP = TRUE)
        trimS <- rowSums(sapply(1:nlevels(hSut@replicates), function(rr) {
            mFSrr <- as.vector(mFS[,rr])
            mFSrr[is.na(mFSrr)] <- (hSut@locLikelihoods[is.na(mFSrr),rr] == 0)
            (mFSrr & hSut@locLikelihoods[,rr] == 0)})) == rowSums(hSut@locLikelihoods == 0)

        trimE <- rowSums(sapply(1:nlevels(hSut@replicates), function(rr) {
            mFErr <- as.vector(mFE[,rr])
            mFErr[is.na(mFErr)] <- (hSut@locLikelihoods[is.na(mFErr),rr] == 0)
            (mFErr & hSut@locLikelihoods[,rr] == 0)})) == rowSums(hSut@locLikelihoods == 0)    
        
        if(sum(trimS) > 0)
            start(hSut@coordinates)[which(trimS)] <-
                start(aMU[match(startHS@coordinates, aMU)[trimS] + 1])
        
        if(sum(trimE) > 0)
            end(hSut@coordinates)[which(trimE)] <-
                end(aMU[match(endHS@coordinates, aMU)[trimE] - 1])
        if(sum(trimS) + sum(trimE) == 0) break
    }
    hSut
}
