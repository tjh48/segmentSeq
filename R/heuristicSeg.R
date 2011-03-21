heuristicSeg <- function(sD, aD, bimodality = TRUE, RKPM = 30, gap = 100, subRegion = NULL, getLikes = TRUE, verbose = TRUE, cl)
  {
    if((missing(aD) || class(aD) != "alignmentData") & getLikes)
      stop("I can't assess the likelihoods of clustered data without an alignmentData object 'aD'")
     
    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)
    
    sD <- sD[order(as.factor(sD@segInfo$chr), sD@segInfo$start, -sD@segInfo$end),]

    dupStarts <- which(fastUniques(cbind(sD@segInfo$chr, sD@segInfo$start)))                                         

    if(verbose) message("Evaluating potential loci...", appendLF = TRUE)

    replicates <- sD@replicates
    
    locDens <- do.call("cbind", lapply(unique(replicates), function(rep) {
      locDens <- (colSums(t(sD@data[,replicates == rep,drop = FALSE]) / sD@libsizes[replicates == rep]) / sum(replicates == rep) / (sD@segInfo$end - sD@segInfo$start + 1))
      locDens
    }))
                                                                              
    if(!bimodality) locM <- locDens * 1e9 > RKPM else {
      message("Finding cutoff values...", appendLF = FALSE)
      locDens <- log10(locDens)
      locCutoff <- apply(locDens, 2, function(x) {
        message(".", appendLF = FALSE)
        nzlocs <- which(x != -Inf)
        bimodalSep(x[nzlocs], bQ = c(0.5, 1))
      })
      message("done!")
      locM <- t(t(locDens) > locCutoff)
    }

    rm(locDens)
    gc()

    selLoc <- rowSums(locM) > 0
    potlociD <- new("postSeg")
    potlociD@libsizes <- sD@libsizes
    potlociD@replicates <- sD@replicates
    potlociD@seglens <- matrix(with(sD@segInfo[selLoc,], end - start + 1), ncol = 1)
    potlociD@data <- sD@data[selLoc,]
    potlociD@annotation <- subset(sD@segInfo, subset = selLoc, select = c(chr, start, end))
    potlociD@posteriors <- log(locM[selLoc,,drop = FALSE])

    emptyNulls <- with(sD@segInfo, data.frame(chr = chr[dupStarts], start = start[dupStarts] - leftSpace[dupStarts], end = start[dupStarts] - 1))
    emptyNulls <- emptyNulls[emptyNulls$start > 1, ]
    emptyNulls <- emptyNulls[emptyNulls$end - emptyNulls$start + 1 > 0,]
    
    if(bimodality) {
      emptyNulls <- emptyNulls[getOverlaps(emptyNulls, potlociD@annotation, overlapType = "within", whichOverlaps = FALSE, cl = NULL),]
      gap <- 10^(bimodalSep(log10(emptyNulls$end - emptyNulls$start + 1)))
    }

    emptyNulls <- emptyNulls[emptyNulls$end - emptyNulls$start + 1 > as.integer(floor(gap)),]

    whSDboth <- with(sD@segInfo,
                     which(leftSpace > 0 & rightSpace > 0))
    seglenBoth <- with(sD@segInfo,
                       (end - start + 1 + leftSpace + rightSpace)[whSDboth])
    whSDboth <- whSDboth[seglenBoth > as.integer(floor(gap))]
    seglenBoth <- seglenBoth[seglenBoth > as.integer(floor(gap))]

    whSDleft <- with(sD@segInfo,
                     which(leftSpace > 0))
    seglenLeft <- with(sD@segInfo,
                       (end - start + 1 + leftSpace)[whSDleft])
    whSDleft <- whSDleft[seglenLeft > as.integer(floor(gap))]
    seglenLeft <- seglenLeft[seglenLeft > as.integer(floor(gap))]

    whSDright <- with(sD@segInfo,
                     which(rightSpace > 0))
    seglenRight <- with(sD@segInfo,
                       (end - start + 1 + rightSpace)[whSDright])
    whSDright <- whSDright[seglenRight > as.integer(floor(gap))]
    seglenRight <- seglenRight[seglenRight > as.integer(floor(gap))]

    gc()

    potnullD <- with(sD@segInfo, new("postSeg",
                                     libsizes = sD@libsizes,
                                     replicates = sD@replicates,
                                     data = matrix(nrow = 0, ncol = ncol(sD))))
    nullData <- rbind(matrix(0L, ncol = ncol(sD), nrow = nrow(emptyNulls)),
                  sD@data[c(whSDboth, whSDleft, whSDright),])    
    nullSeglens = as.integer(c(emptyNulls$end - emptyNulls$start + 1, seglenBoth, seglenLeft, seglenRight))
    
    nullAnnotation = with(sD@segInfo, data.frame(
      chr = levels(chr)[c(emptyNulls$chr, chr[c(whSDboth, whSDleft, whSDright)])],
      start = c(emptyNulls$start, (start[whSDboth] - leftSpace[whSDboth]), (start[whSDleft] - leftSpace[whSDleft]), start[whSDright]),
      end = c(emptyNulls$end, (end[whSDboth] + rightSpace[whSDboth]), end[whSDleft], (end + rightSpace)[whSDright]),
      nullClass = as.factor(rep(c("empty", "segSpace"), c(nrow(emptyNulls), length(whSDboth) + length(whSDleft) + length(whSDright))))
      ))    

    rm(whSDboth, whSDleft, whSDright, seglenBoth, seglenLeft, seglenRight)
    gc()

    if(nrow(nullData) > 0)
      {
        nulDens <- do.call("cbind", lapply(unique(replicates), function(rep) {
          nulDens <- (colSums(t(nullData[,replicates == rep,drop = FALSE]) / potnullD@libsizes[replicates == rep]) / sum(replicates == rep) / (nullSeglens + 1))
          nulDens
        }))
        
        if(!bimodality) nulM <- t(t(nulDens * 1e9 < RKPM) & nullSeglens > gap) else {
          nulM <- t(t(log10(nulDens)) < locCutoff)
        }
      }

    whNull <- which(rowSums(nulM) > 0)
    nulM <- nulM[whNull,,drop = FALSE]
    nullData <- nullData[whNull,,drop = FALSE]
    nullSeglens <- nullSeglens[whNull]
    nullAnnotation <- nullAnnotation[whNull,]

    gc()
        
    nulWithin <- getOverlaps(nullAnnotation, potlociD@annotation, overlapType = "within", whichOverlaps = FALSE, cl = NULL)

    nulM <- nulM[nulWithin,,drop = FALSE]
    nullData <- nullData[nulWithin,,drop = FALSE]
    nullSeglens <- nullSeglens[nulWithin]
    nullAnnotation <- nullAnnotation[nulWithin,,drop = FALSE]

    gc()

    if (is.null(subRegion)) {
      locSub <- 1:nrow(potlociD)
      nulSub <- 1:nrow(nullAnnotation)      
    } else {
      locSub <- sort(unique(c(unlist(apply(subRegion, 1, 
                                           function(sR) which(as.character(potlociD@annotation$chr) == as.character(sR[1]) &
                                                              potlociD@annotation$start >= as.numeric(sR[2]) &
                                                              potlociD@annotation$end <= as.numeric(sR[3])))))))
      nulSub <- sort(unique(c(unlist(apply(subRegion, 1, 
                                           function(sR) which(as.character(nullAnnotation$chr) == as.character(sR[1]) &
                                                              nullAnnotation$start >= as.numeric(sR[2]) &
                                                              nullAnnotation$end <= as.numeric(sR[3])))))))
    }    

    if(nrow(nullAnnotation) == 0) nulSub <- NULL
    if(nrow(potlociD) == 0) locSub <- NULL

    potlociD <- potlociD[locSub,]
    
    potnullD@data <- nullData[nulSub,, drop = FALSE]
    potnullD@annotation <- nullAnnotation[nulSub,]
    potnullD@seglens <- matrix(nullSeglens[nulSub], ncol = 1)
    potnullD@posteriors <- nulM[nulSub,, drop = FALSE]

    rm(nullData, nullAnnotation, nullSeglens, nulM)
    
    gc()
    
    seg <- processPosteriors(potlociD, potnullD, chrs = sD@chrs, aD = aD, lociCutoff = 1, nullCutoff = 1, getLikes = getLikes, cl = cl)
    seg
  }
