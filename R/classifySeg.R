.convertSegToLoci <- function(sD)
  {
    lD <- new("lociData",
              data = matrix(ncol = ncol(sD@data), nrow = length(sD@coordinates)),
              libsizes = sD@libsizes,
              replicates = sD@replicates,
              coordinates = sD@coordinates,
              seglens = width(sD@coordinates))
    if(nrow(sD@data) > 0)
      lD@data <- sapply(1:ncol(sD), function(jj) as.integer(sD@data[,jj]))
    if(nrow(sD@locLikelihoods) > 0)
      lD@locLikelihoods <- sapply(1:ncol(sD@locLikelihoods), function(jj) as.double(sD@locLikelihoods[,jj]))
    lD
  }
    
.getLocLikelihoods <- function(pcD, subset, cl)
  {
    constructWeights <- function(withinCluster = FALSE)
      {
        dupReps <- which(!duplicated(numintSamp[,2]))
        dupCounts <- diff(c(dupReps, nrow(numintSamp) + 1))

        priorWeights <- matrix(0, ncol = ncol(NBpriors), nrow = nrow(NBpriors))
        if(any(dupCounts == 1))
          priorWeights[numintSamp[dupReps[dupCounts == 1],2],] <- numintSamp[dupReps[dupCounts == 1],-(1:2)]
        if(any(dupCounts > 1))
          {
            whichReps <- cbind(c(1, cumsum(dupCounts[-length(dupCounts)]) + 1), cumsum(dupCounts))
            priorWeights[numintSamp[dupReps[dupCounts > 1],2],] <- t(apply(whichReps[dupCounts > 1,], 1, function(z) colSums(numintSamp[z[1]:z[2], -c(1:2)])))
          }

        if(withinCluster) {
          assign("priorWeights", priorWeights, envir = .GlobalEnv)
          return(invisible(NULL))
        } else return(priorWeights)
      }
    
    NBdens <- function(us) {
      `logsum` <-
        function(x)
          max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
      
      
      PDgivenr.NB <- function (number, cts, seglen, libsizes, priors, weights)
        {
          seglen <- rep(seglen, length(cts))

          likes <- rowSums(
                           matrix(
                                  dnbinom(rep(cts, each = nrow(priors)),
                                          size = 1 / priors[,2],
                                          mu = rep(libsizes * seglen, each = nrow(priors)) * priors[,1], log = TRUE),
                                  ncol = length(cts))
                           )
          apply(log(weights) + likes, 2, logsum) - log(colSums(weights))                    
        }

      number <- us[1]
      us <- us[-1]
      cts <- us[-1]
      seglen <- us[1]
      
      PDgivenr.NB(number = number, cts = cts, seglen = seglen, libsizes = libsizes, priors = NBpriors, weights = priorWeights)
    }
    
    if(!is.null(cl))
      {
        clustAssign <- function(object, name)
          {
            assign(name, object, envir = .GlobalEnv)
            NULL
          }
        
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(clustAssign) <- getLikelihoodsEnv
        environment(NBdens) <- getLikelihoodsEnv
      }

    libsizes = pcD@libsizes
    NBpriors <- pcD@priors$priors
    seglens <- pcD@seglens
    numintSamp <- cbind(pcD@priors$sampled, pcD@priors$weights)
    numintSamp <- numintSamp[order(numintSamp[,2]),]

    priorWeights <- constructWeights()
    NBpriors <- NBpriors[rowSums(priorWeights) != 0,]
    priorWeights <- priorWeights[rowSums(priorWeights) != 0,]
    
    if(!is.null(cl))
      {
        clusterCall(cl, clustAssign, NBpriors, "NBpriors")
        clusterCall(cl, clustAssign, numintSamp, "numintSamp")
        clusterCall(cl, clustAssign, libsizes, "libsizes")
        clusterCall(cl, clustAssign, priorWeights, "priorWeights")
      }
    
    if (is.null(cl)) {
      ps <- t(apply(cbind(1:nrow(pcD@data), seglens, pcD@data)[subset,,drop = FALSE], 1, NBdens))
    } else {      
      ps <- parRapply(cl, cbind(1:nrow(pcD@data), seglens, pcD@data)[subset,, drop = FALSE], NBdens)
      ps <- matrix(ps, ncol = 2, byrow = TRUE)
    }

    rps <- matrix(NA, ncol = 2, nrow = nrow(pcD@data))
    rps[subset,] <- ps
    
    return(rps)
  }


.classifyNulls <- function(rep, lociPD, potNulls, lR, nullPriors, nullSampled, nullCutoff, cl) {
  
  getPosts <- function(rep, psD) {
    message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
    
    repCol <- which(levels(psD@replicates) == rep)
    repD <- psD[, psD@replicates == rep]
    repD <- .convertSegToLoci(repD)
    repD@priorType = "NB"
    replicates(repD) <- as.factor(rep(1, ncol(repD)))
    
    groups(repD) <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
        
    repD@priors$priors <- nullPriors[[repCol]]$priors
    repD@priors$sampled <- nullSampled
    repD@priors$weights <- nullPriors[[repCol]]$weights
    repD@priors$sampled[[1]][[1]][,1] <- -1
    repD@priors$sampled[[2]][[1]][,1] <- -1
    lD <- getLikelihoods.NB(cD = repD, bootStraps = 1, 
                            verbose = FALSE, returnPD = TRUE, 
                            subset = NULL, priorSubset = NULL, prs = c(0.5, 0.5), pET = "none", 
                            cl = cl)
    message("...done.", appendLF = TRUE)
    lD
  }
  
  repCol <- which(levels(lociPD@replicates) == rep)
  
  repNull <- rep(NA, nrow(potNulls))
  repLoci <- which(lociPD@locLikelihoods[, repCol])
  repLociDP <- lociPD[repLoci, ]
  postOver <- which(getOverlaps(coordinates = potNulls@coordinates, segments = repLociDP@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL))
    
  overLoci <- unique(c(postOver,
                       unlist(lapply(levels(seqnames(potNulls@coordinates)), function(chrom)
                                     which(seqnames(potNulls@coordinates) == chrom &
                                           (end(potNulls@coordinates) + 1) %in% c(start(repLociDP@coordinates[seqnames(repLociDP@coordinates) == chrom]), seqlengths(repLociDP@coordinates)[levels(seqnames(potNulls@coordinates)) == chrom]) &
                                           (start(potNulls@coordinates) - 1) %in% c(end(repLociDP@coordinates[seqnames(repLociDP@coordinates) == chrom]), 0))
                                     ))))
  
  overLoci <- overLoci[!is.na(overLoci)]
  
  if (length(postOver > 0)) {
    overNulls <- potNulls[overLoci,]
    nullsWithin <- potNulls[postOver, ]
    whichOverlaps <- 1:nrow(nullsWithin)
    
    nullsPD <- getPosts(rep, psD = overNulls)
    if(lR) {
      repNull[overLoci] <- nullsPD[,1] > nullsPD[,2]
    } else {
      
      nullFilter <- .filterSegments(overNulls@coordinates[whichOverlaps,], runif(length(whichOverlaps)))
      nullPosts <- getPosteriors(ps = nullsPD, pET = "BIC", groups = list(rep(1, sum(lociPD@replicates == rep)), rep(1, sum(lociPD@replicates == rep))), priorSubset = nullFilter, prs = c(0.1,0.9), cl = cl)
      repNull[overLoci] <- nullPosts$posteriors[,1] > log(nullCutoff)
    }
    
  }
  repNull
}


.classifyLoci <- function(potlociD, prepD, repWeights, subLoc, locSubset, lR, locps, lociCutoff, cl)
  {
    getLikeLoci <- function(rep, potlociD, priors, repWeights, sampled, ...) {
      message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
      
      repD <- new("lociData",
                  data = sapply(which(potlociD@replicates == rep), function(ii) as.integer(potlociD@data[,ii])),
                  seglens = width(potlociD@coordinates),
                  libsizes = potlociD@libsizes[potlociD@replicates == rep],
                  replicates = as.factor(rep(1, sum(potlociD@replicates == rep))),
                  coordinates = potlociD@coordinates)
      
      repD@priorType = "NB"
      groups(repD) <- list(rep(1, ncol(repD)))
      
      withinWeights <- exp(repWeights[, rep])
      
      if (all(withinWeights == 1, na.rm = TRUE)) {
        repLoci <- rep(TRUE, nrow(repD))
      } else if (all(withinWeights == 0, na.rm = TRUE)) {
        repLoci <- rep(FALSE, nrow(repD))
      } else {        
        repD@priors$priors <- priors[[1]][[which(levels(potlociD@replicates) == rep)]]
        
        repD@priors$weights <- cbind(1-(withinWeights), withinWeights)
        repD@priors$weights[is.na(repD@priors$weights)] <- 0
        repD@priors$weights <- repD@priors$weights * weightFactors
        repD@priors$sampled <- sampled
        
        subset <- intersect(locSubset, which(rowSums(repD@data) > 0))
        
        lD <- .getLocLikelihoods(repD, subset = subset, cl = cl)
        
        if(lR) repLoci <- lD[,2] > lD[,1] else {
          
          ps <- c(1-locps[rep], locps[rep])          
          rps <- t(t(lD) + log(ps))
          
          `logsum` <-
            function(x)
              max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)
          
          rps[subset,] <- rps[subset,,drop = FALSE] - apply(rps[subset,,drop = FALSE], 1, logsum)
          repLoci <- rps[,2] > log(lociCutoff)
        }
      }
      message("...done.")
      repLoci
    }

    weightFactors <- prepD@priors$weights
    priors <- prepD@priors$priors    
    sampled <- prepD@priors$sampled
    sampled[,1] <- subLoc[sampled[,1]]
    
    message("Establishing likelihoods of loci...", appendLF = TRUE)  
    locLikelihoods <- sapply(levels(potlociD@replicates), getLikeLoci, potlociD = potlociD, priors = priors, repWeights = repWeights, sampled = sampled)    
    
    locLikelihoods
  }

.constructNulls <- function(emptyNulls, sDWithin, locDef, withinOnly)
  {
    emptyNulls <- emptyNulls[strand(emptyNulls) == "*",]

    leftRight <- do.call("rbind", lapply(seqlevels(sDWithin@coordinates), function(chr) {
      left <- start(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]) -
        start(emptyNulls[seqnames(emptyNulls) == chr])[match(start(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (end(emptyNulls[seqnames(emptyNulls) == chr,]) + 1))]
      right <- end(emptyNulls[seqnames(emptyNulls) == chr])[match(end(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (start(emptyNulls[seqnames(emptyNulls) == chr,]) - 1))] - end(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),])
      cbind(left, right)
    }))

    emptyInside <- which(getOverlaps(emptyNulls, sDWithin@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL))
    
    if(!withinOnly) {            
      emptyAdjacent <- unlist(lapply(levels(seqnames(sDWithin@coordinates)), function(chr) {
        leftEmpties <- match(start(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (end(emptyNulls[seqnames(emptyNulls) == chr,]) + 1))
        rightEmpties <- match(end(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (start(emptyNulls[seqnames(emptyNulls) == chr,]) - 1))
              which(seqnames(emptyNulls) == chr)[union(leftEmpties, rightEmpties)]
      }))
      empties <- union(emptyInside, emptyAdjacent)
    } else empties <- emptyInside
    empties <- empties[!is.na(empties)]
    
    withinData <- matrix(sapply(1:ncol(sDWithin), function(ii) as.integer(sDWithin@data[,ii])), nrow = nrow(sDWithin))

    potnullD <-
      new("segData",
          data = DataFrame(rbind(matrix(0, ncol = ncol(sDWithin), nrow = length(empties)),
            withinData[!is.na(leftRight[,'left']) & !is.na(leftRight[,'right']),],
            withinData[!is.na(leftRight[,'left']),],
            withinData[!is.na(leftRight[,'right']),])),
          libsizes = sDWithin@libsizes,
          replicates = sDWithin@replicates,
          coordinates = GRanges(
            seqnames = c(seqnames(emptyNulls[empties,]), seqnames(sDWithin@coordinates)[!is.na(leftRight[,'left']) & !is.na(leftRight[,'right'])],
              (seqnames(sDWithin@coordinates))[!is.na(leftRight[,'left'])], (seqnames(sDWithin@coordinates))[!is.na(leftRight[,'right'])]),
            IRanges(
                    start = c(start(emptyNulls[empties,]), (start(sDWithin@coordinates) - leftRight[,"left"])[!is.na(leftRight[,'left']) & !is.na(leftRight[,'right'])],
                      (start(sDWithin@coordinates) - leftRight[,"left"])[!is.na(leftRight[,'left'])], (start(sDWithin@coordinates))[!is.na(leftRight[,'right'])]),
                    end = c(end(emptyNulls[empties,]), (end(sDWithin@coordinates) + leftRight[,"right"])[!is.na(leftRight[,'left']) & !is.na(leftRight[,'right'])],
                      (end(sDWithin@coordinates))[!is.na(leftRight[,'left'])], (end(sDWithin@coordinates) + leftRight[,'right'])[!is.na(leftRight[,'right'])])
                    )
            )
          )

#    potnullD@seglens <- matrix(width(potnullD@coordinates), ncol = 1)
    
    overLoci <- which(getOverlaps(coordinates = potnullD@coordinates, segments = locDef, overlapType = "within", whichOverlaps = FALSE, cl = NULL))
    
    if(!withinOnly)           
      overLoci <- union(1:length(empties), overLoci)

    potnullD <- potnullD[overLoci,]

    potnullD
  }

classifySeg <- function (sD, cD, aD, lociCutoff = 0.9, nullCutoff = 0.9, subRegion = NULL, getLikes = TRUE, lR = FALSE, samplesize = 1e5, cl, ...) 
{  
    fastUniques <- function(x){
      if (nrow(x) > 1) {
        return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x), 
                                  , drop = FALSE]) != ncol(x)))
      } else return(TRUE)
    }
    if (missing(cD)) 
        cD <- NULL
    if (is.null(cD)) 
        cD <- heuristicSeg(sD = sD, aD = aD, subRegion = subRegion, getLikes = TRUE, ..., cl = cl)

    if(lR) {
      pET <- "none"
      lociCutoff <- nullCutoff <- 0.5
    } else {
      locps <- sapply(levels(cD@replicates), function(rep) mean(exp(cD@locLikelihoods[rowSums(cD@data[,cD@replicates == rep, drop = FALSE]) > 0,rep])))
      pET <- "BIC"
    }
    
    if (is.null(subRegion)) {
      locSubset <- 1:nrow(sD)
    } else {
      locSubset <- sort(unique(c(unlist(apply(subRegion, 1, 
                                              function(sR) which(as.character(seqnames(sD@coordinates)) == as.character(sR[1]) &
                                                                 start(sD@coordinates) >= as.numeric(sR[2]) &
                                                              end(sD@coordinates) <= as.numeric(sR[3])))))))
    }

    subLoc <- getOverlaps(coordinates = sD@coordinates,
                          segments = cD@coordinates, overlapType = "within", 
                          whichOverlaps = FALSE, cl = cl)
    subLoc <- which(subLoc)[.filterSegments(sD@coordinates[subLoc,], runif(sum(subLoc)))]
    prepD <- .convertSegToLoci(sD[subLoc,])

    groups(prepD) <- list(prepD@replicates)
    prepD <- getPriors.NB(prepD, samplesize = samplesize, verbose = FALSE, cl = cl)
    weights <- cD@locLikelihoods[unlist(getOverlaps(coordinates = prepD@coordinates[prepD@priors$sampled[,1], ], segments = cD@coordinates, overlapType = "within", cl = cl)),]
    repWeights <- sapply(levels(sD@replicates), function(rep) {
      repWeights <- weights[,levels(sD@replicates) == rep]
      repWeights[rowSums(prepD@data[prepD@priors$sampled[,1], sD@replicates == rep, drop = FALSE]) == 0] <- NA
      repWeights
    })

    locLikelihoods <- .classifyLoci(sD, prepD, repWeights, subLoc, locSubset, lR, locps, lociCutoff, cl = cl)

    selLoc <- which(rowSums(locLikelihoods, na.rm = TRUE) > 0)
    if(length(selLoc) == 0)
      stop("No loci found; maybe your cutoff values are too strict?")    

    sD@locLikelihoods <- DataFrame(locLikelihoods)
    rm(locLikelihoods)
    gc()
    
    withinLoc <- which(getOverlaps(sD@coordinates, sD@coordinates[selLoc,], whichOverlaps = FALSE, overlapType = "within", cl = cl))

    if(length(withinLoc > 0))
    {            
      emptyNulls <- gaps(sD@coordinates[fastUniques(cbind(as.character(seqnames(sD@coordinates)), start(sD@coordinates)))])
      curNullsWithin <- .constructNulls(emptyNulls, sD[withinLoc,], sD@coordinates[selLoc,], withinOnly = TRUE)

      if(nrow(curNullsWithin) > 0) {      
        curNullsWithin <- curNullsWithin[.filterSegments(curNullsWithin@coordinates, runif(nrow(curNullsWithin))),]
        curNullsWithin <- .convertSegToLoci(curNullsWithin)

        message("Finding priors on 'null' regions...")
        
        curNullsWithin@groups <- list(curNullsWithin@replicates)
        curNullsWithin <- getPriors.NB(curNullsWithin, samplesize = samplesize, cl = cl)
        curNullsWithin@priors$weights <- matrix(1, nrow = nrow(curNullsWithin@priors$sampled), ncol = length(unique(curNullsWithin@replicates)))

        nullSegPriors <- cD
        nullSegPriors@locLikelihoods <- matrix(nrow = 0, ncol = 0)
        groups(nullSegPriors) <- list(replicates(nullSegPriors))
        nullSegPriors <- getPriors.NB(nullSegPriors, samplesize = samplesize, verbose = FALSE, cl = cl)       

        nullSampled = list(list(nullSegPriors@priors$sampled), list(curNullsWithin@priors$sampled))
        nullPriors <- lapply(levels(sD@replicates), function(rep) {
          repCol <- which(levels(sD@replicates) == rep)
          list(priors = list(list(nullSegPriors@priors$priors[[1]][[repCol]]), list(curNullsWithin@priors$priors[[1]][[repCol]])),
               weights = list(list((1 - exp(cD@locLikelihoods[nullSegPriors@priors$sampled[,1], repCol])) * nullSegPriors@priors$weights), list(exp(curNullsWithin@priors$weights[,repCol])))
               )
        })
        
        message("Establishing likelihoods of nulls...", appendLF = TRUE)
        potnullD <- .constructNulls(emptyNulls, sD[withinLoc,], sD@coordinates, withinOnly = FALSE)
        potnullD@locLikelihoods <- DataFrame(log(sapply(levels(potnullD@replicates), .classifyNulls, lociPD = sD[selLoc,], potNulls = potnullD, lR = lR, nullPriors = nullPriors, nullSampled = nullSampled, nullCutoff = nullCutoff, cl = cl)))
        
      } else potnullD <- new("lociData")
    } else potnullD <- new("lociData")
    
    sD@locLikelihoods <- DataFrame(sapply(1:ncol(sD@locLikelihoods), function(jj) log(sD@locLikelihoods[,jj])))

#    return(list(lociPD = lociPD, nullPD = potnullD))
    
    locMap <- .processPosteriors(sD[selLoc,], potnullD, aD = aD, lociCutoff = 1, nullCutoff = 1, getLikes = getLikes, cl = cl)
    locMap
  }
