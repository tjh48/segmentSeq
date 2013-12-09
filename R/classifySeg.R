classifySeg <- function(sD, cD, aD, lociCutoff = 0.9, nullCutoff = 0.9, subRegion = NULL, getLikes = TRUE, lR = FALSE, samplesize = 1e5, largeness = 1e8, tempDir = NULL, cl)
  {
    if(missing(aD)) stop("You must supply an aD object.")
    
    if(!is.null(subRegion))
      {
        sD <- sD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(sD@coordinates)) == subRegion$chr[ii] & end(sD@coordinates) >= subRegion$start[ii] & start(sD@coordinates) <= subRegion$end[ii]))),]
        aD <- aD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(aD@alignments)) == subRegion$chr[ii] & end(aD@alignments) >= subRegion$start[ii] & start(aD@alignments) <= subRegion$end[ii]))),]
      }
    
    if(missing(cD))
      cD <- heuristicSeg(sD, aD, largeness = largeness, getLikes = TRUE, verbose = TRUE, cl = cl)
    
    sD <- .massiveClassifyLoci(sD = sD, cD = cD, subRegion = subRegion, samplesize = samplesize, lR = lR, lociCutoff = lociCutoff, tempDir = tempDir, largeness = largeness, cl = cl)    
    
    if(!is.null(tempDir)) save(sD, file = paste(tempDir, "/sD_locLikes.RData", sep = ""))
    
    sD@locLikelihoods <- log(sD@locLikelihoods >= log(lociCutoff))
    
    nullD <- .massiveClassifyNulls(sD = sD, cD = cD, aD = aD, subRegion = subRegion, samplesize = samplesize, lR = lR, nullCutoff = nullCutoff, tempDir = tempDir, largeness = largeness, recoverOld = FALSE, cl = cl)

    if(nrow(nullD) > 0) {
      nullPD <- nullD[!values(nullD@coordinates)$empty,]
      emptyPD = nullD[values(nullD@coordinates)$empty,]
    } else nullPD <- emptyPD <- nullD
    
    lD <- .processPosteriors(lociPD = sD, nullPD = nullPD, emptyPD = emptyPD, getLikes = FALSE, cl = cl)
    
    if(getLikes & !missing(aD)) lD <- lociLikelihoods(lD, aD, cl = cl) else if(getLikes & missing(aD)) warning("I can't calculate locus likelihoods without an aD object. You can run lociLikelihoods on the output of this function for the same result.")      
    lD
  }
    

.massiveClassifyLoci <- function(sD, cD, subRegion = NULL, samplesize = 1e5, lR = FALSE, lociCutoff = 0.9, largeness = 1e8, tempDir = NULL, cl = cl)
  {
    if(!is.null(tempDir)) dir.create(tempDir, showWarnings = FALSE)
          
    if(lR) {
      pET <- "none"
      lociCutoff <- nullCutoff <- 0.5
    } else {
      locps <- sapply(levels(cD@replicates), function(rep) mean(exp(cD@locLikelihoods[rowSums(cD@data[,cD@replicates == rep, drop = FALSE]) > 0, levels(cD@replicates) == rep]), na.rm = TRUE))
      pET <- "BIC"
    }
    
    message("Finding candidate priors...", appendLF = FALSE)
    
    subLoc <- getOverlaps(coordinates = sD@coordinates,
                          segments = cD@coordinates, overlapType = "within", 
                          whichOverlaps = FALSE, cl = NULL)
    subLoc <- which(subLoc)[.filterSegments(sD@coordinates[subLoc,], runif(sum(subLoc)))]
    prepD <- .convertSegToLoci(sD[subLoc,])

    
    message("done.")
    
    groups(prepD) <- list(prepD@replicates)
    if(class(sD) == "segData") {
      prepD <- getPriors.NB(prepD, samplesize = samplesize, verbose = TRUE, cl = cl)
    } else if(class(sD) == "segMeth") prepD <- getPriors.BB(prepD, samplesize = samplesize, verbose = TRUE, cl = cl)

    if(!is.null(tempDir)) save(prepD, file = paste(tempDir, "/prepD.RData", sep = ""))
    
    weights <- cD@locLikelihoods[unlist(getOverlaps(coordinates = prepD@coordinates[prepD@priors$sampled[,1], ], segments = cD@coordinates, overlapType = "within", cl = cl)),,drop = FALSE]
    repWeights <- sapply(levels(sD@replicates), function(rep) {
      repWeights <- weights[,levels(sD@replicates) == rep]
      repWeights[rowSums(prepD@data[prepD@priors$sampled[,1], sD@replicates == rep, drop = FALSE]) == 0] <- NA
      repWeights
    })

    if(!is.null(tempDir)) save(repWeights, subLoc, prepD, file = paste(tempDir, "/prepD.RData", sep = ""))
    
    sDsplit <- .splitSD(sD, largeness)

    message("Segmentation split into ", length(sDsplit), " parts.")
    
    sD@locLikelihoods <- (do.call("rbind", lapply(1:length(sDsplit), function(ii) {  
      message("Establishing likelihoods of loci; Part ", ii, " of ", length(sDsplit))
      x <- sD[sDsplit[[ii]],]
      if (is.null(subRegion)) {
        locSubset <- 1:nrow(x)
      } else {
        locSubset <- sort(unique(c(unlist(apply(subRegion, 1, 
                                                function(sR) which(as.character(seqnames(x@coordinates)) == as.character(sR[1]) &
                                                                   start(x@coordinates) >= as.numeric(sR[2]) &
                                                                   end(x@coordinates) <= as.numeric(sR[3])))))))
      }
      subLoc <- .classifyLoci(potlociD = x, prepD = prepD, repWeights = repWeights, subLoc = subLoc, locSubset = locSubset, lR = lR, locps = locps, cl = cl)

      if(!is.null(tempDir)) save(subLoc, file = paste(tempDir, "/subLocLikes_", ii, ".RData", sep = ""))
      subLoc
    })))
    
    sD
  }


.massiveClassifyNulls <- function(sD, cD, aD, subRegion = NULL, samplesize = 1e5, lR = FALSE, nullCutoff = 0.9, largeness = 1e8, tempDir, recoverOld = FALSE, cl = cl)
  {
    #make nullPriors - again, this copies classifySeg and should be split off as a separate function...    
    selLoc <- which(rowSums(sapply(1:ncol(sD@locLikelihoods), function(jj) sD@locLikelihoods[,jj] == 0), na.rm = TRUE) > 0)
    if(length(selLoc) == 0)
      stop("No loci found; maybe the cutoff values used to generate the sD@locLikelihoods are too strict?")

    withinLoc <- which(getOverlaps(sD@coordinates, sD@coordinates[selLoc,], whichOverlaps = FALSE, overlapType = "within", cl = NULL))

    if(length(withinLoc > 0))
    {
      emptyNulls <- gaps(sD@coordinates[.fastUniques(cbind(as.character(seqnames(sD@coordinates)), start(sD@coordinates), as.character(strand(sD@coordinates))))])
      
      curNullsWithin <- .constructNulls(emptyNulls, sD[withinLoc,], sD@coordinates[selLoc,], forPriors = TRUE, samplesize = samplesize, aD = aD)
      
      gc()

      if(nrow(curNullsWithin) > 0) {
        if(!recoverOld) {        
          message("Finding priors on 'null' regions...")
          
          curNullsWithin <- .convertSegToLoci(curNullsWithin)        
          curNullsWithin@groups <- list(curNullsWithin@replicates)
          if(class(curNullsWithin) == "lociData")
            {
              curNullsWithin <- getPriors.NB(curNullsWithin, samplesize = samplesize, verbose = FALSE, cl = cl)
            } else if(class(curNullsWithin) == "methData")
              curNullsWithin <- getPriors.BB(curNullsWithin, samplesize = samplesize, verbose = FALSE, cl = cl)
          
          curNullsWithin@priors$weights <- matrix(1, nrow = nrow(curNullsWithin@priors$sampled), ncol = length(unique(curNullsWithin@replicates)))
          
          nullSegPriors <- cD
          nullSegPriors@locLikelihoods <- matrix(nrow = 0, ncol = 0)
          groups(nullSegPriors) <- list(replicates(nullSegPriors))
          if(class(nullSegPriors) == "lociData") {
            nullSegPriors <- getPriors.NB(nullSegPriors, samplesize = samplesize, verbose = FALSE, cl = cl)
          } else if(class(nullSegPriors) == "methData")
            nullSegPriors <- getPriors.BB(nullSegPriors, samplesize = samplesize, verbose = FALSE, cl = cl)
          
          nullSampled = list(list(nullSegPriors@priors$sampled), list(curNullsWithin@priors$sampled))        
          nullPriors <- lapply(levels(sD@replicates), function(rep) {
            repCol <- which(levels(sD@replicates) == rep)
            list(priors = list(list(nullSegPriors@priors$priors[[1]][[repCol]]), list(curNullsWithin@priors$priors[[1]][[repCol]])),
                 weights = list(list((1 - exp(cD@locLikelihoods[nullSegPriors@priors$sampled[,1], repCol])) * nullSegPriors@priors$weights), list(exp(curNullsWithin@priors$weights[,repCol])))
                 )                
          })          
          if(!is.null(tempDir)) save(nullPriors, nullSampled, file = paste(tempDir, "/nullPriors.RData", sep = ""))
        } else if(!is.null(tempDir)) load(file = paste(tempDir, "/nullPriors.RData", sep = ""))
        
        splitLoci <- .splitSD(sD, largeness)        
        if(!is.null(tempDir)) save(emptyNulls, splitLoci, file = paste(tempDir, "/nullSplit.RData", sep = ""))
        
        message("Segmentation split into ", length(splitLoci), " parts.")
        
        splitNulls <- lapply(1:length(splitLoci), function(ii) {

          message("Establishing likelihoods of nulls; Part ", ii, " of ", length(splitLoci))
          subSD <- sD[splitLoci[[ii]],]
          subWithinLoc <- getOverlaps(subSD@coordinates, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,],
                                      whichOverlaps = FALSE, overlapType = "within", cl = NULL)
          if(any(subWithinLoc))
            {              
#              if(class(sD) == "segData") {
#                potnullD <- .constructNulls(emptyNulls, subSD[subWithinLoc,], subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,])
#              } else if(class(sD) == "segMeth") {
#                sDPSmall <- .fastUniques(data.frame(chr = as.character(seqnames(subSD@coordinates)), start = as.numeric(start(subSD@coordinates))))
#                empties <- .zeroInMeth(aD = aD, smallSegs = subSD@coordinates[sDPSmall])
#                potnullD <- .constructMethNulls(emptyNulls, subSD[subWithinLoc,], subSD[rowSums(subSD@locLikelihoods == 0) > 0,])
                
#                counts <- getCounts(potnullD@coordinates, aD, cl = cl)        
#                potnullD@Cs <- counts$Cs
#                potnullD@Ts <- counts$Ts                
#              }

              potnullD <- .constructNulls(emptyNulls, subSD[subWithinLoc,], subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,], aD = aD)
              
              potnullD@locLikelihoods <- log(do.call("cbind", ((lapply(levels(potnullD@replicates), .classifyNulls, lociPD = subSD, potNulls = potnullD, lR = lR, nullPriors = nullPriors, nullSampled = nullSampled, nullCutoff = nullCutoff, cl = cl)))))
            } else potnullD <- NULL
          if(class(potnullD) == "segData") {
            emptyD <- rowSums(potnullD@data) == 0
          } else if(class(potnullD) == "segMeth") {
            emptyD <- rowSums(potnullD@Cs) == 0
          }
          
          values(potnullD@coordinates)$empty <- FALSE
          values(potnullD@coordinates)$empty[emptyD] <- TRUE
          potnullD <- potnullD[rowSums(sapply(1:ncol(potnullD@locLikelihoods), function(ii) potnullD@locLikelihoods[,ii] > -Inf), na.rm = TRUE) > 0 | emptyD,]
          if(class(potnullD) == "segData") {
            potnullD@data <- (matrix(ncol = length(potnullD@replicates), nrow = 0))
          } else if(class(potnullD) == "segMeth") {
            potnullD@Cs <- potnullD@Ts <- (matrix(ncol = length(potnullD@replicates), nrow = 0))
          }
          if(!is.null(tempDir)) save(potnullD, file = paste(tempDir, "/subNull_", ii, ".RData", sep = ""))
          
          potnullD
        })

        if(any(!sapply(splitNulls, is.null)))
          potnullD <- .mergeSegData(splitNulls[!sapply(splitNulls, is.null)]) else potnullD <- new(class(sD))
        

      } else potnullD <- new(class(sD))
    } else potnullD <- new(class(sD))

    potnullD
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
      
      PDgivenr.NB <- function (number, cts, seglen, libsizes, priors, weights)
        {
          `logsum` <-
            function(x)
              max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)          
          
          seglen <- rep(seglen, length(cts))
          likes <- rowSums(
                           matrix(
                                  dnbinom(rep(cts, each = nrow(priors)),
                                          size = 1 / priors[,2],
                                          mu = rep(as.numeric(libsizes) * as.numeric(seglen), each = nrow(priors)) * priors[,1], log = TRUE),
                                  ncol = length(cts))
                           )

          c(logsum(log(weights[,1]) + likes), logsum(log(weights[,2]) + likes)) - log(colSums(weights))
        }
      
      number <- us[1]
      us <- us[-1]
      cts <- us[-1]
      seglen <- us[1]
      
      pd <- PDgivenr.NB(number = number, cts = cts, seglen = seglen, libsizes = libsizes, priors = NBpriors, weights = priorWeights)
      pd
    }
    
    BBdens <- function(us) {
      
      PDgivenr.BB <- function (number, cts, secondCts, priors, weights)
        {
          `logsum` <-
            function(x)
              max(x, max(x, na.rm = TRUE) + log(sum(exp(x - max(x, na.rm = TRUE)), na.rm = TRUE)), na.rm = TRUE)

          dbetabinom <- function(x, n, prop, disp) {
              
              smallDisp <- disp < 1e-15 & disp >= 0
              largeDisp <- disp > 1e-15
              
              ps <- matrix(NA, ncol = ncol(prop), nrow = nrow(prop))
              disp <- matrix(disp, ncol = ncol(prop), nrow = nrow(prop))
              x <- matrix(x, ncol = ncol(prop), nrow = nrow(prop), byrow = TRUE)
              n <- matrix(n, ncol = ncol(prop), nrow = nrow(prop), byrow = TRUE)
              
              if(any(largeDisp)) {

                alpha = (1/disp - 1) * prop
                beta = (1/disp - 1) * (1-prop)
                
                ps[largeDisp,] <- lchoose(n[largeDisp,,drop = FALSE], x[largeDisp,,drop = FALSE]) +
                  lbeta((x + alpha)[largeDisp,,drop = FALSE], (n - x + beta)[largeDisp,, drop = FALSE]) -
                    lbeta(alpha[largeDisp,, drop = FALSE], beta[largeDisp,, drop = FALSE]) 
              }
              if(any(smallDisp))
                ps[smallDisp,] <- dbinom(x[smallDisp,,drop = FALSE], n[smallDisp,,drop = FALSE], prob = prop[smallDisp,,drop = FALSE], log = TRUE)
              
              return(ps)
            }

          prop <- matrix(priors[,1], ncol = length(cts), nrow = nrow(priors))          
          
          likes <- rowSums(
                           dbetabinom(cts,
                                      cts + secondCts,
                                      disp = priors[,2],
                                      prop = prop))
          
          c(logsum(log(weights[,1]) + likes), logsum(log(weights[,2]) + likes)) - log(colSums(weights))
        }
      
      number <- us[1]
      cts <- us[2:((length(us) - 1) / 2 + 1)]
      secondCts <- us[((length(us) - 1) / 2 + 2):length(us)]
      
      pd <- PDgivenr.BB(number = number, cts = cts, secondCts = secondCts, priors = NBpriors, weights = priorWeights)
      pd
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

    if(class(pcD) == "lociData")
      {
        if (is.null(cl)) {
          ps <- t(apply(cbind(1:nrow(pcD@data), seglens, pcD@data)[subset,,drop = FALSE], 1, NBdens))
        } else {      
          ps <- parRapply(cl, cbind(1:nrow(pcD@data), seglens, pcD@data)[subset,, drop = FALSE], NBdens)
          ps <- matrix(ps, ncol = 2, byrow = TRUE)
        }
      } else if(class(pcD) == "methData") {
        if (is.null(cl)) {
          ps <- t(apply(cbind(1:nrow(pcD@data), pcD@data, pcD@pairData)[subset,,drop = FALSE], 1, BBdens))
        } else {      
          ps <- parRapply(cl, cbind(1:nrow(pcD@data), pcD@data, pcD@pairData)[subset,, drop = FALSE], BBdens)
          ps <- matrix(ps, ncol = 2, byrow = TRUE)
        }
      }

    rps <- matrix(NA, ncol = 2, nrow = nrow(pcD@data))
    rps[subset,] <- ps
    
    return(rps)
  }


.classifyNulls <- function(rep, lociPD, potNulls, lR, nullPriors, nullSampled, nullCutoff, cl) {
  
  message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
  repCol <- which(levels(potNulls@replicates) == rep)
  repD <- potNulls[, potNulls@replicates == rep]

  repD <- .convertSegToLoci(repD)

  if(class(repD) == "lociData") repD@priorType = "NB" else if(class(repD) == "methData") repD@priorType = "BB"    
    
  replicates(repD) <- as.factor(rep(1, ncol(repD)))
  groups(repD) <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
  
  repD@priors$priors <- nullPriors[[repCol]]$priors
  repD@priors$sampled <- nullSampled
  repD@priors$weights <- nullPriors[[repCol]]$weights
  repD@priors$sampled[[1]][[1]][,1] <- -1
  repD@priors$sampled[[2]][[1]][,1] <- -1
  
  repCol <- which(levels(lociPD@replicates) == rep)
  repNull <- rep(NA, nrow(repD))
  repLoci <- which(lociPD@locLikelihoods[, repCol] == 0)
  repLociDP <- lociPD[repLoci, ]

  postOver <- which(getOverlaps(coordinates = repD@coordinates, segments = repLociDP@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL))
    
  overLoci <- unique(c(postOver,
                       unlist(lapply(levels(seqnames(repD@coordinates)), function(chrom)
                                     which(seqnames(repD@coordinates) == chrom &
                                           (end(repD@coordinates) + 1) %in% c(start(repLociDP@coordinates[seqnames(repLociDP@coordinates) == chrom]), seqlengths(repLociDP@coordinates)[levels(seqnames(repD@coordinates)) == chrom]) &
                                           (start(repD@coordinates) - 1) %in% c(end(repLociDP@coordinates[seqnames(repLociDP@coordinates) == chrom]), 0))
                                     ))))
  
  overLoci <- sort(overLoci[!is.na(overLoci)])
  
  if (length(postOver) > 0) {
    overNulls <- repD[overLoci,]
    nullsWithin <- repD[postOver, ]
    lD <- rep(NA, nrow(overNulls))
    
    if(lR) {
      prs <- c(0.5, 0.5)
      nullCutoff <- 0.5
    } else {
      nullFilter <- .filterSegments(overNulls@coordinates, runif(nrow(overNulls)))              
      if(class(overNulls) == "lociData") {
        fD <- getLikelihoods.NB(cD = overNulls[nullFilter,], bootStraps = 1, 
                                verbose = FALSE, 
                                subset = NULL, priorSubset = NULL, pET = "BIC", 
                                cl = cl)
      } else if(class(overNulls) == "methData") {
        fD <- getLikelihoods.BB(cD = overNulls[nullFilter,], bootStraps = 1, 
                                verbose = FALSE, 
                                subset = NULL, priorSubset = NULL, pET = "BIC", 
                                cl = cl)
      }
      lD[nullFilter] <- fD@posteriors[,1] > log(nullCutoff)
      prs <- fD@estProps      
    }    

    fillInLD <- function(lD, selMin)
      {
        if(any(is.na(lD)))
          {
            if(selMin) subSel <- .fastUniques(cbind(as.character(seqnames(overNulls@coordinates[is.na(lD)])), start(overNulls@coordinates[is.na(lD)]))) else subSel <- TRUE

            nullFilter <- which(is.na(lD))[subSel]

            if(class(overNulls) == "lociData") {
              fD <- getLikelihoods.NB(cD = overNulls[nullFilter,], bootStraps = 1, 
                                      verbose = FALSE, 
                                      subset = NULL, priorSubset = NULL, pET = "none", prs = c(prs[1], 1 - prs[1]), 
                                      cl = cl)
            } else if(class(overNulls) == "methData") {
              fD <- getLikelihoods.BB(cD = overNulls[nullFilter,], bootStraps = 1, 
                                      verbose = FALSE, 
                                      subset = NULL, priorSubset = NULL, pET = "none", prs = c(prs[1], 1 - prs[1]), 
                                      cl = cl)
            }
                            
            lD[nullFilter] <- fD@posteriors[,1] > log(nullCutoff)
          }
        lD
      }

    lD <- fillInLD(lD, selMin = TRUE)

    if(any(lD, na.rm = TRUE))
      lD[getOverlaps(overNulls@coordinates, overNulls@coordinates[which(lD),], overlapType = "contains", whichOverlaps = FALSE, cl = NULL)] <- TRUE

    lD <- fillInLD(lD, selMin = FALSE)

    message("...done.", appendLF = TRUE)

    repNull[overLoci] <- lD
  }
  repNull  
}


.classifyLoci <- function(potlociD, prepD, repWeights, subLoc, locSubset, lR, locps, cl)
  {
    getLikeLoci <- function(rep, potlociD, priors, repWeights, sampled, ...) {
      message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
      
      if(class(potlociD) == "segData")
        {
          repD <- new("lociData",
                      data = sapply(which(potlociD@replicates == rep), function(ii) as.integer(potlociD@data[,ii])),
                      seglens = width(potlociD@coordinates),
                      libsizes = potlociD@libsizes[potlociD@replicates == rep],
                      replicates = as.factor(rep(1, sum(potlociD@replicates == rep))),
                      coordinates = potlociD@coordinates,
                      groups = list(rep(1, sum(potlociD@replicates == rep))))
          repD@priorType = "NB"
        } else if(class(potlociD) == "segMeth") {
          repD <- new("methData",
                      data = round(potlociD@Cs[,which(potlociD@replicates == rep), drop = FALSE]),
                      pairData = round(potlociD@Ts[,which(potlociD@replicates == rep), drop = FALSE]),
                      replicates = as.factor(rep(1, sum(potlociD@replicates == rep))),
                      coordinates = potlociD@coordinates,
                      groups = list(rep(1, sum(potlociD@replicates == rep))),
                      libsizes = potlociD@nonconversion[potlociD@replicates == rep] + 1,
                      pairLibsizes = 1 - potlociD@nonconversion[potlociD@replicates == rep])
          repD@priorType = "BB"
        }     
      
      withinWeights <- exp(repWeights[, levels(potlociD@replicates) == rep])
      
      if (all(withinWeights == 1, na.rm = TRUE)) {
        warning(paste("I didn't find any non-loci to establish weights for the", rep, "replicate group."))
        repLoci <- rep(TRUE, nrow(repD))
      } else if (all(withinWeights == 0, na.rm = TRUE)) {
        warning(paste("I didn't find any loci to establish weights for the", rep, "replicate group."))
        repLoci <- rep(FALSE, nrow(repD))
      } else {        
        repD@priors$priors <- priors[[1]][[which(levels(potlociD@replicates) == rep)]]
        
        repD@priors$weights <- cbind(1-(withinWeights), withinWeights)
        repD@priors$weights[is.na(repD@priors$weights)] <- 0
        repD@priors$weights <- repD@priors$weights * weightFactors
        repD@priors$sampled <- sampled
        
        subset <- intersect(locSubset, which(rowSums(repD@data) > 0))        
        lD <- matrix(NA, ncol = 2, nrow = nrow(repD))
        
        if(length(subset) > 0) {
          if(class(potlociD) == "segData")
            {            
              orddat <- do.call("order", c(lapply(1:ncol(repD@seglens), function(ii) repD@seglens[,ii]), lapply(1:ncol(repD), function(ii) repD@data[,ii])))        
              whunq <- .fastUniques(cbind(repD@seglens, repD@data)[orddat,])
              lD <- .getLocLikelihoods(repD, subset = intersect(orddat[whunq], subset), cl = cl)
              lD[orddat,] <- lD[orddat[rep(which(whunq), diff(c(which(whunq), length(whunq) + 1)))],]
            } else if(class(potlociD) == "segMeth") {
              orddat <- do.call("order", c(lapply(1:ncol(repD@data), function(ii) repD@data[,ii]), lapply(1:ncol(repD@pairData), function(ii) repD@pairData[,ii])))
              whunq <- .fastUniques(cbind(repD@data, repD@pairData)[orddat,])
              lD <- .getLocLikelihoods(repD, subset = intersect(orddat[whunq], subset), cl = cl)
              lD[orddat,] <- lD[orddat[rep(which(whunq), diff(c(which(whunq), length(whunq) + 1)))],]
            }
        }
        
        if(lR) repLoci <- lD[,2] > lD[,1] else {
          
          ps <- c(1-locps[levels(potlociD@replicates) == rep], locps[levels(potlociD@replicates) == rep])
          rps <- t(t(lD) + log(ps))

          logRowSum <- function(z)
            {
              maxes <- do.call(pmax, c(as.list(data.frame(z)), list(na.rm = TRUE)))
              pmax(maxes, maxes + log(rowSums(exp(z - maxes), na.rm = TRUE)), na.rm = TRUE)
            }

          rps <- rps - logRowSum(rps)
          repLoci <- rps[,2]
          repLoci[intersect(locSubset, which(rowSums(repD@data) == 0))] <- -Inf
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
