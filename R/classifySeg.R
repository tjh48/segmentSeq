.getLocLikelihoods <- function(pcD, subset, cl = cl)
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

#          wsInfo <- which(numintSamp[,1] == number)
#          weights[numintSamp[wsInfo,2],] <- weights[numintSamp[wsInfo,2],] - numintSamp[wsInfo,-(1:2)]

#          likes <- rowSums(
#                           matrix(
#                                  dnbinom(rep(cts, each = (nzWts)),
#                                          size = 1 / priors[nzWts,2],
#                                          mu = rep(libsizes * seglen, each = sum(nzWts)) * priors[nzWts,1], log = TRUE),
#                                  ncol = length(cts))
#                           )
#
#          apply(log(weights[nzWts,]) + likes, 2, logsum) - log(colSums(weights[nzWts,]))

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
#        environment(constructWeights) <- getLikelihoodsEnv
#        clusterCall(cl, constructWeights, TRUE)
      }
    
#    message("...", appendLF = FALSE)
    
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
        cD <- heuristicSeg(sD = sD, aD = aD, subRegion = subRegion, ..., cl = cl)
#    if (class(cD) == "countData") 
#        cD <- lociLikelihoods(cD = cD, aD = aD, cl = cl)

    if(lR) {
      pET <- "none"
      lociCutoff <- nullCutoff <- 0.5
    } else {
      locps <- sapply(unique(cD@replicates), function(rep) mean(exp(cD@posteriors[rowSums(cD@data[,cD@replicates == rep, drop = FALSE]) > 0,rep])))
      pET <- "BIC"
    }

    
    smallLoci <- which(fastUniques(cbind(sD@segInfo$chr, sD@segInfo$start)))

    potlociD <- with(sD@segInfo, new("postSeg",
                                     data = sD@data,
                                     seglens = end - start + 1,
                                     libsizes = sD@libsizes,
                                     replicates = sD@replicates,
                                     annotation = data.frame(chr = chr, start = start, end = end, segType = "potloc")))

    if (is.null(subRegion)) {
      locSubset <- 1:nrow(potlociD)
    } else {
      locSubset <- sort(unique(c(unlist(apply(subRegion, 1, 
                                           function(sR) which(as.character(potlociD@annotation$chr) == as.character(sR[1]) &
                                                              potlociD@annotation$start >= as.numeric(sR[2]) &
                                                              potlociD@annotation$end <= as.numeric(sR[3])))))))
    }

    subLoc <- getOverlaps(coordinates = potlociD@annotation,
                          segments = cD@annotation, overlapType = "within", 
                          whichOverlaps = FALSE, cl = cl)
    subLoc <- which(subLoc)[filterSegments(potlociD@annotation[subLoc,], runif(sum(subLoc)))]
    prepD <- potlociD[subLoc,]
    prepD@groups <- list(prepD@replicates)
    prepD <- getPriors.NB(prepD, samplesize = samplesize, verbose = FALSE, cl = cl)
    weights <- cD@posteriors[unlist(getOverlaps(coordinates = prepD@annotation[prepD@priors$sampled[,1], ], segments = cD@annotation, overlapType = "within", cl = cl)),]
  repWeights <- sapply(unique(potlociD@replicates), function(rep) {
    repWeights <- weights[,rep]
    repWeights[rowSums(prepD@data[prepD@priors$sampled[,1], potlociD@replicates == rep, drop = FALSE]) == 0] <- NA
    repWeights
  })

    weightFactors <- prepD@priors$weights
    priors <- prepD@priors$priors    
    sampled <- prepD@priors$sampled
    sampled[,1] <- subLoc[sampled[,1]]
    
    getLikeLoci <- function(rep, psD, priors, repWeights, sampled, smallLoci, ...) {
      message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
      
      repD <- psD[, psD@replicates == rep]
      repD@priorType = "NB"
      repD@replicates <- rep(1, ncol(repD))
      repD@groups <- list(rep(1, ncol(repD)))#, rep(1, ncol(repD)))

#      subset <- which(rowSums(repD@data) > 0)
#      subLoc <- getOverlaps(coordinates = repD@annotation[subset,],
#                            segments = cD@annotation, overlapType = "within", 
#                            whichOverlaps = FALSE, cl = cl)
#      sampSubset <- (subset[subLoc])[filterSegments(repD@annotation[subset[subLoc],], runif(sum(subLoc)))]

#      repWeights <- cD@posteriors[unlist(getOverlaps(coordinates = repD@annotation[sampSubset,],
#                                                     segments = cD@annotation, overlapType = "within", 
#                                                     whichOverlaps = TRUE, cl = cl)),rep]
      
#      repD <- getPriors.NB(repD, verbose = FALSE, cl = cl, samplingSubset = sampSubset)
#      withinWeights <- cD@posteriors[unlist(getOverlaps(coordinates = repD@annotation[repD@priors$sampled[,1], ], segments = cD@annotation, overlapType = "within", cl = cl)),rep]


      withinWeights <- exp(repWeights[, rep])

      if (all(withinWeights == 1, na.rm = TRUE)) {
        repLoci <- rep(TRUE, nrow(repD))
      } else if (all(withinWeights == 0, na.rm = TRUE)) {
        repLoci <- rep(FALSE, nrow(repD))
      } else {        
#        repD@priors$priors <- repD@priors$priors[[1]][[1]]
        repD@priors$priors <- priors[[1]][[rep]]

        repD@priors$weights <- cbind(1-(withinWeights), withinWeights)
        repD@priors$weights[is.na(repD@priors$weights)] <- 0
        repD@priors$weights <- repD@priors$weights * weightFactors
        repD@priors$sampled <- sampled
                
        subset <- intersect(locSubset, which(rowSums(repD@data) > 0))

#        subset <- c(1:1000, 8261628)
        
        lD <- .getLocLikelihoods(repD, subset = subset, cl = cl)

        if(lR) repLoci <- lD[,2] > lD[,1] else {
          
          #repSmall <- smallLoci[rowSums(repD@data[smallLoci,,drop = FALSE]) > 0]
          #whichPriorSubset <- unlist(lapply(unique(repD@annotation$chr), function(chr) {
          #  chrrep <- repD@annotation$chr == chr
          #  which(chrrep & repD@annotation$start %in% repD@annotation$start[intersect(which(chrrep), repSmall)] & repD@annotation$end %in% repD@annotation$end[intersect(which(chrrep), repSmall)])
          #}))                      
          #priorSubset <- whichPriorSubset[filterSegments(repD@annotation[whichPriorSubset, ], runif(length(whichPriorSubset)))]
                  
          #priorLoc <- sum(lD[priorSubset,2] > lD[priorSubset,1], na.rm = TRUE) / length(priorSubset)
          #ps <- c(1-priorLoc, priorLoc)
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
  
  message("Establishing likelihoods of loci...", appendLF = TRUE)  
  potlociD@posteriors <- sapply(unique(potlociD@replicates), getLikeLoci, psD = potlociD, priors = priors, repWeights = repWeights, sampled = sampled, smallLoci = smallLoci)
  lociPD <- potlociD[which(rowSums(potlociD@posteriors, na.rm = TRUE) > 0),]

  withinLoc <- which(getOverlaps(sD@segInfo, lociPD@annotation, whichOverlaps = FALSE, overlapType = "within", cl = cl))

  if(length(withinLoc > 0))
    {      
      constructNulls <- function(withinLoc, locDef, withinOnly)
        {
          smallIn <- intersect(smallLoci, withinLoc)
          emptyNulls <- with(sD@segInfo, data.frame(chr = chr[smallIn], start = start[smallIn] - leftSpace[smallIn], end = start[smallIn] - 1))
          emptyNulls <- emptyNulls[emptyNulls$start > 1, ]
          emptyNulls <- emptyNulls[emptyNulls$end - emptyNulls$start + 1 > 0,]
          
          sDWithin <- sD[withinLoc,]
          
          potnullD <- with(sDWithin@segInfo, new("postSeg",
                                             data = rbind(matrix(0, ncol = ncol(sDWithin), nrow = nrow(emptyNulls)),
                                               sDWithin@data[leftSpace > 0 & rightSpace > 0,],
                                               sDWithin@data[leftSpace > 0,],
                                               sDWithin@data[rightSpace > 0,]),
                                             seglens = c(emptyNulls$end - emptyNulls$start + 1,
                                               (end - start + 1 + leftSpace + rightSpace)[leftSpace > 0 & rightSpace > 0],
                                               (end - start + 1 + leftSpace)[leftSpace > 0],
                                               (end - start + 1 + rightSpace)[rightSpace > 0]),
                                                 libsizes = sDWithin@libsizes,
                                                 replicates = sDWithin@replicates,   
                                                 annotation = data.frame(chr = (sDWithin@segInfo$chr)[c(emptyNulls$chr, sDWithin@segInfo$chr[sDWithin@segInfo$leftSpace > 0 & sDWithin@segInfo$rightSpace > 0],
                                                                           sDWithin@segInfo$chr[sDWithin@segInfo$leftSpace > 0], sDWithin@segInfo$chr[sDWithin@segInfo$rightSpace > 0])],
                                                   start = c(emptyNulls$start, (start - leftSpace)[leftSpace > 0 & rightSpace > 0], (start - leftSpace)[leftSpace > 0], start[rightSpace > 0]),
                                                   end = c(emptyNulls$end, (end + rightSpace)[leftSpace > 0 & rightSpace > 0], end[leftSpace > 0], (end + rightSpace)[rightSpace > 0]),
                                                   nullClass = as.factor(rep(c("empty", "expBoth", "expLeft", "expRight"), c(nrow(emptyNulls), sum(leftSpace > 0 & rightSpace > 0), sum(leftSpace > 0), sum(rightSpace > 0)))))
                                                 ))

          overLoci <- which(getOverlaps(coordinates = potnullD@annotation, segments = locDef, overlapType = "within", whichOverlaps = FALSE, cl = cl))
          
          if(!withinOnly)           
            overLoci <- unique(c(overLoci,
                                 unlist(lapply(unique(potnullD@annotation$chr), function(chrom)
                                               which(potnullD@annotation$chr == chrom)[with(potnullD@annotation[potnullD@annotation$chr == chrom,],
                                                       which(nullClass == "empty" &
                                                             (end + 1) %in% c(locDef$start[locDef$chr == chrom], aD@chrs$len[aD@chrs$chr == chrom]) &
                                                             (start - 1) %in% c(locDef$end[locDef$chr == chrom], 0)))]
                                               ))))
          
          potnullD <- potnullD[overLoci,]
          potnullD
        }

      curNullsWithin <- constructNulls(withinLoc, lociPD@annotation, withinOnly = TRUE)
      curNullsWithin <- curNullsWithin[filterSegments(curNullsWithin@annotation, runif(nrow(curNullsWithin))),]
      curNullsWithin@groups <- list(curNullsWithin@replicates)
      curNullsWithin <- getPriors.NB(curNullsWithin, samplesize = samplesize, cl = cl)
      curNullsWithinWeights <- matrix(1, nrow = nrow(curNullsWithin@priors$sampled), ncol = length(unique(curNullsWithin@replicates)))#cD@posteriors[unlist(getOverlaps(curNullsWithin@annotation[curNullsWithin@priors$sampled[,1],], cD@annotation, overlapType = "within", cl = cl)),]

      potnullD <- constructNulls(withinLoc, lociPD@annotation, withinOnly = FALSE)
      
      
      getNullPosteriors <- function(rep, lociPD, potNulls) {
        
        getPosts <- function(rep, psD) {
          message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
          
          repD <- psD[, psD@replicates == rep]
          repD@priorType = "NB"
          repD@replicates <- rep(1, ncol(repD))
          
          repD@groups <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
          
          nullPriors <- cD
          nullPriors@posteriors <- matrix(nrow = 0, ncol = 0)
          nullPriors <- nullPriors[,nullPriors@replicates == rep]
          nullPriors@replicates <- rep(1, ncol(nullPriors))
          nullPriors <- getPriors.NB(nullPriors, samplesize = samplesize, verbose = FALSE, cl = cl)
          
          repD@priors$priors <- list(list(nullPriors@priors$priors[[1]][[1]]), list(curNullsWithin@priors$priors[[1]][[rep]]))
          repD@priors$sampled <- list(list(nullPriors@priors$sampled), list(curNullsWithin@priors$sampled))
          repD@priors$weights <- list(list((1 - exp(cD@posteriors[nullPriors@priors$sampled[,1], rep])) * nullPriors@priors$weights), list(exp(curNullsWithinWeights[,rep])))
          repD@priors$sampled[[1]][[1]][,1] <- -1
          repD@priors$sampled[[2]][[1]][,1] <- -1
          lD <- getLikelihoods.NB(cD = repD, bootStraps = 1, 
                                  verbose = FALSE, returnPD = TRUE, 
                                  subset = NULL, priorSubset = NULL, prs = c(0.5, 0.5), pET = "none", 
                                  cl = cl)
          message("...done.", appendLF = TRUE)
          lD
        }

        repNull <- rep(NA, nrow(potNulls))
        repLoci <- which(lociPD@posteriors[, rep])
        repLociDP <- lociPD[repLoci, ]
        postOver <- which(getOverlaps(coordinates = potNulls@annotation, segments = repLociDP@annotation, overlapType = "within", whichOverlaps = FALSE, cl = cl))
        overLoci <- unique(c(postOver,
                             unlist(lapply(unique(potNulls@annotation$chr), function(chrom)
                                           which(potNulls@annotation$chr == chrom)[with(potNulls@annotation[potNulls@annotation$chr == chrom,],
                                                   which(nullClass == "empty" &
                                                         (end + 1) %in% c(repLociDP@annotation$start[repLociDP@annotation$chr == chrom], aD@chrs$len[aD@chrs$chr == chrom]) &
                                                         (start - 1) %in% c(repLociDP@annotation$end[repLociDP@annotation$chr == chrom], 0)))]
                                           ))
                             ))
        overLoci <- overLoci[!is.na(overLoci)]
        
        if (length(postOver > 0)) {
          overNulls <- potNulls[overLoci,]
          nullsWithin <- potNulls[postOver, ]
          whichOverlaps <- 1:nrow(nullsWithin)
          nullsPD <- getPosts(rep, psD = overNulls)
          if(lR) {
            repNull[overLoci] <- nullsPD[,1] > nullsPD[,2]
          } else {

            nullFilter <- filterSegments(overNulls@annotation[whichOverlaps,], runif(length(whichOverlaps)))
            nullPosts <- getPosteriors(ps = nullsPD, pET = "BIC", groups = list(rep(1, sum(lociPD@replicates == rep)), rep(1, sum(lociPD@replicates == rep))), priorSubset = nullFilter, prs = c(0.1,0.9), cl = cl)
            repNull[overLoci] <- nullPosts$posteriors[,1] > log(nullCutoff)
          }
            
#            message(" Refining priors...", appendLF = FALSE)
#            for (ii in 1:100) {
#
                                        #              nullFilter <- filterSegments(overNulls@annotation[whichOverlaps,], runif(length(whichOverlaps)))
#              prs <- sum(nullsPD[nullFilter,1] < nullsPD[nullFilter,2]) / length(nullFilter)
#              prs <- c(1 - prs, prs)
#              
                                        #              nullPosts <- getPosteriors(ps = nullsPD, pET = "none", groups = list(rep(1, sum(lociPD@replicates == rep)), rep(1, sum(lociPD@replicates == rep))), priorSubset = NULL, prs = prs, cl = cl)
                                        #              nullsDP <- overNulls[which(nullPosts$posteriors[,1] > log(nullCutoff)), ]
#              if (nrow(nullsDP) > 0) {
#                newLoci <- which(!getOverlaps(coordinates = lociPD@annotation, segments = nullsDP@annotation, overlapType = "contains", whichOverlaps = FALSE, cl = NULL) & lociPD@posteriors[, rep])
#                if ((all(newLoci %in% repLoci) & all(repLoci %in% newLoci)))
#                  break()
#                repLoci <- newLoci
#                repLociDP <- lociPD[repLoci, ]
#                whichOverlaps <- which(getOverlaps(coordinates = overNulls@annotation, 
#                                                   segments = repLociDP@annotation, overlapType = "within", 
#                                                   whichOverlaps = FALSE, cl = cl))
#        }
#              else break()
#              message(".", appendLF = FALSE)
#            }
#          message("done.", appendLF = TRUE)
#          repNull[overLoci] <- nullPosts$posteriors[, 1] > nullCutoff
        }
        repNull
      }
      message("Establishing likelihoods of nulls...", appendLF = TRUE)
      potnullD@posteriors <- sapply(unique(potnullD@replicates), getNullPosteriors, lociPD = lociPD, potNulls = potnullD)
      potnullD@posteriors <- log(potnullD@posteriors)
      } else potnullD <- NULL

  lociPD@posteriors <- log(lociPD@posteriors)
  
    locMap <- .processPosteriors(lociPD, potnullD, chrs = sD@chrs, aD = aD, lociCutoff = 1, nullCutoff = 1, getLikes = getLikes, cl = cl)
  locMap
      }
