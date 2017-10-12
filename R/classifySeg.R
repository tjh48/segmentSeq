# modification on git from copied files
.constructPriorNulls <- function(sD, aD, samplesize, cl)
    {
        if(inherits(sD, "lociData")) sD <- list(sD)
        nullPL <- lapply(sD, function(ssD) {
            if(is.character(ssD)) {
                message(ssD)
                load(ssD)
            } else subSD <- ssD
            
            subWithinLoc <- .getOverlaps(subSD@coordinates, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,], whichOverlaps = FALSE, overlapType = "within", cl = NULL)
                                        #subWithinLoc <- !is.na(findOverlaps(subSD@coordinates, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0], type = "within", select = "first"))
            
            if(any(subWithinLoc))
                {              
                    subSDWithin <- subSD[subWithinLoc,]
                    if(nrow(subSDWithin@data) == 0) subSDWithin@data <- getCounts(subSDWithin@coordinates, aD, cl = cl)
                    
                    emptyNulls <- gaps(subSD@coordinates[.fastUniques(cbind(as.character(seqnames(subSD@coordinates)), start(subSD@coordinates), as.character(strand(subSD@coordinates))))])
                
                    potnullPriors <- .constructNulls(emptyNulls, subSDWithin, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,], aD = aD, cl = cl, forPriors = TRUE, samplesize = samplesize)
                    seglens(potnullPriors) <- width(potnullPriors@coordinates)
                }
            potnullPriors
        })
        nullPL
    }
               
               
                                        #        chunkNull <- function(subSD) {
                                        #            chunks <- subSD@coordinates$chunk
                                        #            z <- do.call("c", sapply(unique(chunks), function(chunk) {
#                chunkCod <- subSD@coordinates[chunks == chunk]
#                selCC <- which(rowSums(subSD@locLikelihoods[which(chunks == chunk),,drop = FALSE] == 0, na.rm = TRUE) > 0)
#                withinLoc <- which(.getOverlaps(chunkCod, chunkCod[selCC,], whichOverlaps = FALSE, overlapType = "within", cl = #NULL))
#                if(length(withinLoc > 1))
#                    {
#                        emptyNulls <- gaps(chunkCod[.fastUniques(cbind(as.character(seqnames(chunkCod)), start(chunkCod), as.cha#racter(strand(chunkCod))))])
#                        emptyNulls <- emptyNulls[end(emptyNulls) > min(start(chunkCod))]
#                        if(length(emptyNulls) == 0) return(GRanges())
#                        sDWithin <- new("lociData", coordinates = chunkCod[withinLoc])
#                        curNullsWithin <- .constructNulls(emptyNulls, sDWithin, locDef = chunkCod[selCC,], forPriors = NA, sampl#esize = 1, aD = new("alignmentData"), cl = cl)
#                    } else return(GRanges())
#                return(curNullsWithin)
#            }))
#        }
#        
#        if(is.character(sD)) {
#            priorNulls <- do.call("c", lapply(sD, function(sDfile) {
#                load(sDfile)
                                        #                message(sDfile)
#                chunkNull(subSD)
#            }))
#        } else if(inherits(sD, "lociData")) priorNulls <- chunkNull(sD)#
#
#        return(priorNulls)
#    }
    


classifySeg <- function(sD, cD, aD, lociCutoff = 0.9, nullCutoff = 0.9, subRegion = NULL, getLikes = TRUE, lR = FALSE, samplesize = 1e5, largeness = 1e8, tempDir = NULL, recoverFromTemp = FALSE, cl)
    {
        if(missing(aD)) stop("You must supply an aD object.")
        
        if(!is.null(subRegion))
            {
                sD <- sD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(sD@coordinates)) == subRegion$chr[ii] & end(sD@coordinates) >= subRegion$start[ii] & start(sD@coordinates) <= subRegion$end[ii]))),]
                aD <- aD[unlist(lapply(1:nrow(subRegion), function(ii) which(as.character(seqnames(aD@alignments)) == subRegion$chr[ii] & end(aD@alignments) >= subRegion$start[ii] & start(aD@alignments) <= subRegion$end[ii]))),]
            }
        
        if(missing(cD))
            cD <- heuristicSeg(sD, aD, largeness = largeness, getLikes = TRUE, verbose = TRUE, cl = cl)
    
        if(!is.null(tempDir))
            saveSD <- paste0(tempDir, "/sD_locLikes.RData")
        
        if(recoverFromTemp && !is.null(tempDir) && file.exists(saveSD)) {
            load(saveSD)        
        } else {
            sD <- .massiveClassifyLoci(sD = sD, cD = cD, aD = aD, subRegion = subRegion, samplesize = samplesize, lR = lR, lociCutoff = lociCutoff, tempDir = tempDir, largeness = largeness, recoverFromTemp = recoverFromTemp, cl = cl)            
        }

        if(!is.null(tempDir)) nullPriorsFile <- paste(tempDir, "/nullPriors.RData", sep = "") 
        if(recoverFromTemp && file.exists(nullPriorsFile)) {
            load(nullPriorsFile)
        } else {
                                        # construct null priors based on locus classification analysis
            sampleNulls <- .constructPriorNulls(sD, aD, samplesize = samplesize, cl = cl)        
            if(!is.null(tempDir)) save(sampleNulls, file = paste(tempDir, "/sampleNulls.RData", sep = ""))
            
            curNullsWithin <- do.call("c", sampleNulls)
            
            if(length(sampleNulls) == 0) {
                warning("No null regions found within loci; this may indicate a problem with the analysis.")
            }
            
            curNullsWithin@groups = list(replicates(aD))        
            curNullsWithin <- getPriors.NB(curNullsWithin, samplesize = samplesize, verbose = FALSE, cl = cl)
        
            curNullsWithin@priors$weights <- matrix(0, nrow = nrow(curNullsWithin@priors$sampled), ncol = length(unique(curNullsWithin@replicates)))
            
            nullSegPriors <- cD
            nullSegPriors@locLikelihoods <- matrix(nrow = 0, ncol = 0)
            groups(nullSegPriors) <- list(replicates(nullSegPriors))
            if(class(nullSegPriors) == "lociData") {
                nullSegPriors <- getPriors.NB(nullSegPriors, samplesize = samplesize, verbose = FALSE, cl = cl)
            } #else if(class(nullSegPriors) == "methData")
                                        #    nullSegPriors <- getPriors.BB(nullSegPriors, samplesize = samplesize, verbose = FALSE, cl = cl)
            
                                        #nullSampled = list(list(nullSegPriors@priors$sampled), list(curNullsWithin@priors$sampled))
            
            if(!is.null(tempDir)) save(curNullsWithin, nullSegPriors, file = paste(tempDir, "/nullPriorInput.RData", sep = ""))        
            nullPriors <- lapply(levels(aD@replicates), function(rep) {
                repCol <- which(levels(aD@replicates) == rep)
                list(priors = list(list(nullSegPriors@priors$priors[[1]][[repCol]]), list(curNullsWithin@priors$priors[[1]][[repCol]])),
                     weights = list(list((1 - exp(cD@locLikelihoods[nullSegPriors@priors$sampled[,1], repCol]))), list(exp(curNullsWithin@priors$weights[,repCol]))),
                     sampled = list(list(nullSegPriors@priors$sampled), list(curNullsWithin@priors$sampled))
                     )                
            })
            if(!is.null(tempDir)) save(nullPriors, file = nullPriorsFile)
        }
                                        #sD@locLikelihoods <- do.call("DataFrame", lapply(as.list(sD@locLikelihoods), function(x) log(x >= log(lociCutoff))))


        if(inherits(sD, "lociData")) sD <- list(sD)
        lDList <- lapply(seq_along(sD), function(ii) {
            lDfile <- paste0(tempDir, "/lD_", ii, ".RData", sep = "")            
            if(!file.exists(lDfile) | !recoverFromTemp) {
                message("Establishing likelihoods of nulls; Part ", ii, " of ", length(sD))
                subSD <- sD[[ii]]
                if(is.character(subSD)) load(subSD)
                
                nullFile = paste0(tempDir, "/nullD_", ii, ".RData", sep = "")
                if(!file.exists(nullFile) | !recoverFromTemp) {
                    nullD <- .massiveClassifyNulls(subSD = subSD, aD = aD, nullPriors = nullPriors, subRegion = subRegion, samplesize = samplesize, lR = lR, nullCutoff = nullCutoff, tempDir = tempDir, largeness = largeness, recoverOld = FALSE, cl = cl)
                    if(!is.null(tempDir)) save(nullD, file = nullFile)
                } else load(nullFile)
                           
                
            if(nrow(nullD) > 0) {
                nullPD <- nullD[!values(nullD@coordinates)$empty,]
                emptyPD = nullD[values(nullD@coordinates)$empty,]
            } else nullPD <- emptyPD <- nullD
                
                gc()
                lD <- .processPosteriors(lociPD = subSD, nullPD = nullPD, emptyPD = emptyPD, getLikes = FALSE, extendLoci = TRUE, cl = cl)
                if(!is.null(tempDir)) {
                    save(lD, file = lDfile)
                    return(lDfile)
                } else return(lD)
            } else return(lDfile)
        })

        for(ll in 1:length(lDList)) {
            if(is.character(lDList[[ll]])) load(lDList[[ll]]) else lD <- lDList[[ll]]
            if(ll == 1) mlD <- lD
            if(ll > 1) mlD <- c(mlD, lD)
        }

        lD <- mlD
        if(all(is.na(lD@data))) lD@data <- getCounts(lD@coordinates, aD, cl = cl)
        
        if(getLikes & !missing(aD)) lD <- lociLikelihoods(lD, aD, cl = cl) else if(getLikes & missing(aD)) warning("I can't calculate locus likelihoods without an aD object. You can run lociLikelihoods on the output of this function for the same result.")      
        lD
  }
    

.massiveClassifyLoci <- function(sD, cD, aD, subRegion = NULL, samplesize = 1e5, lR = FALSE, lociCutoff = 0.9, largeness = 1e8, tempDir = NULL, recoverFromTemp, cl = cl)
  {
      if(!is.null(tempDir)) dir.create(tempDir, showWarnings = FALSE) else recoverFromTemp = FALSE
          
    if(lR) {
      pET <- "none"
      lociCutoff <- nullCutoff <- 0.5
    } else {
      locps <- sapply(levels(cD@replicates), function(rep) mean(exp(cD@locLikelihoods[rowSums(cD@data[,cD@replicates == rep, drop = FALSE]) > 0, levels(cD@replicates) == rep]), na.rm = TRUE))
      pET <- "BIC"
    }
      prepD_file <- paste(tempDir, "/prepD.RData", sep = "")
      if(!recoverFromTemp | !file.exists(prepD_file)) {
          message("Finding candidate priors...", appendLF = FALSE)
          
          subLoc <- !is.na(findOverlaps(sD@coordinates, cD@coordinates, type = "within", select = "first"))
                                        #    subLoc <- getOverlaps(coordinates = sD@coordinates,
                                        #                          segments = cD@coordinates, overlapType = "within", 
                                        #                          whichOverlaps = FALSE, cl = NULL)
          subLoc <- which(subLoc)[.filterSegments(sD@coordinates[subLoc,], runif(sum(subLoc)))]
          prepD <- sD[subLoc,]
                                        #prepD <- .convertSegToLoci(sD[subLoc,])
          
          if(all(is.na(prepD@data))) prepD@data <- getCounts(prepD@coordinates, aD, cl = cl)
          message("done.")
          
          groups(prepD) <- list(prepD@replicates)
          if(class(aD) == "alignmentData") {
              prepD <- getPriors.NB(prepD, samplesize = samplesize, verbose = TRUE, cl = cl)
          } else if(class(aD) == "alignmentMeth") {
              densityFunction(prepD) <- bbDensity
              libsizes(prepD) <- matrix(1, nrow = ncol(prepD), ncol = 2)
              prepD@data <- round(prepD@data)
              prepD <- getPriors(prepD, samplesize = samplesize, verbose = TRUE, cl = cl)
          }
                
          weights <- cD@locLikelihoods[subjectHits(findOverlaps(prepD@coordinates[prepD@priors$sampled[,1], ], cD@coordinates, type = "within")),,drop = FALSE]
          repWeights <- sapply(levels(sD@replicates), function(rep) {
              repWeights <- weights[,levels(sD@replicates) == rep]
              repWeights[rowSums(prepD[prepD@priors$sampled[,1], which(sD@replicates == rep)]@data) == 0] <- NA
              repWeights
          })     
          
          if(!is.null(tempDir)) save(repWeights, subLoc, prepD, file = prepD_file)
      } else load(prepD_file)
      
      sDsplit <- .splitSD(sD, sD@coordinates$chunk, largeness)      
      message("Segmentation split into ", length(sDsplit), " parts.")
      
      chunkLL <- lapply(1:length(sDsplit), function(ii) {
          sDfile <- paste(tempDir, "/subSD_", ii, ".RData", sep = "")
          if(!is.null(tempDir) && file.exists(sDfile) && recoverFromTemp) return(sDfile)
          
          message("Establishing likelihoods of loci; Part ", ii, " of ", length(sDsplit))
          x <- sD[sDsplit[[ii]],]
          subLoc <- x@locLikelihoods
          if (is.null(subRegion)) {
              locSubset <- 1:nrow(x)
          } else {
              locSubset <- sort(unique(c(unlist(apply(subRegion, 1, 
                                                      function(sR) which(as.character(seqnames(x@coordinates)) == as.character(sR[1]) &
                                                                             start(x@coordinates) >= as.numeric(sR[2]) &
                                                                                 end(x@coordinates) <= as.numeric(sR[3])))))))
          }
          
          if(nrow(x@data) == 0) {
              if(inherits(aD, "alignmentData")) x@data <- getCounts(x@coordinates, aD, cl = cl)
              if(inherits(aD, "alignmentMeth")) x@data <- array(do.call("c", getCounts(x@coordinates, aD, cl = cl)), dim = c(dim(x), 2))
          }
          
          subLoc <- .classifyLoci(potlociD = x, prepD = prepD, repWeights = repWeights, subLoc = subLoc, locSubset = locSubset, lR = lR, locps = locps, cl = cl)
          
          #if(!is.null(tempDir)) save(subLoc, file = paste(tempDir, "/subLocLikes_", ii, ".RData", sep = ""))
          subLoc <- log(subLoc >= log(lociCutoff))

          if(!is.null(tempDir)) {
              chunkSum <- sapply(split(subLoc == 0, x@coordinates$chunk), sum)
              chunkLoci <- unique(x@coordinates$chunk)[which(chunkSum > 0)]
              subSD <- x[which(x@coordinates$chunk %in% chunkLoci),,drop = FALSE]
              subSD@locLikelihoods <- subLoc[which(x@coordinates$chunk %in% chunkLoci),,drop = FALSE]
              save(subSD, file = sDfile)
              return(sDfile)
          } else return(subLoc)
      })

      if(is.null(tempDir)) {
          sD@locLikelihoods <- do.call("rbind", chunkLL)
          return(sD)
      } else return(do.call("c", chunkLL))
  }


.massiveClassifyNulls <- function(subSD, aD, nullPriors, subRegion = NULL, samplesize = 1e5, lR = FALSE, nullCutoff = 0.9, largeness = 1e8, tempDir, recoverOld = FALSE, cl = cl)
  {   
      
      subWithinLoc <- .getOverlaps(subSD@coordinates, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,], whichOverlaps = FALSE, overlapType = "within", cl = NULL)
                                        #subWithinLoc <- !is.na(findOverlaps(subSD@coordinates, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0], type = "within", select = "first"))
      
      if(any(subWithinLoc))
          {              
              subSDWithin <- subSD[subWithinLoc,]
              if(nrow(subSDWithin@data) == 0) subSDWithin@data <- getCounts(subSDWithin@coordinates, aD, cl = cl)
              
              emptyNulls <- gaps(subSD@coordinates[.fastUniques(cbind(as.character(seqnames(subSD@coordinates)), start(subSD@coordinates), as.character(strand(subSD@coordinates))))])

              #potnullPriors <- .constructNulls(emptyNulls, subSDWithin, subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,], aD = aD, cl = cl, forPriors = TRUE, samplesize = samplesize)
              
              potnullD <- .constructNulls(emptyNulls, sDWithin = subSDWithin, locDef = subSD@coordinates[rowSums(subSD@locLikelihoods == 0) > 0,], aD = aD, cl = cl)

              postOver <- which(.getOverlaps(coordinates = potnullD@coordinates, segments = subSD@coordinates, overlapType = "within", whichOverlaps = FALSE, cl = NULL))
              
              potnullD@locLikelihoods <- (do.call("cbind", ((lapply(levels(potnullD@replicates), .classifyNulls, lociPD = subSD, potNulls = potnullD, lR = lR, nullPriors = nullPriors, postOver = postOver, nullCutoff = nullCutoff, cl = cl)))))
          } else potnullD <- NULL
      if(class(aD) == "alignmentData") {
          emptyD <- rowSums(potnullD@data) == 0
      } else if(class(aD) == "alignmentMeth") {
          emptyD <- rowSums(potnullD@Cs) == 0
      }
        
      values(potnullD@coordinates)$empty <- FALSE
      values(potnullD@coordinates)$empty[emptyD] <- TRUE

#      potnullD <- potnullD[which(.rowSumDF(do.call("cbind", lapply(1:ncol(potnullD@locLikelihoods), function(ii) potnullD@locLikelihoods[,ii] > -Inf)), na.rm = TRUE) > 0 | emptyD),]
      potnullD <- potnullD[which(rowSums(potnullD@locLikelihoods == 0, na.rm = TRUE) > 0 | emptyD),]
      
#      if(class(potnullD) == "segData") {
#          potnullD@data <- (matrix(ncol = length(potnullD@replicates), nrow = 0))
#      } else if(class(potnullD) == "segMeth") {
#          potnullD@Cs <- potnullD@Ts <- (matrix(ncol = length(potnullD@replicates), nrow = 0))
                                        #      }
      
#      if(!is.null(tempDir)) save(potnullD, file = paste(tempDir, "/subNull_", ii, ".RData", sep = ""))
      
#    } else potnullD <- new(class(sD))
      
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
    

    libsizes = libsizes(pcD)
    NBpriors <- pcD@priors$priors
    seglens <- seglens(pcD)
    numintSamp <- cbind(pcD@priors$sampled, pcD@priors$weights)
    numintSamp <- numintSamp[order(numintSamp[,2]),]

    priorWeights <- constructWeights()
    NBpriors <- NBpriors[rowSums(priorWeights) != 0,]
    priorWeights <- priorWeights[rowSums(priorWeights) != 0,]
    
    if(!is.null(cl))
      {
        getLikelihoodsEnv <- new.env(parent = .GlobalEnv)
        environment(NBdens) <- getLikelihoodsEnv
        environment(BBdens) <- getLikelihoodsEnv
        
        clusterExport(cl, c("NBpriors", "numintSamp", "libsizes", "priorWeights"), envir = environment())
      }

    if(length(dim(pcD@data)) == 2)
      {
        if (is.null(cl)) {
          ps <- t(apply(cbind(1:nrow(pcD@data), seglens, pcD@data)[subset,,drop = FALSE], 1, NBdens))
        } else {      
          ps <- parRapply(cl, cbind(1:nrow(pcD@data), seglens, pcD@data)[subset,, drop = FALSE], NBdens)
          ps <- matrix(ps, ncol = 2, byrow = TRUE)
        }
      } else if(length(dim(pcD@data)) == 3) {
        if (is.null(cl)) {
          ps <- t(apply(cbind(1:nrow(pcD@data), pcD@data[,,1], pcD@data[,,2])[subset,,drop = FALSE], 1, BBdens))
        } else {      
          ps <- parRapply(cl, cbind(1:nrow(pcD@data), pcD@data[,,1], pcD@data[,,2])[subset,, drop = FALSE], BBdens)
          ps <- matrix(ps, ncol = 2, byrow = TRUE)
        }
      }

    rps <- matrix(NA, ncol = 2, nrow = nrow(pcD@data))
    rps[subset,] <- ps
    
    return(rps)
  }


.classifyNulls <- function(rep, lociPD, potNulls, lR, nullPriors, postOver,
                           nullCutoff, cl) {
  
    message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
    repCol <- which(levels(potNulls@replicates) == rep)
    repD <- potNulls[, potNulls@replicates == rep]
    
  #repD <- .convertSegToLoci(repD)

  if(class(repD) == "lociData") {
      repD@priorType = "NB"
      densityFunction(repD) <- nbinomDensity
      seglens(repD) <- width(repD@coordinates)
  } else if(class(repD) == "methData") repD@priorType = "BB"    
    
  replicates(repD) <- as.factor(rep(1, ncol(repD)))
  groups(repD) <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
  
  repD@priors$priors <- nullPriors[[repCol]]$priors
  
  repD@priors$sampled <- nullPriors[[repCol]]$sampled
  repD@priors$weights <- nullPriors[[repCol]]$weights
  repD@priors$sampled[[1]][[1]][,1] <- -1
  repD@priors$sampled[[2]][[1]][,1] <- -1
  
  repCol <- which(levels(lociPD@replicates) == rep)
  repNull <- rep(NA, nrow(repD))
  repLoci <- which(lociPD@locLikelihoods[, repCol] == 0)
  repLociDP <- lociPD[repLoci, ]    
  
#postOver <- which(!is.na(findOverlaps(repD@coordinates, repLociDP@coordinates, type = "within", select = "first")))
    
  overLoci <- unique(c(postOver,
                       unlist(lapply(levels(seqnames(repD@coordinates)), function(chrom)
                                     which(seqnames(repD@coordinates) == chrom &
                                           (end(repD@coordinates) + 1) %in% c(start(repLociDP@coordinates[seqnames(repLociDP@coordinates) == chrom]), seqlengths(repLociDP@coordinates)[levels(seqnames(repD@coordinates)) == chrom]) &
                                           (start(repD@coordinates) - 1) %in% c(end(repLociDP@coordinates[seqnames(repLociDP@coordinates) == chrom]), 0))
                                     ))))
  
    overLoci <- sort(overLoci[!is.na(overLoci)])
  
    if (length(postOver) > 0) {
        overNulls <- repD[overLoci,]
        overNulls <- overNulls[order(overNulls@coordinates),]
        nullsWithin <- repD[postOver, ]
        lD <- rep(NA, nrow(overNulls))
    
        if(lR) {
            prs <- c(0.5, 0.5)
            nullCutoff <- 0.5
        } else {
            nullFilter <- .filterSegments(overNulls@coordinates, runif(nrow(overNulls)))              
            if(class(overNulls) == "lociData") {
                fD <- getLikelihoods(cD = overNulls[nullFilter,], bootStraps = 1, 
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
            prs <- colMeans(exp(fD@posteriors))
        }    

        fillInLD <- function(lD, selMin, prs, nullTest)
            {
                if(any(is.na(lD)))
                    {
                        if(missing(nullTest)) nullTest = rep(TRUE, length(lD))
                        if(selMin & sum(is.na(lD) & nullTest) > 1e4)
                            subSel <- .fastUniques(cbind(as.character(seqnames(overNulls@coordinates[is.na(lD) & nullTest])), start(overNulls@coordinates[is.na(lD) & nullTest]))) else subSel <- TRUE
                        
                        nullFilter <- which(is.na(lD) & nullTest)[subSel]
                        
                        if(class(overNulls) == "lociData") {
                            fD <- getLikelihoods(cD = overNulls[nullFilter,], bootStraps = 1, 
                                                 verbose = FALSE, 
                                                 subset = NULL, priorSubset = NULL, pET = "none", prs = c(prs[1], 1 - prs[1]), 
                                                 cl = cl)
                        } else if(class(overNulls) == "methData") {
                                        #              fD <- getLikelihoods.BB(cD = overNulls[nullFilter,], bootStraps = 1, 
                                        #                                      verbose = FALSE, 
                                        #                                      subset = NULL, priorSubset = NULL, pET = "none", prs = c(prs[1], 1 - prs[1]), 
                                        #                                      cl = cl)
                            stop("classifySeg not currently available for methylation data.")
                        }
                            
                        lD[nullFilter] <- fD@posteriors[,1] > log(nullCutoff)
                    }
                lD
            }
        
#        if(any(lD, na.rm = TRUE)) {
            repeat {
                locOverNull <- .getOverlaps(lociPD@coordinates, overNulls@coordinates[which(lD),], overlapType = "contains", whichOverlaps = FALSE, cl = NULL)
                nullWithinLoc <- .getOverlaps(overNulls@coordinates, lociPD@coordinates[which(!locOverNull)], overlapType = "within", whichOverlaps = FALSE, cl = NULL)
                if(any(nullWithinLoc) && any(is.na(lD[nullWithinLoc]))) {
                    lD <- fillInLD(lD, selMin = TRUE, prs = prs, nullTest = nullWithinLoc)
                } else break
            }
#        }

    message("...done.", appendLF = TRUE)
    
    repNull[overLoci] <- lD
  }
                                        #Rle(log(repNull))
  log(repNull)
}


.classifyLoci <- function(potlociD, prepD, repWeights, subLoc, locSubset, lR, locps, cl)
  {
    getLikeLoci <- function(rep, potlociD, priors, repWeights, sampled, ...) {
      message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)
      
#      if(class(potlociD) == "segData")
#        {
#          repD <- new("lociData",
#                      data = sapply(which(potlociD@replicates == rep), function(ii) as.integer(potlociD@data[,ii])),
#                      seglens = width(potlociD@coordinates),
#                      libsizes = potlociD@libsizes[potlociD@replicates == rep],
#                      replicates = as.factor(rep(1, sum(potlociD@replicates == rep))),
#                      coordinates = potlociD@coordinates,
#                      groups = list(rep(1, sum(potlociD@replicates == rep))))
#          repD@priorType = "NB"
#        } else if(class(potlociD) == "segMeth") {
#          repD <- new("methData",
#                      data = round(potlociD@Cs[,which(potlociD@replicates == rep), drop = FALSE]),
#                      pairData = round(potlociD@Ts[,which(potlociD@replicates == rep), drop = FALSE]),
#                      replicates = as.factor(rep(1, sum(potlociD@replicates == rep))),
#                      coordinates = potlociD@coordinates,
#                      groups = list(rep(1, sum(potlociD@replicates == rep))),
#                      libsizes = potlociD@nonconversion[potlociD@replicates == rep] + 1,
#                      pairLibsizes = 1 - potlociD@nonconversion[potlociD@replicates == rep])
#          repD@priorType = "BB"
                                        #        }

      repD <- potlociD[,which(potlociD@replicates == rep)]
      groups(repD) <- list(rep(1, ncol(repD)))
      
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
          if(is.null(weightFactors)) weightFactors <- rep(1, nrow(sampled))
          repD@priors$weights <- sampled[,"weights"] * repD@priors$weights * weightFactors
          repD@priors$sampled <- sampled[,1:2]
          
          subset <- intersect(locSubset, which(rowSums(repD@data) > 0))        
          lD <- matrix(NA, ncol = 2, nrow = nrow(repD))
          
          if(length(subset) > 0) {
              if(length(dim(potlociD@data)) == 2)
                  {              
                      orddat <- do.call("order", c(lapply(1:ncol(as.matrix(seglens(repD))), function(ii) seglens(repD)), lapply(1:ncol(repD), function(ii) repD@data[,ii])))                      
                      whunq <- .fastUniques(cbind(seglens(repD), repD@data)[orddat,])
                      lD <- .getLocLikelihoods(repD, subset = intersect(orddat[whunq], subset), cl = cl)
                      
                      lD[orddat,] <- lD[orddat[rep(which(whunq), diff(c(which(whunq), length(whunq) + 1)))],]
              
                  } else if(length(dim(potlociD@data)) == 3) {
                      orddat <- do.call("order", c(lapply(1:ncol(repD@data), function(ii) repD@data[,ii,1]), lapply(1:ncol(repD@data), function(ii) repD@data[,ii,2])))
                      whunq <- .fastUniques(cbind(repD@data[,,1], repD@data[,,2])[orddat,])
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
      repLoci #Rle(repLoci)
    }

    weightFactors <- prepD@priors$weights
    priors <- prepD@priors$priors    
    sampled <- prepD@priors$sampled
    sampled[,1] <- subLoc[sampled[,1]]
    
    message("Establishing likelihoods of loci...", appendLF = TRUE)  
    locLikelihoods <- do.call("cbind", lapply(levels(potlociD@replicates), getLikeLoci, potlociD = potlociD, priors = priors, repWeights = repWeights, sampled = sampled))
    colnames(locLikelihoods) <- levels(potlociD@replicates)
    
    locLikelihoods
  }
