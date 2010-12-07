classifySeg <- function (sD, cD, aD, lociCutoff = 0.5, nullCutoff = 0.9, subRegion = NULL, getLikes = TRUE, cl) 
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
        cD <- clustSeg(sD = sD, aD = aD, cl = cl)
#    if (class(cD) == "countData") 
#        cD <- lociLikelihoods(cD = cD, aD = aD, cl = cl)
    
    dupStarts <- which(fastUniques(cbind(sD@segInfo$chr, sD@segInfo$start)))                                         
    emptyNulls <- with(sD@segInfo, data.frame(chr = chr[dupStarts], start = start[dupStarts] - leftSpace[dupStarts], end = start[dupStarts] - 1))
    emptyNulls <- emptyNulls[emptyNulls$start > 1, ]

    potnullD <- with(sD@segInfo, new("postSeg",
                                              data = rbind(matrix(0, ncol = ncol(sD), nrow = nrow(emptyNulls)),
                                                sD@data[leftSpace > 0 & rightSpace > 0,],
                                                sD@data[leftSpace > 0,],
                                                sD@data[rightSpace > 0,]),
                                     seglens = c(emptyNulls$end - emptyNulls$start + 1,
                                       (end - start + 1 + leftSpace + rightSpace)[leftSpace > 0 & rightSpace > 0],
                                       (end - start + 1 + leftSpace)[leftSpace > 0],
                                       (end - start + 1 + rightSpace)[rightSpace > 0]),
                                     libsizes = sD@libsizes,
                                     replicates = sD@replicates,
                                     annotation = data.frame(chr = c(emptyNulls$chr, chr[leftSpace > 0 & rightSpace > 0], chr[leftSpace > 0], chr[rightSpace > 0]),
                                       start = c(emptyNulls$start, (start - leftSpace)[leftSpace > 0 & rightSpace > 0], (start - leftSpace)[leftSpace > 0], start[rightSpace > 0]),
                                       end = c(emptyNulls$end, (end + rightSpace)[leftSpace > 0 & rightSpace > 0], end[leftSpace > 0], (end + rightSpace)[rightSpace > 0]),
                                       segType = "potnul",
                                       nullClass = c(rep("empty", nrow(emptyNulls)), rep("expBoth", sum(leftSpace > 0 & rightSpace > 0)), rep("expLeft", sum(leftSpace > 0)), rep("expRight", sum(rightSpace > 0))))
                                     ))

    potlociD <- with(sD@segInfo, new("postSeg",
                                     data = sD@data,
                                     seglens = end - start + 1,
                                     libsizes = sD@libsizes,
                                     replicates = sD@replicates,
                                     annotation = data.frame(chr = chr, start = start, end = end, segType = "potloc")))

    getPosts <- function(rep, psD, nulPost, prs, ...) {
        message(paste("\t\t...for replicate group ", rep, "...", sep = ""), appendLF = FALSE)

        repD <- psD[, psD@replicates == rep]

        repD@priorType = "NB"
        repD@replicates <- rep(1, ncol(repD))

        if (is.null(subRegion)) {
          subset <- 1:nrow(repD)
        } else {
          subset <- sort(unique(c(unlist(apply(subRegion, 1, 
                                               function(sR) which(as.character(repD@annotation$chr) == as.character(sR[1]) &
                                                                  repD@annotation$start >= as.numeric(sR[2]) &
                                                                  repD@annotation$end <= as.numeric(sR[3])))))))
        }
        
        if (!nulPost) prepD <- repD[rowSums(repD@data) > 0, ] else prepD <- repD[repD@annotation$nullClass %in% c("empty", "expBoth"),]

        prepD <- prepD[getOverlaps(coordinates = prepD@annotation, 
                                   segments = cD@annotation, overlapType = "within", 
                                   whichOverlaps = FALSE, cl = cl), ]
        prepD <- prepD[filterSegments(prepD@annotation, runif(nrow(prepD))), ]
        prepD@groups <- list(rep(1, ncol(prepD)))
        prepD <- getPriors.NB(prepD, verbose = FALSE, cl = cl)
        
        withinWeights <- exp(cD@posteriors[unlist(getOverlaps(coordinates = prepD@annotation[prepD@priors$sampled[,1], ], segments = cD@annotation, overlapType = "within", cl = cl)), rep])
        prepD@priors$sampled[, 1] <- -1
                
        repD@groups <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))        

        if (!nulPost)
          {
            priorSubset <- which(rowSums(repD@data) > 0)[filterSegments(repD@annotation[which(rowSums(repD@data) > 0), ], runif(sum(rowSums(repD@data) > 0)))]
            subset <- subset[rowSums(repD@data[subset, , drop = FALSE]) > 0]
            repD@priors$priors <- list(list(prepD@priors$priors[[1]][[1]]), list(prepD@priors$priors[[1]][[1]]))
            repD@priors$sampled <- prepD@priors$sampled

            if (all(withinWeights == 1)) {
              repD@posteriors <- matrix(NA, nrow = nrow(repD), ncol = 2, byrow = TRUE)
              repD@posteriors[subset,] <- matrix(c(-Inf, 0), nrow = length(subset), ncol = 2, byrow = TRUE)
            } else if (all(withinWeights == 0)) {
              repD@posteriors <- matrix(NA, nrow = nrow(repD), ncol = 2, byrow = TRUE)
              repD@posteriors[subset,] <- matrix(c(-Inf, 0), nrow = length(subset), ncol = 2, byrow = TRUE)
            } else {
              repD@priors$weights <- list(list(1 - withinWeights), list(withinWeights))
              repD <- getLikelihoods.NB(cD = repD, bootStraps = 1, 
                                        verbose = FALSE, prs = c(prs, 1 - prs), returnPD = nulPost, 
                                        subset = subset, priorSubset = priorSubset, pET = "BIC", 
                                        cl = cl)
            }
          } else {
            nullPriors <- cD
            nullPriors@posteriors <- matrix(nrow = 0, ncol = 0)
            nullPriors <- nullPriors[,nullPriors@replicates == rep]
            nullPriors@replicates <- rep(1, ncol(nullPriors))
            nullPriors <- getPriors.NB(nullPriors, verbose = FALSE, cl = cl)
            
            repD@priors$priors <- list(list(nullPriors@priors$priors[[1]][[1]]), list(prepD@priors$priors[[1]][[1]]))
            repD@priors$sampled <- list(list(nullPriors@priors$sampled), list(prepD@priors$sampled))
            repD@priors$weights <- list(list((1 - exp(cD@posteriors[nullPriors@priors$sampled[,1], rep])) * nullPriors@priors$weights), list(withinWeights))
            repD@priors$sampled[[1]][[1]][,1] <- -1
            repD <- getLikelihoods.NB(cD = repD, bootStraps = 1, 
                                      verbose = FALSE, prs = c(prs, 1 - prs), returnPD = nulPost, 
                                      subset = subset, priorSubset = NULL, pET = "BIC", 
            cl = cl)
          }
        message("...done.", appendLF = !nulPost)
        if (!nulPost) return(repD@posteriors[, 2]) else return(repD)
      }
        
    message("Establishing likelihoods of loci...", appendLF = TRUE)
    potlociD@posteriors <- sapply(unique(potlociD@replicates), getPosts, psD = potlociD, nulPost = FALSE)
    lociPD <- potlociD[which(rowSums(potlociD@posteriors > log(lociCutoff), na.rm = TRUE) > 0), ]
    
    overLoci <- unique(c(which(getOverlaps(coordinates = potnullD@annotation, segments = lociPD@annotation, overlapType = "within", whichOverlaps = FALSE, cl = cl)),
                         unlist(lapply(unique(potnullD@annotation$chr), function(chrom)
                                        which(potnullD@annotation$chr == chrom)[with(potnullD@annotation[potnullD@annotation$chr == chrom,],
                                             which(nullClass == "empty" &
                                                   (end + 1) %in% c(lociPD@annotation$start[lociPD@annotation$chr == chrom], aD@chrlens[aD@chrs == chrom]) &
                                                   (start - 1) %in% c(lociPD@annotation$end[lociPD@annotation$chr == chrom], 0)))]
                                        ))))
                       
    nullPD <- potnullD[overLoci,]

    getNullPosteriors <- function(rep, lociPD, potNulls) {
        posteriors <- rep(NA, nrow(potNulls))
        repLoci <- which(lociPD@posteriors[, rep] > log(lociCutoff))
        repLociDP <- lociPD[repLoci, ]
        postOver <- which(getOverlaps(coordinates = potNulls@annotation, segments = repLociDP@annotation, overlapType = "within", whichOverlaps = FALSE, cl = cl))
        overLoci <- unique(c(postOver,
                             unlist(lapply(unique(potNulls@annotation$chr), function(chrom)
                                           which(potNulls@annotation$chr == chrom)[with(potNulls@annotation[potNulls@annotation$chr == chrom,],
                                                   which(nullClass == "empty" &
                                                         (end + 1) %in% c(repLociDP@annotation$start[repLociDP@annotation$chr == chrom], aD@chrlens[aD@chrs == chrom]) &
                                                         (start - 1) %in% c(repLociDP@annotation$end[repLociDP@annotation$chr == chrom], 0)))]
                                           ))
                             ))
        
        if (length(postOver > 0)) {
            overNulls <- potNulls[overLoci,]
            nullsWithin <- potNulls[postOver, ]
            whichOverlaps <- 1:nrow(nullsWithin)
            nullsPD <- getPosts(rep, psD = overNulls, nulPost = TRUE)
            message(" Refining priors...", appendLF = FALSE)
            for (ii in 1:100) {
#              priorSubset <- whichOverlaps[filterSegments(nullsWithin@annotation, runif(nrow(nullsWithin)))]

              prs <- mean(exp(repLociDP[filterSegments(repLociDP@annotation, repLociDP@seglens, decreasing = TRUE),]@posteriors[,rep]))
              prs <- c(1 - prs, prs)
                
              nullPosts <- getPosteriors(ps = nullsPD, pET = "none", groups = list(rep(1, sum(lociPD@replicates == rep)), rep(1, sum(lociPD@replicates == rep))), priorSubset = NULL, prs = prs, cl = cl)
              nullsDP <- overNulls[which(nullPosts$posteriors[,1] > log(nullCutoff)), ]
                if (nrow(nullsDP) > 0) {
                  newLoci <- which(!getOverlaps(coordinates = lociPD@annotation, segments = nullsDP@annotation, overlapType = "contains", whichOverlaps = FALSE, cl = NULL) & lociPD@posteriors[, rep] > log(lociCutoff))
                  if (all(newLoci %in% repLoci) & all(repLoci %in% newLoci)) 
                    break()
                  repLoci <- newLoci
                  repLociDP <- lociPD[repLoci, ]
                  whichOverlaps <- which(getOverlaps(coordinates = overNulls@annotation, 
                                                     segments = repLociDP@annotation, overlapType = "within", 
                                                     whichOverlaps = FALSE, cl = cl))
                  nullsWithin <- overNulls[whichOverlaps, ]
                }
                else break()
                message(".", appendLF = FALSE)
            }
            message("done.", appendLF = TRUE)
            posteriors[overLoci] <- nullPosts$posteriors[, 1]
        }
        posteriors
    }
    message("Establishing likelihoods of nulls...", appendLF = TRUE)
    nullPD@posteriors <- sapply(unique(nullPD@replicates), getNullPosteriors, lociPD = lociPD, potNulls = nullPD)
                                     

    .processPosteriors(lociPD, nullPD, aD = aD, lociCutoff = lociCutoff, nullCutoff = nullCutoff, getLikes = getLikes, cl = cl)
  }
