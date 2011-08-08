lociLikelihoods <- function(cD, aD, newCounts = FALSE, bootStraps = 1, inferNulls = TRUE, nasZero = FALSE, usePosteriors = TRUE, cl)
  {    
    loci <- cD@annotation
    lociLens <- loci$end - loci$start + 1
    if(newCounts) countLoci <- getCounts(aD = aD, segments = loci, cl = cl) else countLoci <- cD@data

    chrs <- aD@chrs
    
    if(inferNulls)
      {
        nulls <- matrix(unlist(sapply(1:nrow(chrs), function(ii) with(loci, rbind(ii, c(1, end[chr == chrs$chr[ii]] + 1),
                                                                                  c(start[chr == chrs$chr[ii]] - 1, chrs$len[chrs$chr == chrs$chr[ii]]),
                                                                                  c(NA, start[chr == chrs$chr[ii]]),
                                                                                  c(end[chr == chrs$chr[ii]], NA))))), ncol = 5, byrow = TRUE)
        nulls <- data.frame(chr = chrs$chr[nulls[,1]], start = as.integer(nulls[,2]), end = as.integer(nulls[,3]), leftLocusExtend = nulls[,4], rightLocusExtend = nulls[,5])
        nullLens <- nulls$end - nulls$start + 1
        nulls[nullLens <= 0 | is.na(nullLens),] <- NA
        nullLens[nullLens <= 0 | is.na(nullLens)] <- NA
        nulls <- nulls[rowSums(is.na(nulls)) == 0,1:3,drop = FALSE]
        nullLens <- nulls$end - nulls$start + 1    
        countNulls <- getCounts(aD = aD, segments = nulls, cl = cl)
        
        mD <- new("postSeg", data = rbind(countNulls, countLoci),
                  seglens = c(nullLens, lociLens),
                  libsizes = aD@libsizes,
                  replicates = aD@replicates, 
                  annotation = data.frame(chr = I(c(as.character(nulls$chr), as.character(cD@annotation$chr))),
                    start = c(nulls$start, cD@annotation$start),
                    end = c(nulls$end, cD@annotation$end),
                    segType = c(rep("null", nrow(nulls)), rep("locus", nrow(cD)))))
        posteriors = rbind(matrix(-Inf, nrow = nrow(countNulls), ncol = length(unique(aD@replicates))), cD@posteriors)
      } else {
        posteriors <- cD@posteriors
        mD <- new("postSeg", data = countLoci,
                  seglens = lociLens,
                  libsizes = aD@libsizes,
                  replicates = aD@replicates, 
                  annotation = loci)
      }

    mD <- mD[rowSums(is.na(mD@seglens)) == 0,]    
    mD@groups <- list(mD@replicates)
    mD <- getPriors.NB(mD, verbose = FALSE, cl = cl)

    if(usePosteriors)
      {
        lociWeights <- sapply(unique(mD@replicates), function(rep) {
#          locW <- as.numeric(mD@annotation$segType[mD@priors$sampled[,1]] == "null")
#          locW[mD@annotation$segType[mD@priors$sampled[,1]] == "locus"] <- 1 - exp(posteriors[mD@priors$sampled[mD@annotation$segType[mD@priors$sampled[,1]] == "locus",1],rep])
          locW <- 1 - exp(posteriors[mD@priors$sampled[,1],rep])
          locW <- list(list(locW), list(1 - locW))
          
          if(nasZero) {
            locW <- lapply(locW, function(x) {
              y <- x[[1]]
              y[is.na(y)] <- 0
              list(y * mD@priors$weights)
            })} else locW <- lapply(1:2, function(ii) {
              y <- locW[[ii]][[1]]
              y[is.na(y)] <- as.numeric(ii == 1)
              list(y * mD@priors$weights)
            })
            
          message(paste("Getting likelihoods for replicate group", rep), appendLF = FALSE)
          repD <- mD[,mD@replicates == rep]
          repD@replicates <- rep(1, ncol(repD))
          repD@groups <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[rep]]), list(mD@priors$priors[[1]][[rep]]))
          repD@priors$weights <- locW
          
          repD <- getLikelihoods.NB(cD = repD, bootStraps = bootStraps, verbose = FALSE, cl = cl)
          
          message("done!", appendLF = TRUE)
          repD@posteriors[,2]
        })
      } else {
        mD <- getPriors.NB(mD, verbose = FALSE, cl = cl)
        
        lociWeights <- sapply(unique(mD@replicates), function(rep) {
          
          message(paste("Getting likelihoods for replicate group", rep), appendLF = FALSE)
          repD <- mD[,mD@replicates == rep]
          repD@replicates <- rep(1, ncol(repD))
          repD@groups <- list(rep(1, ncol(repD)))
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[rep]]))
          repD@priors$weights <- NULL
          
          repD <- getLikelihoods.NB(cD = repD, bootStraps = bootStraps, nullData = TRUE, verbose = FALSE, cl = cl)
          
          message("done!", appendLF = TRUE)
          repD@posteriors[,1]
        })
      }
          

    mD@posteriors <- lociWeights
    mD@annotation <- subset(mD@annotation, select = c("chr", "start", "end"))

    mD <- mD[with(mD@annotation, order(chr, start, end)),]
    
    mD
  }
