lociLikelihoods <- function(cD, aD, bootStraps = 5, inferNulls = TRUE, cl)
  {    
    loci <- cD@annotation
    lociLens <- loci$end - loci$start + 1
    countLoci <- getCounts(aD = aD, segments = loci, cl = cl)

    if(inferNulls)
      {
        nulls <- matrix(unlist(sapply(aD@chrs, function(x) with(loci, rbind(x, c(1, end[chr == x] + 1), c(start[chr == x] - 1, aD@chrlens[aD@chrs == x]), c(NA, start[chr == x]), c(end[chr == x], NA))))), ncol = 5, byrow = TRUE)
        nulls <- data.frame(chr = I(nulls[,1]), start = as.integer(nulls[,2]), end = as.integer(nulls[,3]), leftLocusExtend = nulls[,4], rightLocusExtend = nulls[,5])
        nullLens <- nulls$end - nulls$start + 1
        nulls[nullLens <= 0 | is.na(nullLens),] <- NA
        nullLens[nullLens <= 0 | is.na(nullLens)] <- NA
        nulls <- nulls[rowSums(is.na(nulls)) == 0,1:3]
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

        mD <- mD[!is.na(mD@seglens),]    
        mD@groups <- list(mD@replicates)
        mD <- getPriors.NB(mD, verbose = FALSE, cl = cl)
        
        locW <- as.numeric(mD@annotation$segType[mD@priors$sampled[,1]] == "null")
        
        lociWeights <- sapply(unique(mD@replicates), function(rep, locW) {
          
          message(paste("Getting likelihoods for replicate group", rep), appendLF = FALSE)
          repD <- mD[,mD@replicates == rep]
          repD@replicates <- rep(1, ncol(repD))
          repD@groups <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[rep]]), list(mD@priors$priors[[1]][[rep]]))
          repD@priors$weights <- list(list(locW), list(1 - locW))
          
          repD <- getLikelihoods.NB(cD = repD, bootStraps = bootStraps, verbose = FALSE, cl = cl)
          
          message("done!", appendLF = TRUE)
          repD@posteriors[,2]
        }, locW = locW)

      } else {
        mD <- cD
        class(mD) == "postSeg"

        mD@groups <- list(mD@replicates)
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
    mD@annotation <- subset(mD@annotation, select = c(chr, start, end))

    mD <- mD[with(mD@annotation, order(chr, start, end)),]
    
    mD
  }
