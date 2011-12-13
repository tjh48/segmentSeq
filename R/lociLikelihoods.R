lociLikelihoods <- function(cD, aD, newCounts = FALSE, bootStraps = 1, inferNulls = TRUE, nasZero = FALSE, usePosteriors = TRUE, cl)
  {    
    loci <- cD@coordinates
    lociLens <- width(loci)
    if(newCounts) countLoci <- getCounts(aD = aD, segments = loci, as.matrix = TRUE, cl = cl) else countLoci <- cD@data

    if(inferNulls)
      {                
        nulls <- gaps(loci)
        nulls <- nulls[seqnames(nulls) %in% unique(as.character(seqnames(cD@coordinates))),]
        nulls <- nulls[strand(nulls) == "*",]

        nullLens <- width(nulls)
#        nulls[nullLens <= 0 | is.na(nullLens),] <- NA
#        nullLens[nullLens <= 0 | is.na(nullLens)] <- NA
#        nulls <- nulls[rowSums(is.na(nulls)) == 0,1:3,drop = FALSE]
#        nullLens <- nulls$end - nulls$start + 1    
        countNulls <- getCounts(aD = aD, segments = nulls, as.matrix = TRUE, cl = cl)
        
        mD <- new("lociData", data = rbind(countNulls, countLoci),
                  seglens = c(nullLens, lociLens),
                  libsizes = as.double(aD@libsizes),
                  replicates = aD@replicates,
                  coordinates = c(nulls, loci),                  
                  locLikelihoods = rbind(matrix(-Inf, nrow = nrow(countNulls), ncol = length(levels(aD@replicates))), cD@locLikelihoods))
        
      } else {
        
        mD <- new("lociData", data = countLoci,
                  seglens = lociLens,
                  libsizes = as.double(cD@libsizes),
                  replicates = cD@replicates, 
                  coordinates = loci,
                  locLikelihoods = cD@locLikelihoods)
      }

    mD@groups <- list(mD@replicates)
    mD <- getPriors.NB(mD, verbose = TRUE, cl = cl)
    
    if(usePosteriors)
      {
        lociWeights <- sapply(levels(mD@replicates), function(rep) {
          repCol <- which(levels(mD@replicates) == rep)
          locW <- 1 - exp(mD@locLikelihoods[mD@priors$sampled[,1],repCol])
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
          replicates(repD) <- as.factor(rep(1, ncol(repD)))
          groups(repD) <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[repCol]]), list(mD@priors$priors[[1]][[repCol]]))
          repD@priors$weights <- locW
          
          repD <- getLikelihoods.NB(cD = repD, bootStraps = bootStraps, verbose = FALSE, cl = cl)
          
          message("...done!", appendLF = TRUE)
          repD@posteriors[,2]
        })
      } else {        
        lociWeights <- sapply(levels(mD@replicates), function(rep) {
          
          message(paste("Getting likelihoods for replicate group", rep), appendLF = FALSE)
          repD <- mD[,mD@replicates == rep]
          replicates(repD) <- rep(1, ncol(repD))
          groups(repD) <- list(rep(1, ncol(repD)))
          repCol <- which(levels(mD@replicates) == rep)
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[repCol]]))
          repD@priors$weights <- NULL
          
          repD <- getLikelihoods.NB(cD = repD, bootStraps = bootStraps, nullData = TRUE, verbose = FALSE, cl = cl)
          
          message("done!", appendLF = TRUE)
          repD@posteriors[,1]
        })
      }    

    mD@locLikelihoods <- lociWeights
    colnames(mD@locLikelihoods) <- levels(mD@replicates)

    #mD@annotation <- subset(mD@annotation, select = c("chr", "start", "end"))

    mD <- mD[order(as.character(seqnames(mD@coordinates)), start(mD@coordinates), end(mD@coordinates)),]

    mD
  }
