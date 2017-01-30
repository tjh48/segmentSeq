lociLikelihoods <- function(cD, aD, newCounts = FALSE, bootStraps = 3, inferNulls = TRUE, nasZero = FALSE, usePosteriors = TRUE, cl)
  {
    if(class(aD) == "alignmentMeth")
      {
          mD <- .methLikelihoods(cD = cD, aD = aD, newCounts = newCounts, bootStraps = bootStraps, inferNulls = inferNulls, usePosteriors = usePosteriors, cl = cl)
      } else if(class(aD) == "alignmentData")
            mD <- .lociLikelihoods(cD = cD, aD = aD, newCounts = newCounts, bootStraps = bootStraps, inferNulls = inferNulls, nasZero = nasZero, usePosteriors = usePosteriors, cl = cl)
    mD
  }

.inferNulls <- function(cD, aD, cl = NULL) {
    loci <- cD@coordinates
    lociLens <- width(loci)
    
    nulls <- gaps(loci)
    nulls <- nulls[seqnames(nulls) %in% unique(as.character(seqnames(cD@coordinates))),]
    nulls <- nulls[strand(nulls) == "*",]
    
    nullLens <- width(nulls)
    countNulls <- getCounts(aD = aD, segments = nulls, cl = cl)

    mD <- new("lociData", data = abind(countNulls, cD@data, along = 1),
              replicates = aD@replicates,
              coordinates = c(nulls, loci),
              seglens = c(nullLens, lociLens))
    
    if(nrow(cD@locLikelihoods) == nrow(cD))
        mD@locLikelihoods <- rbind(matrix(-Inf, length(nulls), ncol(cD@locLikelihoods)), cD@locLikelihoods)
    
    if(nrow(cD@posteriors) == nrow(cD))
        mD@posteriors <- rbind(matrix(NA, length(nulls), ncol(cD@posteriors)), cD@posteriors)
    mD@sampleObservables <- c(mD@sampleObservables)
    mD
}

.lociLikelihoods <- function(cD, aD, newCounts = FALSE, bootStraps = 3, inferNulls = TRUE, nasZero = FALSE, usePosteriors = TRUE, cl)
  {    
    loci <- cD@coordinates
    if(newCounts) cD@data <- getCounts(aD = aD, segments = loci, useChunk = TRUE, cl = cl)
    if(inferNulls)
        {
            mD <- .inferNulls(cD, aD)
        } else {
            mD <- cD
            #mD@data <- countLoci
        }

    if(all(sapply(split(replicates(mD), levels(replicates(mD))), length) == 1))
        replicates(mD) <- rep(1, ncol(mD))
    densityFunction(mD) <- nbinomDensity
    mD@groups <- list(mD@replicates)
    libsizes(mD) <- libsizes(cD)
    oldPosteriors <- mD@posteriors
    mD@posteriors <- matrix(nrow = 0, ncol = length(groups(mD)))
    mD@estProps <- numeric(0)
    
    mD <- getPriors.NB(mD, verbose = TRUE, cl = cl)
#    print("stop here!")
#    return(mD)
    
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
              list(y)
            })} else locW <- lapply(1:2, function(ii) {
              y <- locW[[ii]][[1]]
              y[is.na(y)] <- as.numeric(ii == 1)
              list(y)
            })
            
          message(paste("Getting likelihoods for replicate group", rep), appendLF = FALSE)
          repD <- mD[,which(mD@replicates == rep)]
          replicates(repD) <- as.factor(rep(1, ncol(repD)))
          groups(repD) <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[repCol]]), list(mD@priors$priors[[1]][[repCol]]))
          repD@priors$weights <- locW
          
          repD <- getLikelihoods(cD = repD, bootStraps = bootStraps, verbose = FALSE, cl = cl)
          
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
          
          repD <- getLikelihoods(cD = repD, bootStraps = bootStraps, nullData = TRUE, verbose = FALSE, cl = cl)
          
          message("done!", appendLF = TRUE)
          repD@posteriors[,1]
        })
      }    

    mD@locLikelihoods <- lociWeights
    colnames(mD@locLikelihoods) <- levels(mD@replicates)

    #mD@annotation <- subset(mD@annotation, select = c("chr", "start", "end"))

    mD <- mD[order(as.factor(seqnames(mD@coordinates)), start(mD@coordinates), end(mD@coordinates)),]

    mD@groups <- cD@groups
    mD@estProps <- cD@estProps
    #mD@posteriors <- oldPosteriors
    
    
    mD
  }


.inferMethNulls <- function(cD, aD, cl = NULL) {

    loci <- cD@coordinates
    lociLens <- width(loci)
    countLoci <- cD@data

    nulls <- gaps(loci)
    
    nulls <- nulls[seqnames(nulls) %in% unique(as.character(seqnames(cD@coordinates))),]        
    if(all(strand(loci) == "*")) {
        nulls <- nulls[strand(nulls) == "*",]
    }
    removeNull <- c(which(strand(nulls) == "*")[getOverlaps(nulls[strand(nulls) == "*",], loci, whichOverlaps = FALSE)],
                    which(strand(nulls) == "+")[getOverlaps(nulls[strand(nulls) == "+",], loci[strand(loci) %in% c("*", "+"),], whichOverlaps = FALSE)],
                    which(strand(nulls) == "-")[getOverlaps(nulls[strand(nulls) == "-",], loci[strand(loci) %in% c("*", "-"),], whichOverlaps = FALSE)])
    if(length(removeNull) > 0) nulls <- nulls[-removeNull,]
    
                                        #nulls <- nulls[strand(nulls) == "*",]
    countNulls <- getCounts(segments = nulls, aD, cl = cl)
    
    mD <- new("lociData",
              data = array(data = c(rbind(countNulls$Cs, countLoci[,,1]), rbind(countNulls$Ts, countLoci[,,2])), dim = c(nrow(countNulls$Cs) + nrow(countLoci), ncol(countLoci), 2)),
              replicates = cD@replicates,
              coordinates = c(nulls, loci),
              sampleObservables = cD@sampleObservables
              )
    
    if(nrow(cD@locLikelihoods) == nrow(cD))
        mD@locLikelihoods <- rbind(matrix(-Inf, nrow(countNulls$Cs), ncol(cD@locLikelihoods)), cD@locLikelihoods)
    
    if(nrow(cD@posteriors) == nrow(cD))
        mD@posteriors <- rbind(matrix(NA, length(nulls), ncol(cD@posteriors)), cD@posteriors)

    mD <- mD[which(rowSums(mD@data) > 0),]
    mD
}

.methLikelihoods <- function(cD, aD, newCounts = FALSE, bootStraps = 1, inferNulls = TRUE, usePosteriors = TRUE, cl)
  {
    loci <- cD@coordinates
    lociLens <- width(loci)
    countLoci <- cD@data

    if(inferNulls)
        {
            mD <- .inferMethNulls(cD, aD)
        } else mD <- cD

    mD <- mD[which(rowSums(mD@data) > 0),]
    mD@groups <- list(mD@replicates)
    oldPosteriors <- mD@posteriors
    mD@posteriors <- matrix(nrow = 0, ncol = length(groups(mD)))
    mD@estProps <- numeric(0)

    mD <- mD[order(as.factor(seqnames(mD@coordinates)), start(mD@coordinates), end(mD@coordinates)),]
    mD@data <- round(mD@data)

    if(all(mD@sampleObservables$nonconversion == 0)) {
        densityFunction(mD) <- bbDensity
        libsizes(mD) <- matrix(1, ncol = 2, nrow = ncol(mD))
    } else {
        densityFunction(mD) <- bbNCDist
        mD <- methObservables(mD)
    }
    
    mD <- getPriors(mD, samplesize = 1e4, verbose = TRUE, cl = cl)        
    

    if(usePosteriors)
      {
        lociWeights <- sapply(levels(mD@replicates), function(rep) {
          repCol <- which(levels(mD@replicates) == rep)
          locW <- 1 - exp(mD@locLikelihoods[mD@priors$sampled[,1],repCol])
          locW <- list(list(locW), list(1 - locW))
          
          locW <- lapply(locW, function(x) {
            y <- x[[1]]
            y[is.na(y)] <- 0
            list(y)
          })
            
          message(paste("Getting likelihoods for replicate group", rep), appendLF = FALSE)
          repD <- mD[,mD@replicates == rep]
          replicates(repD) <- as.factor(rep(1, ncol(repD)))
          groups(repD) <- list(rep(1, ncol(repD)), rep(1, ncol(repD)))
          repD@priors$priors <- list(list(mD@priors$priors[[1]][[repCol]]), list(mD@priors$priors[[1]][[repCol]]))
          repD@priors$weights <- locW
          
          repD <- getLikelihoods(cD = repD, bootStraps = bootStraps, verbose = FALSE, cl = cl)
          
          
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
          
          repD <- getLikelihoods(cD = repD, bootStraps = bootStraps, nullData = TRUE, verbose = FALSE, cl = cl)
          
          message("done!", appendLF = TRUE)
          repD@posteriors[,1]
        })
      }    

    mD@locLikelihoods <- lociWeights
    colnames(mD@locLikelihoods) <- levels(mD@replicates)

    #mD@annotation <- subset(mD@annotation, select = c("chr", "start", "end"))    

    mD@groups <- cD@groups
    mD@estProps <- cD@estProps
    mD@posteriors <- oldPosteriors
    
    mD
  }
