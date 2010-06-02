segmentSequences <-
function(sDP,
                             pcut = 0.5, priorDE = 1e-2, verbose = TRUE, topOnly = FALSE, ...,
                             cl)
  {
    logsum <- function(x) {
      if(all(is.na(x))) return(NA)
      if(all(x == -Inf)) return(-Inf)
      y <- x[!is.na(x)]
      max(y, max(y) + log(sum(exp(y - max(y)))))
    }

    if(!is.null(cl))
      {
        segmentEnv <- new.env(parent = .GlobalEnv)
        environment(logsum) <- segmentEnv
      }
    seglens <- sDP@segInfo$end - sDP@segInfo$start + 1
    
    probSames <- function(rep, sDP, ...)
      {
        priorType <- sDP@priorType
        priors <- sDP@priors
        data <- sDP@data
        replicates <- sDP@replicates
        
        chrs <- unique(sDP@segInfo$chr)
        
        if(verbose) message("Replicate group: ", rep, appendLF = FALSE)

        sideSame <- function(sideData, sideSpace)
          {
            if(nrow(sideData) == 0) {
              sideData <- matrix(0, ncol = sum(replicates == rep), nrow = nrow(data))
              smallerThan <- TRUE
            } else {
              sideData <- sideData[,replicates == rep,drop = FALSE]
              smallerThan <- colSums(t(sideData / sideSpace) / sDP@libsizes[replicates == rep]) <= colSums(t(sDP@data[,replicates == rep, drop = FALSE] / seglens) / sDP@libsizes[replicates == rep])
            }
            
            countSide <- new("countData",
                             data = cbind(sideData, data[,replicates == rep, drop = FALSE]),
                             seglens = cbind(matrix(sideSpace, ncol = sum(replicates == rep), nrow = nrow(data)),
                               matrix(seglens, ncol = sum(replicates == rep), nrow = nrow(data))),
                             libsizes = rep(sDP@libsizes[replicates == rep], 2),
                             groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                             priorType = priorType,
                             priors = list(priors = list(NDE = list(priors[[rep]]), DE = list(priors[[rep]], priors[[rep]])), sampled = sDP@priors$sampled))
            if(any(smallerThan & rowSums(countSide@seglens == 0) == 0))
              {
                countSide <- getLikelihoods(countSide, prs = c(1-priorDE, priorDE), pET = "none",
                                            subset = which(smallerThan & rowSums(countSide@seglens == 0) == 0),
                                            verbose = FALSE, cl = cl)
              } else countSide@posteriors <- matrix(c(0, -Inf), nrow = nrow(countSide@data), ncol = 2, byrow = TRUE)
            
            if(any(!smallerThan[!is.na(smallerThan)]))
              {
                countSide@posteriors[!smallerThan,1] <- 0
                countSide@posteriors[!smallerThan,2] <- -Inf
              }
            
            if(verbose) message(".", appendLF = FALSE)

            countSide@posteriors
          }

        countL <- sideSame(sDP@leftData, sDP@segInfo$leftSpace)
        countR <- sideSame(sDP@rightData, sDP@segInfo$rightSpace)
        
        if(verbose) message(".", appendLF = FALSE)
        
        countN <- new("countData",
                      data = cbind(matrix(0, ncol = sum(replicates == rep), nrow = nrow(data)), data[,replicates == rep, drop = FALSE]),
                      seglens = seglens,
                      libsizes = rep(sDP@libsizes[replicates == rep], 2),
                      groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(NDE = list(priors[[rep]]), DE = list(priors[[rep]], priors[[rep]])), sampled = sDP@priors$sampled))
        nullLike <- getLikelihoods(countN, prs = c(0.9, 0.1), pET = "none", verbose = FALSE, cl = cl)@posteriors
        
        if(verbose) message(".")

        if(is.null(cl))
          {
            sumSame <- apply(cbind(nullLike[,1L], countL[,1L], countL[,1L], nullLike[,1L] + countL[,1L] + countR[,1L]), 1, logsum)
            minSame <- apply(cbind(nullLike[,1L] + countL[,1L], nullLike[,1L] + countR[,1L], countL[,1L] + countR[,1L]), 1, logsum)
          } else {
            sumSame <- parRapply(cl = cl, cbind(nullLike[,1L], countL[,1L], countR[,1L], nullLike[,1L] + countL[,1L] + countR[,1L]), logsum)
            minSame <- parRapply(cl = cl, cbind(nullLike[,1L] + countL[,1L], nullLike[,1L] + countR[,1L], countL[,1L] + countR[,1L]), logsum)
          }
        
        PSames <- sumSame + log(1 - exp(minSame - sumSame))
        PSames
      }

    if(verbose) message("Evaluating likelihoods of similarity for each segment...")
    
    pSames <- matrix(unlist(lapply(unique(sDP@replicates), probSames, sDP = sDP)), nrow = nrow(sDP@data), byrow = FALSE)
    sSames <- rowSums(pSames)
      
    
    if(verbose) message("Filtering loci...", appendLF = FALSE)

    ssDP <- sDP[exp(sSames) < pcut,]
    sSames <- sSames[exp(sSames) < pcut]
      
    
    if(nrow(ssDP) > 0) {
      if(topOnly) filsegs <- which.min(sSames) else filsegs <- filterSegments(segs = subset(ssDP@segInfo, select = c(chr, start, end)), orderOn = sSames, decreasing = FALSE)
    } else filsegs <- NULL
    
    if(verbose) message("done!")
    
    ssDP <- ssDP[filsegs,]

    if(length(filsegs) > 0)
      {
        sD <- new("countData", data = ssDP@data, replicates = ssDP@replicates, annotation = cbind(subset(ssDP@segInfo, select = c(chr, start, end)), PSame = exp(sSames[filsegs])), seglens = seglens[filsegs], libsizes = ssDP@libsizes)
        sD <- sD[order(sD@annotation$chr, sD@annotation$start, sD@annotation$end, decreasing = FALSE),]
      } else sD <- new("countData")
    
    sD
    
  }

