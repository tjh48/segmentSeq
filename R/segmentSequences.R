segmentSequences <-
function(sDP,
                             pcut = 0.5, estimatePriors = FALSE, priorDE = 1e-2, verbose = TRUE, ...,
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

        if(nrow(sDP@leftData) == 0) {
          sideData <- matrix(0, ncol = sum(replicates == rep), nrow = nrow(data))
          smallerThan <- TRUE
        } else {
          sideData <- sDP@leftData[,replicates == rep]
          smallerThan <- colSums(t(sideData / sDP@segInfo$leftSpace) / sDP@libsizes[replicates == rep]) <= colSums(t(sDP@data[,replicates == rep] / seglens) / sDP@libsizes[replicates == rep])
        }
        
        countL <- new("countData",
                      data = cbind(sideData, data[,replicates == rep, drop = FALSE]),
                      seglens = cbind(matrix(sDP@segInfo$leftSpace, ncol = sum(replicates == rep), nrow = nrow(data)),
                        matrix(seglens, ncol = sum(replicates == rep), nrow = nrow(data))),
                      libsizes = rep(sDP@libsizes[replicates == rep], 2),
                      groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(NDE = list(priors[[rep]]), DE = list(priors[[rep]], priors[[rep]])), sampled = sDP@priors$sampled))
        leftLike <- getLikelihoods(countL, prs = c(1-priorDE, priorDE), pET = "none",
                                   subset = which(smallerThan & rowSums(countL@seglens == 0) == 0),
                                   verbose = FALSE, ..., cl = cl)
        if(any(!smallerThan))
          {
            leftLike@posteriors[!smallerThan,1] <- 0
            leftLike@posteriors[!smallerThan,2] <- -Inf
          }

        
        if(verbose) message(".", appendLF = FALSE)

        if(nrow(sDP@rightData) == 0) {
          sideData <- matrix(0, ncol = sum(replicates == rep), nrow = nrow(data))
          smallerThan <- TRUE
        } else {
          sideData <- sDP@rightData[,replicates == rep]
          smallerThan <- colSums(t(sideData / sDP@segInfo$rightSpace) / sDP@libsizes[replicates == rep]) <= colSums(t(sDP@data[,replicates == rep] / seglens) / sDP@libsizes[replicates == rep])
        }
        
        countR <- new("countData",
                      data = cbind(sideData, data[,replicates == rep, drop = FALSE]),
                      seglens = cbind(matrix(sDP@segInfo$rightSpace, ncol = sum(replicates == rep), nrow = nrow(data)),
                        matrix(seglens, ncol = sum(replicates == rep), nrow = nrow(data))),
                      libsizes = rep(sDP@libsizes[replicates == rep], 2),
                      groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(NDE = list(priors[[rep]]), DE = list(priors[[rep]], priors[[rep]])), sampled = sDP@priors$sampled))
        rightLike <- getLikelihoods(countR, prs = c(1-priorDE, priorDE), pET = "none",
                                    subset = which(smallerThan & rowSums(countR@seglens == 0) == 0),
                                    verbose = FALSE, ..., cl = cl)
        if(any(!smallerThan))
          {
            rightLike@posteriors[!smallerThan,1] <- 0
            rightLike@posteriors[!smallerThan,2] <- -Inf
          }

        if(verbose) message(".", appendLF = FALSE)
        
        countN <- new("countData",
                      data = data[,replicates == rep, drop = FALSE],
                      seglens = seglens,
                      libsizes = sDP@libsizes[replicates == rep],
                      groups = list(c(rep(1, sum(replicates == rep))), c(rep(1, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(null = list(priors[[rep]]), loc = list(priors[[rep]])), sampled = sDP@priors$sampled))
        nullLike <- getLikelihoods(countN, prs = c(0.1, 0.9), pET = "none", verbose = FALSE, ..., cl = cl)
        
        if(verbose) message(".")

        if(is.null(cl))
          {
            sumSame <- apply(cbind(nullLike@posteriors[,1L], leftLike@posteriors[,1L], rightLike@posteriors[,1L], nullLike@posteriors[,1L] + leftLike@posteriors[,1L] + rightLike@posteriors[,1L]), 1, logsum)
            minSame <- apply(cbind(nullLike@posteriors[,1L] + leftLike@posteriors[,1L], nullLike@posteriors[,1L] + rightLike@posteriors[,1L], leftLike@posteriors[,1L] + rightLike@posteriors[,1L]), 1, logsum)
          } else {
            sumSame <- parRapply(cl = cl, cbind(nullLike@posteriors[,1L], leftLike@posteriors[,1L], rightLike@posteriors[,1L], nullLike@posteriors[,1L] + leftLike@posteriors[,1L] + rightLike@posteriors[,1L]), logsum)
            minSame <- parRapply(cl = cl, cbind(nullLike@posteriors[,1L] + leftLike@posteriors[,1L], nullLike@posteriors[,1L] + rightLike@posteriors[,1L], leftLike@posteriors[,1L] + rightLike@posteriors[,1L]), logsum)
          }
        
        PSames <- sumSame + log(1 - exp(minSame - sumSame))
        PSames
      }

    if(verbose) message("Evaluating likelihoods of similarity for each segment...")
    
    pSames <- sapply(unique(sDP@replicates), probSames, sDP = sDP, ...)
    sSames <- rowSums(pSames)

    if(verbose) message("Filtering loci...", appendLF = FALSE)

    ssDP <- sDP[exp(sSames) < pcut,]
    sSames <- sSames[exp(sSames) < pcut]

    if(nrow(sDP) > 0) {
      filsegs <- filterSegments(segs = subset(ssDP@segInfo, select = c(chr, start, end)), orderOn = sSames, decreasing = FALSE)
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

