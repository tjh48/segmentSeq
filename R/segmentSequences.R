segmentSequences <-
function(sDP,
                             pcut = 0.5, estimatePriors = FALSE, verbose = TRUE, ...,
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

    priorType <- sDP@priorType
    priors <- sDP@priors
    data <- sDP@data
    replicates <- sDP@replicates

    seglens <- sDP@segInfo$end - sDP@segInfo$start + 1

    chrs <- unique(sDP@segInfo$chr)

    probSames <- function(rep,...)
      {
        if(verbose) message("Replicate group: ", rep, appendLF = FALSE)
        
        countL <- new("countData",
                      data = cbind(matrix(0, ncol = sum(replicates == rep), nrow = nrow(data)), data[,replicates == rep, drop = FALSE]),
                      seglens = cbind(matrix(sDP@segInfo$leftSpace, ncol = sum(replicates == rep), nrow = nrow(data)),
                        matrix(seglens, ncol = sum(replicates == rep), nrow = nrow(data))),
                      libsizes = rep(sDP@libsizes[replicates == rep], 2),
                      groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(NDE = list(priors[[2]][[rep]]), DE = list(priors[[2]][[rep]], priors[[1]][[rep]])), sampled = sDP@priors$sampled))
        leftLike <- getLikelihoods(countL, prs = c(0.9, 0.1), pET = "none", subset = which(rowSums(countL@seglens != 0) == ncol(countL@seglens)), verbose = FALSE, ..., cl = cl)
                
        if(verbose) message(".", appendLF = FALSE)
        
        countR <- new("countData",
                      data = cbind(matrix(0, ncol = sum(replicates == rep), nrow = nrow(data)), data[,replicates == rep, drop = FALSE]),
                      seglens = cbind(matrix(sDP@segInfo$rightSpace, ncol = sum(replicates == rep), nrow = nrow(data)),
                        matrix(seglens, ncol = sum(replicates == rep), nrow = nrow(data))),
                      libsizes = rep(sDP@libsizes[replicates == rep], 2),
                      groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(NDE = list(priors[[2]][[rep]]), DE = list(priors[[2]][[rep]], priors[[1]][[rep]])), sampled = sDP@priors$sampled))
        rightLike <- getLikelihoods(countR, prs = c(0.9, 0.1), pET = "none", subset = which(rowSums(countR@seglens != 0) == ncol(countR@seglens)), verbose = FALSE, ..., cl = cl)

        if(verbose) message(".", appendLF = FALSE)
                
        countN <- new("countData",
                      data = data[,replicates == rep, drop = FALSE],
                      seglens = seglens,
                      libsizes = sDP@libsizes[replicates == rep],
                      groups = list(c(rep(1, sum(replicates == rep))), c(rep(1, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(null = list(priors[[1]][[rep]]), loc = list(priors[[2]][[rep]])), sampled = sDP@priors$sampled))
        nullLike <- getLikelihoods(countN, prs = c(0.1, 0.9), pET = "none", verbose = FALSE, ..., cl = cl)

        if(verbose) message(".")
        
        if(estimatePriors)
          {
            filleft <- filterSegments(segs = subset(sDP@segInfo, select = c(chr, start, end)), orderOn = leftLike@posteriors[,1], decreasing = FALSE)
            leftLike <- getLikelihoods(countL, prs = colSums(exp(leftLike@posteriors), na.rm = TRUE) / sum(exp(leftLike@posteriors), na.rm = TRUE), pET = "BIC", subset = which(rowSums(countL@seglens != 0) == ncol(countL@seglens)), priorSubset = filleft, verbose = FALSE, ..., cl = cl)
            filright <- filterSegments(segs = subset(sDP@segInfo, select = c(chr, start, end)), orderOn = rightLike@posteriors[,1], decreasing = FALSE)
            rightLike <- getLikelihoods(countR, prs = colSums(exp(rightLike@posteriors), na.rm = TRUE) / sum(exp(rightLike@posteriors), na.rm = TRUE), pET = "BIC", subset = which(rowSums(countR@seglens != 0) == ncol(countR@seglens)), priorSubset = filright, verbose = FALSE, ..., cl = cl)
            filnull <- filterSegments(segs = subset(sDP@segInfo, select = c(chr, start, end)), orderOn = nullLike@posteriors[,1], decreasing = FALSE)
            nullLike <- getLikelihoods(countN, prs = colSums(exp(nullLike@posteriors), na.rm = TRUE) / sum(exp(nullLike@posteriors), na.rm = TRUE), pET = "BIC", subset = which(rowSums(countN@seglens != 0) == ncol(countN@seglens)), priorSubset = filnull, verbose = FALSE, ..., cl = cl)
          }
        
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
    
    pSames <- sapply(unique(replicates), probSames, ...)
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

