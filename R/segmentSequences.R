segmentSequences <-
function(sDP,
                             pcut = 0.5, estimatePriors = FALSE, verbose = TRUE, priorDE = 1e-2, ...,
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
                      priors = list(priors = list(NDE = list(priors[[1]][[rep]]), DE = list(priors[[1]][[rep]], priors[[1]][[rep]])), sampled = sDP@priors$sampled))
        leftLike <- getLikelihoods(countL, prs = c(1-priorDE, priorDE), pET = "none", subset = which(rowSums(countL@seglens != 0) == ncol(countL@seglens)), verbose = FALSE, ..., cl = cl)
                
        if(verbose) message(".", appendLF = FALSE)

        countR <- new("countData",
                      data = cbind(matrix(0, ncol = sum(replicates == rep), nrow = nrow(data)), data[,replicates == rep, drop = FALSE]),
                      seglens = cbind(matrix(sDP@segInfo$rightSpace, ncol = sum(replicates == rep), nrow = nrow(data)),
                        matrix(seglens, ncol = sum(replicates == rep), nrow = nrow(data))),
                      libsizes = rep(sDP@libsizes[replicates == rep], 2),
                      groups = list(c(rep(1, sum(replicates == rep) * 2)), c(rep(1, sum(replicates == rep)), rep(2, sum(replicates == rep)))),
                      priorType = priorType,
                      priors = list(priors = list(NDE = list(priors[[1]][[rep]]), DE = list(priors[[1]][[rep]], priors[[1]][[rep]])), sampled = sDP@priors$sampled))
        rightLike <- getLikelihoods(countR, prs = c(1-priorDE, priorDE), pET = "none", subset = which(rowSums(countR@seglens != 0) == ncol(countR@seglens)), verbose = FALSE, ..., cl = cl)

        if(verbose) message(".")
        
        if(is.null(cl))
          {
            sumSame <- apply(cbind(leftLike@posteriors[,1L], rightLike@posteriors[,1L]), 1, logsum)
          } else {
            sumSame <- parRapply(cl = cl, cbind(leftLike@posteriors[,1L], rightLike@posteriors[,1L]), logsum)
          }

        minSame <- leftLike@posteriors[,1L] + rightLike@posteriors[,1L]
        PSames <- sumSame + log(1 - exp(minSame - sumSame))
        PSames[is.na(leftLike@posteriors[,1L])] <- rightLike@posteriors[is.na(leftLike@posteriors[,1L]),1L]
        PSames[is.na(rightLike@posteriors[,1L])] <- leftLike@posteriors[is.na(rightLike@posteriors[,1L]),1L]
        
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

