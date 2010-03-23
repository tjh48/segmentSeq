getPriors <- function(sD, type = "Pois", verbose = TRUE, ...)
  {
    chrs <- unique(sD@segInfo$chr)

    if(verbose) message("Finding priors...")
    
    reppriors <- lapply(unique(sD@replicates), function(rr, ...)
                        {
                          getSubPriors <- function(subD, type, ...)
                            { 
                              filsegs <- filterSegments(subset(subD@segInfo, select = c(chr, start, end)), runif(nrow(subD)))
                              
                              cD <- new("countData", data = subD@data[filsegs,, drop = FALSE], seglens = (subD@segInfo$end - subD@segInfo$start + 1L)[filsegs],
                                        libsizes = subD@libsizes, groups = list(rep(1, ncol(subD))))
                              
                              cD <- switch(type,
                                           Pois = (getPriors.Pois(cD, verbose = FALSE, ...)),
                                           NB = (getPriors.NB(cD, verbose = FALSE, ...)))
                              
                              list(priors = cD@priors, priorType = cD@priorType)
                            }     

                          if(verbose) message(paste("Replicate group: ", rr, "...", sep = ""), appendLF = FALSE)
                          
                          rsD <- sD[,sD@replicates == rr]
                          
                          sumD <- colSums(t(rsD@data) / rsD@libsizes)  / (sD@segInfo$end - sD@segInfo$start + 1)
                          
                          NZs <- sumD != 0
                          NZD <- log(sumD[NZs])
                          
                          dsortp <- sort(NZD, decreasing = TRUE)
                          dcummeans <- cumsum(dsortp) / 1:length(dsortp)
                          dcumsquare <- cumsum(dsortp ^ 2) / (1:length(dsortp) - 1)
                          dcumvar <- dcumsquare - (dcummeans^2) * 1:length(dsortp) / (1:length(dsortp) - 1)
                          dcumse <- sqrt(dcumvar / 1:length(dsortp))
                          
                          isortp <- sort(NZD, decreasing = FALSE)
                          icummeans <- cumsum(isortp) / 1:length(isortp)
                          icumsquare <- cumsum(isortp ^ 2) / (1:length(isortp) - 1)
                          icumvar <- icumsquare - (icummeans^2) * 1:length(isortp) / (1:length(isortp) - 1)
                          icumse <- sqrt(icumvar / 1:length(isortp))
                          
                          meanvars <- (c(rev(dcumvar)[-1], NA) + icumvar) / 2
                          
                          varrise <- which(meanvars[-1] > meanvars[-length(meanvars)])
                          sep <- min(varrise[varrise > length(NZD) / 100])
        
                          kstar <- c(mean(isortp[1:sep]), mean(isortp[(sep + 1):length(isortp)]))
                          
                          kmNZ <- kmeans(NZD, centers = matrix(kstar, ncol = 1))
                          km <- rep(1, length(sumD))
                          km[NZs] <- kmNZ$cluster
                          
                          nulPriors <- getSubPriors(rsD[km == 1,], type, ...)
                          locPriors <- getSubPriors(rsD[km == 2,], type, ...)
                          
                          priors <-
                            list(priors = list(nulPriors$priors$priors[[1]][[1]], locPriors$priors$priors[[1]][[1]]), priorType = nulPriors$priorType)
                          if(verbose) message(paste("done."))
                          
                          priors
                        }, ...)

    sD@priors <- list(lapply(reppriors, function(x) x$priors[[1]]), lapply(reppriors, function(x) x$priors[[2]]))
    sD@priorType <- unique(sapply(reppriors, function(x) x$priorType))

    sD
  }
