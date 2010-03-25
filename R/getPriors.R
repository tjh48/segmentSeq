getPriors <- function(sD, type = "Pois", verbose = TRUE, ...)
  {
    if(verbose) message("Finding priors...")

    reppriors <- lapply(unique(sD@replicates), function(rep,...)
                        {
                          if(verbose) message(paste("Replicate group: ", rep, "...", sep = ""), appendLF = FALSE)
                          
                          filsegs <- filterSegments(subset(sD@segInfo, select = c(chr, start, end)), runif(nrow(sD)))
                          
                          cD <- new("countData",
                                    data = sD@data[filsegs, sD@replicates == rep, drop = FALSE],
                                    seglens = sD@segInfo$end[filsegs] - sD@segInfo$start[filsegs] + 1L,
                                    libsizes = sD@libsizes[sD@replicates == rep],
                                    groups = list(rep(1, sum(sD@replicates == rep))))
                          
                          cD <- switch(type,
                                       Pois = (getPriors.Pois(cD, verbose = FALSE, ...)),
                                       NB = (getPriors.NB(cD, verbose = FALSE, ...)))
                          
                              message("done.")
                          
                          list(priors = cD@priors$priors[[1]], priorType = cD@priorType)
                        }, ...)        
                      
    sD@priors <- list(lapply(reppriors, function(x) x$priors[[1]]))
    sD@priorType <- unique(sapply(reppriors, function(x) x$priorType))

    sD
  }
