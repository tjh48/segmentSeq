getPriors <- function(sD, type = "Pois", verbose = TRUE, ...)
  {
    chrs <- unique(sD@segInfo$chr)

    if(verbose) message("Constructing sample of non-overlapping segments...", appendLF = FALSE)
    
    filsegs <- filterSegments(subset(sD@segInfo, select = c(chr, start, end)), runif(nrow(sD)))

    if(verbose) message("done.")
    
    cD <- new("countData", data = sD@data[filsegs,, drop = FALSE], seglens = (sD@segInfo$end - sD@segInfo$start + 1L)[filsegs],
              libsizes = sD@libsizes, groups = list(sD@replicates))

    cD <- switch(type,
                 Pois = (getPriors.Pois(cD, ...)),
                 NB = (getPriors.NB(cD, ...)))
           
    sD@priors <- cD@priors
    sD@priorType <- cD@priorType

    sD
  }
