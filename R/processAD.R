
processAD <-
function(aD, maxgaplen = 500, maxloclen = NULL, verbose = TRUE, cl = cl)
  {    
    cTags <- aD@alignments
    cTags$tag <- as.numeric(as.factor(cTags$tag))
    chrs <- aD@chrs
    chrlens <- aD@chrlens
    libnames <- aD@libnames
    libsizes <- aD@libsizes
    replicates <- aD@replicates
    tagData <- aD@data

    segInfo <- NULL
    data <- NULL

    if(is.null(maxloclen)) maxloclen <- Inf
    
    for(cc in 1:length(chrs))
      {
        if(verbose){
          message("Chromosome: ", chrs[cc])
          message("Finding start-stop co-ordinates...")
        }

        if(any(cTags$chr == chrs[cc]))
          {
            chrTags <- cTags[cTags$chr == chrs[cc],, drop = FALSE]
            chrTagData <- tagData[cTags$chr == chrs[cc],, drop = FALSE]
            
            if(any(chrTags[,3L] > chrlens[cc]))
              warning(paste("Chromosome", chrs[cc], "has tags which extend over the given chromsome length."))
            
            chrSS <- chrTags[, 2:3, drop = FALSE]
            if(nrow(chrSS) > 1)
              {
                chrSS <- chrSS[order(chrSS[,1L], chrSS[,2L]),]
                chrmax <- apply(chrSS, 2, cummax)
                ch <- which(chrmax[-1L,1L] > chrmax[-nrow(chrmax),2L] + 1)
                startstop <- cbind(starts = as.numeric(c(min(chrmax[,1L]), chrmax[ch + 1L,1L])),
                                   ends = as.numeric(c(chrmax[ch,2L], max(chrmax[,2L]))))
              } else if(nrow(chrSS) == 1) startstop <- cbind(starts = chrSS$start, ends = chrSS$end) else startstop <- matrix(nrow = 0, ncol = 0)
            
            gaps <- c(which(startstop[-nrow(startstop),2L] < startstop[-1L, 1L] - maxgaplen), nrow(startstop))

            if(verbose)
              message("Defining potential subsegments...", appendLF = FALSE)
            gapped <- lapply(1:nrow(startstop), function(x) x:min(gaps[gaps >= x]))
            csegs <- cbind(start = rep(startstop[,1L], lapply(gapped, length)), end = unlist(lapply(gapped, function(x) startstop[x,2L])))

            if(any(csegs[,2L] - csegs[,1L] + 1L > maxloclen))
              {
                rsegs <- csegs[csegs[,2L] - csegs[,1L] + 1 <= maxloclen,]
                rsegs <- rbind(rsegs, startstop[(!(startstop[,1L] %in% rsegs[,1L])),])
              } else rsegs <- csegs
            
            rsegs <- rsegs[order(rsegs[,1L], rsegs[,2L]),,drop = FALSE]

            if(verbose){
              message(".", nrow(rsegs), " found.")
              message("Getting count data for each potential subsegment...", appendLF = FALSE)
            }
            
            winsize <- 1e5
            
            countSegs <- matrix(unlist(lapply(1:ceiling(nrow(rsegs) / winsize), function(rr)
                                {
                                  winSegs <- rsegs[((rr - 1) * winsize + 1):min(rr * winsize, nrow(rsegs)),,drop = FALSE]
                                  seltags <- which(chrTags$start <= max(winSegs) & chrTags$end >= min(winSegs))
                                  winTags <- chrTags[seltags,,drop = FALSE]
                                  winTagData <- chrTagData[seltags,,drop = FALSE]
                                  
                                  t(getCounts(segments = data.frame(chr = chrs[cc], start = winSegs[,1L], end = winSegs[,2L]),
                                            aD = new("alignmentData", libnames = libnames, libsizes = libsizes, alignments = winTags, data = winTagData, chrs = chrs[cc], chrlens = chrlens[cc], replicates = replicates),                                      cl = cl))
                                })), nrow = nrow(rsegs), byrow = TRUE)
            
            if(verbose)
              message("done!")

            rsegs <- cbind(rsegs,
                           leftSpace = (rsegs[,1L] - c(0L, startstop[,2])[findInterval(rsegs[,1L], startstop[,2L]) + 1L]) - 1L,
                           rightSpace = c(startstop[,1L], chrlens[cc] + 1L)[findInterval(rsegs[,2L], startstop[,1L]) + 1L] - rsegs[,2L] - 1L)
            rownames(rsegs) <- rep("", nrow(rsegs))
            segInfo <- rbind(segInfo, data.frame(chr = I(chrs[cc]), rsegs))
            data <- rbind(data, countSegs)
          } else if(verbose) message("No tags found for this chromosome.")
      }

    colnames(data) <- aD@libnames
    
    new("segData", segInfo = segInfo, data = data, libsizes = libsizes, replicates = aD@replicates)
  }

