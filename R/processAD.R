.windowing <- function(cs, winsize = 5e5)
  {
    chunkWindows <- cbind(1:length(runLength(cs$csegChunk)), cumsum(runLength(cs$csegChunk)))
    
    dupWin <- which(!duplicated(round(chunkWindows[,2] / winsize)))
    dupWin <- cbind(dupWin, c(dupWin[-1] - 1, nrow(chunkWindows)))
    windowChunks <- lapply(1:nrow(dupWin), function(ii) unique(cs$csegChunk)[chunkWindows[dupWin[ii,1]:dupWin[ii,2],1L]])
    windowChunks
  }

.chrProcessing <- function(cTags, cc, verbose)
  {
    chrSS <- cTags[cTags@seqnames == cc,]
    
    if(!is.na(seqlengths(cTags)[seqlevels(cTags) == cc]) && any(end(chrSS) > seqlengths(cTags)[seqlevels(cTags) == cc]))
      warning(paste("Chromosome", cc, "has tags which extend over the given chromsome length."))
    
    if(verbose){
      message("Chromosome: ", cc)
      message("Finding start-stop co-ordinates...", appendLF = FALSE)
    }
    
    startstop <- reduce(chrSS@ranges, min.gapwidth = 0)
    ssChunk <- values(chrSS)$chunk[match(start(startstop), start(chrSS))]
    chunkDiff <- runLength(ssChunk)
    
    startRep <- (rep(cumsum(chunkDiff), chunkDiff) - 1:length(startstop) + 1)
    endRep <- cbind(1:length(startstop), rep(cumsum(chunkDiff) + 1, chunkDiff) - 1)
    
    csegs <- IRanges(start = rep(start(startstop), startRep),
                     end = end(startstop)[unlist(lapply(1:nrow(endRep), function(ii) endRep[ii,1]:endRep[ii,2]))])
    
    csegChunk <- Rle(as.integer(rep(ssChunk, startRep)))

    list(csegs = csegs, csegChunk = csegChunk)
  }


processAD <-
function(aD, gap = NULL, verbose = TRUE, cl)
  {
    cTags <- aD@alignments

    if("tag" %in% colnames(values(cTags))) {      
      values(cTags)$tag <- as.integer(as.factor(values(cTags)$tag))
    } else values(cTags)$tag <- 1:length(cTags)

    if(!is.null(gap))
      cTags <- findChunks(cTags, gap)

    coordinates <- GRanges()
    Cs <- Ts <- data <- DataFrame()

    for(cc in seqlevels(cTags))
      {
        if(any(cTags@seqnames == cc))
          {
            cs <- .chrProcessing(cTags, cc, verbose = verbose)
            
            if(verbose){
              message(".", length(cs$csegs), " found.")
              message("Getting count data for each potential subsegment...", appendLF = FALSE)
            }

            windowChunks <- .windowing(cs)
            
            waD <- new("alignmentData")            
            waD@libnames <- aD@libnames            
            waD@replicates = aD@replicates
            
            partCounts <- function(chunks, data) {
              chad <- which(as.integer(values(cTags)$chunk) %in% chunks)
              waD@alignments <- cTags[chad,,drop = FALSE]
              waD@data <- data[chad,,drop = FALSE]
              message(".", appendLF = FALSE)                                              
              segments <- GRanges(seqnames = cc, cs$csegs[which(as.integer(cs$csegChunk) %in% chunks),])
              getCounts(segments = segments, 
                        aD = waD, preFiltered = FALSE, cl = cl)
            }

            if(class(aD) == "alignmentData") {
              data <- rbind(data, do.call("rbind",
                                          lapply(windowChunks, partCounts, data = aD@data)))
            } else if(class(aD) == "methAlignment") {
              Cs <- rbind(Cs, do.call("rbind",
                                      lapply(windowChunks, partCounts, data = aD@Cs)))                                     
              Ts <- rbind(Ts, do.call("rbind",
                                      lapply(windowChunks, partCounts, data = aD@Ts)))
            }
            
            if(verbose)
              message("...done!")
            
            coordinates <- suppressWarnings(c(coordinates, GRanges(seqnames = cc, cs$csegs)))
            
            gc()

          } else if(verbose) message("No tags found for this chromosome.")
      }

    if(class(aD) == "alignmentData") {
      tD <- new("segData", coordinates = coordinates, data = data, libsizes = aD@libsizes, replicates = aD@replicates)
      seqlengths(tD@coordinates) <- seqlengths(aD@alignments)
      colnames(tD@data) <- aD@libnames
    } else if(class(aD) == "methAlignment") {
      tD <- new("methSegs", coordinates = coordinates, Cs = Cs, Ts = Ts, replicates = aD@replicates)
      colnames(tD@Cs) <- colnames(tD@Ts) <- aD@libnames
    }

    tD
  }
