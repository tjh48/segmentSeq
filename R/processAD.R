processAD <-
function(aD, gap = NULL, verbose = TRUE, cl)
  {
    cTags <- aD@alignments

    if("tag" %in% colnames(values(cTags))) {      
      values(cTags)$tag <- as.integer(as.factor(values(cTags)$tag))
    } else values(cTags)$tag <- 1:length(cTags)

    if(!is.null(gap))
      cTags <- findChunks(cTags, gap)

    libnames <- aD@libnames
    libsizes <- aD@libsizes
    replicates <- aD@replicates
    tagData <- aD@data

    coordinates <- GRanges()
    data <- DataFrame()

    for(cc in seqlevels(cTags))
      {
        if(any(cTags@seqnames == cc))
          {
            chrSS <- cTags[cTags@seqnames == cc,]
            
            if(any(end(chrSS) > seqlengths(cTags)[seqlevels(cTags) == cc]))
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

            if(verbose){
              message(".", length(csegs), " found.")
              message("Getting count data for each potential subsegment...", appendLF = FALSE)
            }
            
            winsize <- 5e5
            chunkWindows <- cbind(1:length(runLength(csegChunk)), cumsum(runLength(csegChunk)))

            dupWin <- which(!duplicated(round(chunkWindows[,2] / winsize)))
            dupWin <- cbind(dupWin, c(dupWin[-1] - 1, nrow(chunkWindows)))
            windowChunks <- lapply(1:nrow(dupWin), function(ii) unique(csegChunk)[chunkWindows[dupWin[ii,1]:dupWin[ii,2],1L]])

            waD <- new("alignmentData")
            waD@libnames <- libnames
            waD@libsizes <- libsizes
            waD@replicates = replicates

            gc()

            data <- rbind(data, do.call("rbind",
                                        lapply(windowChunks, function(chunks) {
                                          chad <- which(as.integer(values(cTags)$chunk) %in% chunks)
                                          
                                          waD@alignments <- cTags[chad,,drop = FALSE]
                                          waD@data <- aD@data[chad,,drop = FALSE]
                                          message(".", appendLF = FALSE)
                                          
                                          segments <- GRanges(seqnames = cc, csegs[which(as.integer(csegChunk) %in% chunks),])
                                          
                                          getCounts(
                                                    segments = segments, 
                                                    aD = waD, preFiltered = FALSE, cl = cl)
                                        })
                                        ))

            if(verbose)
              message("...done!")

            coordinates <- suppressWarnings(c(coordinates, GRanges(seqnames = cc, csegs)))
            
            gc()

          } else if(verbose) message("No tags found for this chromosome.")
      }

    tD <- new("segData", coordinates = coordinates, data = data, libsizes = libsizes, replicates = replicates)

    seqlengths(tD@coordinates) <- seqlengths(aD@alignments)
    colnames(tD@data) <- aD@libnames

    tD
  }

