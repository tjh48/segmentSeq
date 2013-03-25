.windowing <- function(cs, winsize = 1e5)
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

.squeezeAlign <- function(aD, squeeze = 10)
  {     
    if(class(aD) == "alignmentData") {
      chunkRed <- findChunks(aD@alignments, gap = squeeze, checkDuplication = FALSE)
    } else if(class(aD) == "methAlignment") {
      aD <- aD[rowSums(sapply(1:ncol(aD), function(ii) as.integer(aD@Cs[,ii]))) > 0,]
      chunkRed <- findMethChunks(aD, gap = squeeze)
    }
    
    chunks <- values(chunkRed@alignments)$chunk
    chunks <- unique(chunks[!is.na(chunks)])
    
    chunkSeqs <- chunkRed@alignments[match(chunks, as.integer(values(chunkRed@alignments)$chunk)),]
    end(chunkSeqs) <- end(chunkRed@alignments)[c(match(chunks, as.integer(values(chunkRed@alignments)$chunk))[-1] - 1, nrow(chunkRed))]
    strand(chunkSeqs) = "*"
    chunkSeqs

    chunkSeqs <- do.call("c", lapply(seqlevels(aD@alignments), function(chr) {
      chrchunk <- chunkSeqs[seqnames(chunkSeqs) == chr,]
      chrmd <- aD@alignments[seqnames(aD@alignments) == chr,]
      matchmd <- match(start(chrchunk), start(chrmd))
      values(chrchunk)$chunk <- values(chrmd)$chunk[matchmd]
      chrchunk
    }))

    chunkSeqs
  }



processAD <-
function(aD, gap = NULL, squeeze = 0, # filterProp,
         verbose = TRUE, cl)
  {
#    if("tag" %in% colnames(values(aD@alignments))) {      

                                        #      values(aD@alignments)$tag <- as.integer(as.factor(values(aD@alignments)$tag))
#    } else values(aD@alignments)$tag <- 1:nrow(aD)
    
    if(!is.null(gap)) {
      if(class(aD) == "alignmentData") {
        aD@alignments <- findChunks(aD@alignments, gap)
      } #else if(class(aD) == "methAlignment") {
#        aD <- findMethChunks(aD, gap)    
#      }
    }

#    if(!missing(filterProp)) {
#      filaD <- aD[rowSums(.methFunction(aD, prop = filterProp, locCutoff = 0.5), na.rm = TRUE) > 0,]
#      filaD <- findMethChunks(filaD, gap)
#    } else  {
    filaD <- aD
                                        #  }
    
    if(squeeze > 0) {
      sqaD <- .squeezeAlign(filaD, squeeze = squeeze)
      cTags <- sqaD
    } else cTags <- filaD@alignments
    
    coordinates <- GRanges()
    seqinfo(coordinates) <- seqinfo(aD@alignments)

    data <- Cs <- Ts <- DataFrame()
    
    partCounts <- function(chunks) {
      chad <- which(as.integer(values(cTags)$chunk) %in% chunks)
      
      waD <- aD[(seqnames(aD@alignments) == cc & end(aD@alignments) >= min(start(cTags[chad])) & start(aD@alignments) <= max(end(cTags[chad]))),]
      
                                        #              waD@alignments <- cTags[chad,,drop = FALSE]
                                        #              waD@data <- data[chad,,drop = FALSE]
                                        #if(chunks %% 100 == 0)
      message(".", appendLF = FALSE)
      
      segments <- GRanges(seqnames = cc, cs$csegs[which(as.integer(cs$csegChunk) %in% chunks),])
      
      if(class(waD) == "alignmentData") {
        counts <- getCounts(segments = segments, 
                            aD = waD, preFiltered = FALSE, cl = list(NULL, cl)[[as.integer(length(segments) > 1) + 1]])
      } else if(class(aD) == "methAlignment") {
        counts <- getMethylatedCounts(segments = segments, 
                                      mD = waD, preFiltered = FALSE,
                                      cl = cl)
      }
      return(counts)
    }
    
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
            windata <- lapply(windowChunks, partCounts)

#            save(windata, file = paste("temp/", cc, "_count.RData", sep = ""))
            
            if(class(aD) == "alignmentData") {
              data <- rbind(data, do.call("rbind", windata))                                          
            } else if(class(aD) == "methAlignment") {
              Cs <- rbind(Cs, do.call("rbind", lapply(windata, function(x) x$Cs)))
#              for(ww in 1:length(windata)) windata[[ww]]$Cs <- NULL
#              gc()
              Ts <- rbind(Ts, do.call("rbind", lapply(windata, function(x) x$Ts)))
              rm(windata)
              gc()
            }
            
            if(verbose)
              message("...done!")

            gcoord <- GRanges(seqnames = cc, cs$csegs, seqlengths = seqlengths(coordinates))
            seqinfo(gcoord, force = TRUE) <- seqinfo(coordinates)
            coordinates <- suppressWarnings(c(coordinates, gcoord))
            
            gc()

          } else if(verbose) message("No tags found for this chromosome.")
      }

    if(class(aD) == "alignmentData") {
      tD <- new("segData", coordinates = coordinates, data = data, libsizes = aD@libsizes, replicates = aD@replicates)
      colnames(tD@data) <- aD@libnames
    } else if(class(aD) == "methAlignment") {
      tD <- new("methSegs", coordinates = coordinates, Cs = Cs, Ts = Ts, replicates = aD@replicates)
      colnames(tD@Cs) <- colnames(tD@Ts) <- aD@libnames
    }
    tD
  }
