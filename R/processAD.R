.findMethChunks <- function(mD, gap)
{
  zeroCs <- rowSums(sapply(1:ncol(mD), function(ii) as.integer(mD@Cs[,ii]))) == 0  
  alignments <- mD@alignments[!zeroCs,]
  chunks <- Rle(rep(NA, length(alignments)))
  chunkID <- maxChunk <- 0L
  
  for(chr in seqlevels(alignments)) {
    chraD <- which(seqnames(alignments) == chr)
    if(length(chraD) > 0) {
      if(length(chraD) == 1) {
        chunks[chraD] <- maxChunk + 1
      } else if(length(chraD) > 1)
        {
          chral <- ranges(alignments[chraD,])
          chunkNum <- c(0L, which(start(chral)[-1] - cummax(end(chral)[-length(chral)]) > gap))
          chunkID <- rep(1L:as.integer(length(chunkNum)) + maxChunk, diff(c(chunkNum, length(chral))))
          chunks[chraD] <- chunkID
        }
      maxChunk <- max(chunks, na.rm = TRUE)
    }
  }

  allchunk <- Rle(NA, nrow(mD))
  allchunk[!zeroCs] <- chunks
  values(mD@alignments)$chunk <- allchunk
  
  mD
}


.windowing <- function(cs, winsize = 1e5)
  {
    chunks <- sort(cs$csegChunk)
    chunkWindows <- cbind(1:length(runLength(chunks)), cumsum(runLength(chunks)))
    
    dupWin <- which(!duplicated(round(chunkWindows[,2] / winsize)))
    dupWin <- cbind(dupWin, c(dupWin[-1] - 1, nrow(chunkWindows)))
    windowChunks <- lapply(1:nrow(dupWin), function(ii) runValue(chunks)[chunkWindows[dupWin[ii,1]:dupWin[ii,2],1L]])
    windowChunks
  }

.chrProcessing <- function(cTags, cc, strand, verbose)
  {
    if(!missing(strand)) {
      chrSS <- cTags[cTags@seqnames == cc & strand(cTags) == strand,]
      if(verbose) message(paste("Strand:", strand))
    } else chrSS <- cTags[cTags@seqnames == cc,]
    
    if(!is.na(seqlengths(chrSS)[seqlevels(chrSS) == cc]) && any(end(chrSS) > seqlengths(chrSS)[seqlevels(chrSS) == cc]))
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
    
    csegs <- GRanges(seqnames = cc, IRanges(start = rep(start(startstop), startRep),
                             end = end(startstop)[unlist(lapply(1:nrow(endRep), function(ii) endRep[ii,1]:endRep[ii,2]))]), seqinfo = seqinfo(cTags))

    if(!missing(strand)) strand(csegs) <- strand
    
    csegs$csegChunk <- Rle(as.integer(rep(ssChunk, startRep)))
    if(verbose) message("done!")
    csegs
  }

.squeezeAlign <- function(aD, squeeze = 10, strand)
  {
    if(!missing(strand)) aD <- aD[strand(aD@alignments) == strand,]
    
    if(class(aD) == "alignmentData") {
      chunkRed <- findChunks(aD@alignments, gap = squeeze, checkDuplication = FALSE)
    } else if(class(aD) == "alignmentMeth") {
      aD <- aD[rowSums(sapply(1:ncol(aD), function(ii) as.integer(aD@Cs[,ii]))) > 0,]
      chunkRed <- .findMethChunks(aD, gap = squeeze)
    }
    
    chunks <- values(chunkRed@alignments)$chunk
    chunks <- unique(chunks[!is.na(chunks)])
    
    chunkSeqs <- chunkRed@alignments[match(chunks, as.integer(values(chunkRed@alignments)$chunk)),]
    end(chunkSeqs) <- end(chunkRed@alignments)[c(match(chunks, as.integer(values(chunkRed@alignments)$chunk))[-1] - 1, nrow(chunkRed))]
    chunkSeqs

    chunkSeqs <- do.call("c", lapply(seqlevels(aD@alignments), function(chr) {
      chrchunk <- chunkSeqs[seqnames(chunkSeqs) == chr,]
      chrmd <- aD@alignments[seqnames(aD@alignments) == chr,]
      matchmd <- match(start(chrchunk), start(chrmd))
      values(chrchunk)$chunk <- values(chrmd)$chunk[matchmd]
      chrchunk
    }))

    if(missing(strand)) strand(chunkSeqs) = "*"

    chunkSeqs
  }



processAD <-
function(aD, gap = NULL, squeeze = 0, filterProp = 0.1,
         strandSplit = FALSE, verbose = TRUE, cl)
  {
#    if("tag" %in% colnames(values(aD@alignments))) {      

                                        #      values(aD@alignments)$tag <- as.integer(as.factor(values(aD@alignments)$tag))
#    } else values(aD@alignments)$tag <- 1:nrow(aD)

    if(!is.null(gap)) {
      if(class(aD) == "alignmentData") {
        aD@alignments <- findChunks(aD@alignments, gap)
      } else if(class(aD) == "alignmentMeth") {
        aD <- .findMethChunks(aD, gap)    
      }
    }

    if(!missing(filterProp) & class(aD) == "alignmentMeth") {
        filaD <- aD[which(rowSums(.methFunction(aD, prop = filterProp, locCutoff = NA), na.rm = TRUE) > 0),]
        filaD <- .findMethChunks(filaD, gap)
    } else  {
      filaD <- aD
    }
    
    if(squeeze > 0) {
      if(strandSplit) {
        cTags <- c(.squeezeAlign(filaD, squeeze = squeeze, strand = "+"),                  
                  .squeezeAlign(filaD, squeeze = squeeze, strand = "-"))
      } else cTags <- .squeezeAlign(filaD, squeeze = squeeze)
    } else cTags <- filaD@alignments
    
    coordinates <- GRanges()
    seqinfo(coordinates) <- seqinfo(aD@alignments)

    data <- Cs <- Ts <- NULL
    
    partCounts <- function(chunks) {
      segments <- cs[which(as.integer(cs$csegChunk) %in% chunks),]             
      chad <- which(as.integer(values(cTags)$chunk) %in% chunks)
      
      waD <- aD[which(seqnames(aD@alignments) == cc & end(aD@alignments) >= min(start(cTags[chad])) & start(aD@alignments) <= max(end(cTags[chad]))),]
      
                                        #              waD@alignments <- cTags[chad,,drop = FALSE]
                                        #              waD@data <- data[chad,,drop = FALSE]
                                        #if(chunks %% 100 == 0)
      message(".", appendLF = FALSE)
      
      counts <- getCounts(segments = segments, 
                          aD = waD, preFiltered = FALSE, useChunk = TRUE, cl = list(NULL, cl)[[as.integer(length(segments) > 1) + 1]])
      return(counts)
    }
    
    for(cc in seqlevels(cTags))
      {
        if(any(cTags@seqnames == cc))
          {
            if(strandSplit) {
              cs <- c(.chrProcessing(cTags, cc, strand = "+", verbose = verbose),
                      .chrProcessing(cTags, cc, strand = "-", verbose = verbose))              
            } else {
              cs <- .chrProcessing(cTags, cc, verbose = verbose)
            }
            cs <- cs[order(as.integer(cs$csegChunk), start(cs), end(cs)),]
            
#            cs <- .chrProcessing(cTags, cc, verbose = verbose)
#            if(strandSplit) {
#              cs <- c(cs, cs)
#              strand(cs) <- c(rep(c("+", "-"), each = length(cs) / 2))
#              cs <- cs[c(which(strand(cs) == "+")[getOverlaps(cs[strand(cs) == "+",], cTags[strand(cTags) == "+",], whichOverlaps = FALSE)],
#                         which(strand(cs) == "-")[getOverlaps(cs[strand(cs) == "-",], cTags[strand(cTags) == "-",], whichOverlaps = FALSE)]),]
#            }
            
            if(verbose){
              message(length(cs), " candidate loci found.")
              message("Getting count data for each candidate locus...", appendLF = FALSE)
            }
            
            windowChunks <- .windowing(cs)
            windata <- lapply(windowChunks, partCounts)

            if(class(aD) == "alignmentData") {
              data <- rbind(data, do.call("rbind", windata))                                          
            } else if(class(aD) == "alignmentMeth") {
              Cs <- rbind(Cs, do.call("rbind", lapply(windata, function(x) x$Cs)))
              Ts <- rbind(Ts, do.call("rbind", lapply(windata, function(x) x$Ts)))
              rm(windata)
              gc()
            }
            
            if(verbose)
              message("...done!")

            coordinates <- suppressWarnings(c(coordinates, cs))            
            gc()

          } else if(verbose) message("No tags found for chromosome ", cc)
      }

    values(coordinates) <- NULL

    if(class(aD) == "alignmentData") {
      tD <- new("segData", coordinates = coordinates, data = data, libsizes = aD@libsizes, replicates = aD@replicates)
      colnames(tD@data) <- aD@libnames
    } else if(class(aD) == "alignmentMeth") {
      colnames(Cs) <- colnames(Ts) <- aD@libnames
      tD <- new("segMeth")
      tD@nonconversion <- aD@nonconversion
      tD@replicates = aD@replicates
      tD@coordinates = coordinates
      tD@Cs = Cs
      tD@Ts = Ts
    }
    tD
  }
