% modification on git from copied files
.findMethChunks <- function(alignments, gap)
{
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

  values(alignments)$chunk <- chunks

  alignments
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
    if(length(chrSS) > 0) {
        startstop <- reduce(chrSS@ranges, min.gapwidth = 0)
        ssChunk <- values(chrSS)$chunk[match(start(startstop), start(chrSS))]
        chunkDiff <- runLength(ssChunk)
      
        startRep <- (rep(cumsum(chunkDiff), chunkDiff) - 1:length(startstop) + 1)
        endRep <- cbind(1:length(startstop), rep(cumsum(chunkDiff) + 1, chunkDiff) - 1)
      
        csegs <- GRanges(seqnames = cc, IRanges(start = rep(start(startstop), startRep),
                         end = end(startstop)[unlist(lapply(1:nrow(endRep), function(ii) endRep[ii,1]:endRep[ii,2]))]), seqinfo = seqinfo(cTags))
      
      if(!missing(strand)) strand(csegs) <- strand
    
      csegs$csegChunk <- Rle(as.integer(rep(ssChunk, startRep)))
    } else csegs <- GRanges()
    if(verbose) message("done!")
    csegs
  }

.squeezeAlign <- function(sqAD, squeeze = 10, strand)
  {
    if(!missing(strand)) sqAD <- sqAD[strand(sqAD) == strand,]
    if(length(sqAD) > 1) {
      chunkRed <- findChunks(sqAD, gap = squeeze, checkDuplication = FALSE)      
      chunks <- values(chunkRed)$chunk
      chunks <- unique(chunks[!is.na(chunks)])

      chunkSeqs <- chunkRed[match(chunks, as.integer(values(chunkRed)$chunk)),]
      end(chunkSeqs) <- end(chunkRed)[c(match(chunks, as.integer(values(chunkRed)$chunk))[-1] - 1, length(chunkRed))]
      
      chunkSeqs <- do.call("c", lapply(seqlevels(sqAD), function(chr) {
        chrchunk <- chunkSeqs[seqnames(chunkSeqs) == chr,]
        chrmd <- sqAD[seqnames(sqAD) == chr,]
        matchmd <- match(start(chrchunk), start(chrmd))
        values(chrchunk)$chunk <- values(chrmd)$chunk[matchmd]
        chrchunk
      }))      
      if(missing(strand)) strand(chunkSeqs) = "*"
    } else chunkSeqs <- sqAD
    chunkSeqs$multireads <- NULL
    
    chunkSeqs
  }

.filterChunks <- function(aD, gap, filterProp)
  {
    if(class(aD) == "alignmentData") {
      filaD <- findChunks(aD@alignments, gap, checkDuplication = FALSE)      
    } else if(class(aD) == "alignmentMeth") {
      if(!missing(filterProp)) {
          filaD <- aD@alignments[which(rowSums(.methFunction(aD, prop = filterProp, locCutoff = NA), na.rm = TRUE) > 0),]        
      } else filaD <- aD[rowSums(aD@Cs) > 0,]
      filaD <- .findMethChunks(filaD, gap)
    }
    filaD <- filaD[!duplicated(filaD),]
    filaD
  }


processAD <- function(aD, gap = 300, squeeze = 2, filterProp = 0.05, strandSplit = FALSE, verbose = TRUE, getCounts = FALSE, cl)         
  {
#    if("tag" %in% colnames(values(aD@alignments))) {      

                                        #      values(aD@alignments)$tag <- as.integer(as.factor(values(aD@alignments)$tag))
#    } else values(aD@alignments)$tag <- 1:nrow(aD)

      filaD <- .filterChunks(aD, gap, filterProp)
    
    if(squeeze > 0) {
      if(strandSplit) {
        cTags <- c(.squeezeAlign(filaD, squeeze = squeeze, strand = "+"),                  
                  .squeezeAlign(filaD, squeeze = squeeze, strand = "-"))
      } else cTags <- .squeezeAlign(filaD, squeeze = squeeze)
    } else cTags <- filaD

    coordinates <- GRanges()
    seqinfo(coordinates) <- seqinfo(aD@alignments)
    
    data <- Cs <- Ts <- NULL
    
    partCounts <- function(chunks, cs, cc) {
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
    
    chrDat <- lapply(seqlevels(cTags), function(cc) {
      if(any(cTags@seqnames == cc))
        {
          if(strandSplit) {
            cs <- c(.chrProcessing(cTags, cc, strand = "+", verbose = verbose),
                    .chrProcessing(cTags, cc, strand = "-", verbose = verbose))              
          } else {
            cs <- .chrProcessing(cTags, cc, verbose = verbose)
          }
          cs <- cs[order(as.integer(cs$csegChunk), start(cs), end(cs)),]
          
          if(verbose) message(length(cs), " candidate loci found.")
          
          if(getCounts) {
            if(verbose) message("Getting count data for each candidate locus...", appendLF = FALSE)
            windowChunks <- .windowing(cs)
            windata <- lapply(windowChunks, partCounts, cs = cs, cc = cc)
            
            if(class(aD) == "alignmentData") {
              data <- do.call("rbind", windata)
            } else if(class(aD) == "alignmentMeth") {
              Cs <- do.call("rbind", lapply(windata, function(x) x$Cs))
              Ts <- do.call("rbind", lapply(windata, function(x) x$Ts))
              rm(windata)
              gc()
            }

            if(verbose)
              message("...done!")            
            
          } else data <- Cs <- Ts <- matrix(nrow = 0, ncol = ncol(aD))
                    
          gc()
          
        } else {
          if(verbose) message("No tags found for chromosome ", cc)
          cs <- GRanges()
          data <- Cs <- Ts <- matrix(nrow = 0, ncol = ncol(aD))
        }
      if(class(aD) == "alignmentData") {
        return(list(data = data, coordinates = cs))
      } else if(class(aD) == "alignmentMeth") return(list(Cs = Cs, Ts = Ts, coordinates = cs))                                              
    })

    coordinates <- do.call("c", lapply(chrDat, function(x) x$coordinates))
    values(coordinates) <- NULL

    if(class(aD) == "alignmentData") {
      data <- do.call("rbind", lapply(chrDat, function(x) x$data))
      tD <- new("lociData", coordinates = coordinates, data = data, libsizes = aD@libsizes, replicates = aD@replicates)
      colnames(tD@data) <- aD@libnames
    } else if(class(aD) == "alignmentMeth") {
      Cs <- do.call("rbind", lapply(chrDat, function(x) x$Cs))
      Ts <- do.call("rbind", lapply(chrDat, function(x) x$Ts))
      colnames(Cs) <- colnames(Ts) <- aD@libnames
      tD <- new("lociData")
      tD@sampleObservables$nonconversion <- aD@nonconversion
      tD@replicates = aD@replicates
      tD@coordinates = coordinates
      tD@data <- array(c(Cs, Ts), c(dim(Cs), 2))
#      tD@Cs = Cs
#      tD@Ts = Ts
  }
      tD@coordinates$chunk <- findChunks(tD@coordinates, gap = gap, checkDuplication = FALSE, justChunks = TRUE)
    tD
  }
