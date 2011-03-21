processAD <-
function(aD, gap = NULL, verbose = TRUE, cl)
  {

    if("tag" %in% colnames(aD@alignments)) aD@alignments$tag <- as.integer(as.factor(aD@alignments$tag)) else aD@alignments$tag <- 1:nrow(aD@alignments)
    
    if(!is.null(gap))
      aD <- findChunks(aD, gap)

    cTags <- aD@alignments    
    chrs <- aD@chrs$chr
    chrlens <- aD@chrs$len
    libnames <- aD@libnames
    libsizes <- aD@libsizes
    replicates <- aD@replicates
    tagData <- aD@data
    
    tD <- new("segData", segInfo = data.frame(), data = matrix(nrow = 0, ncol = length(replicates)), libsizes = libsizes, replicates = replicates, chrs = aD@chrs)

    for(cc in 1:length(chrs))
      {
        if(any(cTags$chr == chrs[cc]))
          {            
            chrSS <- subset(cTags, subset = cTags$chr == chrs[cc], select = c(start, end, chunk, chunkDup))

            if(any(chrSS$end > chrlens[cc]))
              warning(paste("Chromosome", chrs[cc], "has tags which extend over the given chromsome length."))

            if(verbose){
              message("Chromosome: ", chrs[cc])
              message("Finding start-stop co-ordinates...", appendLF = FALSE)
            }
            
            if(nrow(chrSS) > 1)
              {
                chrSS <- chrSS[order(chrSS$chunk, chrSS$start, chrSS$end),,drop = FALSE]
                chrmax <- cbind(cummax(chrSS$start), cummax(chrSS$end))
                ch <- which(chrmax[-1L,1L] > chrmax[-nrow(chrmax),2L])
                startstop <- cbind(starts = as.integer(c(min(chrmax[,1L]), chrmax[ch + 1L,1L])),
                                   ends = as.integer(c(chrmax[ch,2L], max(chrmax[,2L]))),
                                   chunk = chrSS$chunk[c(1, ch + 1L)])
              } else if(nrow(chrSS) == 1) startstop <- cbind(starts = chrSS$start, ends = chrSS$end, chunk = chrSS$chunk) else startstop <- matrix(nrow = 0, ncol = 0)

            chunkDups <- which(!duplicated(startstop[,"chunk"]))
            chunkDiff <- diff(c(chunkDups, nrow(startstop) + 1))
            startRep <- rep(c(chunkDups[-1] - 1, nrow(startstop)), chunkDiff) - 1:nrow(startstop) + 1
            endRep <- cbind(1:nrow(startstop), rep(c(chunkDups[-1], nrow(startstop) + 1), chunkDiff) - 1)
            csegs <- data.frame(start = rep(startstop[,1L], startRep),
                                end = startstop[unlist(lapply(1:nrow(endRep), function(ii) endRep[ii,1]:endRep[ii,2])),2L],
                                chunk = rep(startstop[,"chunk"], startRep))

            if(verbose){
              message(".", nrow(csegs), " found.")
              message("Getting count data for each potential subsegment...", appendLF = FALSE)
            }
            
            winsize <- 5e5

            chunkWindows <- which(!duplicated(csegs$chunk))
            chunkWindows <- cbind(csegs$chunk[chunkWindows], cumsum(diff(c(which(!duplicated(csegs$chunk)), nrow(csegs) + 1))))

            dupWin <- which(!duplicated(round(chunkWindows[,2] / winsize)))
            dupWin <- cbind(dupWin, c(dupWin[-1] - 1, nrow(chunkWindows)))
            windowChunks <- lapply(1:nrow(dupWin), function(ii) chunkWindows[dupWin[ii,1]:dupWin[ii,2],1L])

            waD <- new("alignmentData")
            waD@libnames <- libnames
            waD@libsizes <- libsizes
            waD@chrs <- data.frame(chr = chrs[cc], len = chrlens[cc])
            waD@replicates = replicates

            gc()
            
            tD@data <- rbind(tD@data, do.call("rbind",
                                              lapply(windowChunks, function(chunks) {
                                                chad <- which(aD@alignments$chunk %in% chunks)
                                                waD@alignments <- aD@alignments[chad,]
                                                waD@data <- aD@data[chad,]
                                                message(".", appendLF = FALSE)
                                                 getCounts(
                                                          segments = data.frame(chr = chrs[cc], subset(csegs, csegs$chunk %in% chunks, select = c(start, end)))
                                                          , aD = waD, preFiltered = TRUE, cl = cl)
                                              })))
                                                            
            if(verbose)
              message("...done!")

            csegs <- cbind(csegs,
                           leftSpace = (csegs[,1L] - c(0L, startstop[,2])[findInterval(csegs[,1L], startstop[,2L]) + 1L]) - 1L,
                           rightSpace = c(startstop[,1L], chrlens[cc] + 1L)[findInterval(csegs[,2L], startstop[,1L]) + 1L] - csegs[,2L] - 1L)
            
            tD@segInfo <- rbind(tD@segInfo, data.frame(chr = I(chrs[cc]), csegs))

            gc()

          } else if(verbose) message("No tags found for this chromosome.")
      }
      
    colnames(tD@data) <- aD@libnames
    
    rownames(tD@segInfo) <- rownames(tD@data) <- 1:nrow(tD@data)

    tD
  }

