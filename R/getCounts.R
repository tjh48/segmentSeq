getCounts <- function(segments, aD, preFiltered = FALSE, adjustMultireads = TRUE, useChunk = FALSE, cl = cl)
  {    
    if(class(aD) == "alignmentData") {
      counts <- matrix(NA, ncol = ncol(aD), nrow = length(segments))
      for(ss in levels(strand(segments)))
        if(any(strand(segments) == ss))
          counts[which(strand(segments) == ss),] <- .getCounts(segments[strand(segments) == ss,],
                               aD = aD[which(strand(aD@alignments) %in% list(c("+", "*"), c("-", "*"), c("+", "-", "*"))[[which(c("+", "-", "*") == ss)]]),],
                               preFiltered = preFiltered, as.matrix = TRUE, useChunk = useChunk, cl = list(NULL, cl)[[as.integer(length(segments) > 1) + 1]])
    } else if(class(aD) == "alignmentMeth") {
      Cs <- Ts <- matrix(NA, ncol = ncol(aD), nrow = length(segments))
      for (ss in levels(strand(segments))) {
        if(any(strand(segments) == ss)) {
          strandCounts <- .getMethylatedCounts(
                            segments = segments[strand(segments) == ss,],
                            mD = aD[which(strand(aD@alignments) %in% list("+", "-", c("+", "-", "*"))[[which(c("+", "-", "*") == ss)]]),],
                            preFiltered = preFiltered, as.matrix = TRUE, adjustMultireads = adjustMultireads, cl = list(NULL, cl)[[as.integer(length(segments) > 1) + 1]])
          Cs[which(strand(segments) == ss),] <- strandCounts$Cs
          Ts[which(strand(segments) == ss),] <- strandCounts$Ts
        }
      }
      counts <- list(Cs = Cs, Ts = Ts)
    }
    counts
  }
  
.getMethylatedCounts <- function(segments, mD, preFiltered = FALSE, as.matrix = FALSE, adjustMultireads = TRUE, cl)
  {     
    if(!is.null(cl))
      {
        clusterAssign <- function(assignList)
          {
            lapply(assignList, function(x) {
              assign(x$name, x$data, envir = .GlobalEnv)
              return(NULL)
            })
            return(NULL)
          }

        getCountEnv <- new.env(parent = .GlobalEnv)
        environment(clusterAssign) <- getCountEnv
      }

    alover <- getOverlaps(mD@alignments, segments, whichOverlaps = FALSE, cl = NULL)
    mD <- mD[alover,]
        
    alignments <- mD@alignments
    dataCs <- mD@Cs
    dataTs <- mD@Ts

    if(length(segments) == 0) {
      Tnas <- Cnas <- matrix(ncol = ncol(dataCs), nrow = 0)
      return(list(Cs = Cnas, Ts = Cnas))
#      if(as.matrix) return(list(Cs = Cnas, Ts = Cnas))
#      if(!as.matrix) {
#        Cdat <- do.call("DataFrame", lapply(1:ncol(Cnas), function(jj) Rle(Cnas[,jj])))
#        Tdat <- do.call("DataFrame", lapply(1:ncol(Tnas), function(jj) Rle(Tnas[,jj])))
#        return(list(Cs = Cdat, Ts = Tdat))
#      }
    }      

    if(!preFiltered)
      {
        segnas <- as.vector(is.na(seqnames(segments)) | is.na(start(segments)) | is.na(end(segments)))
        segments <- segments[!segnas,,drop = FALSE]

        rodering <- order(as.integer(seqnames(segments)), start(segments), end(segments))
        rodsegs <- segments[rodering,, drop = FALSE]
        dup <- which(!duplicated(rodsegs))
        reps <- c(dup[-1], length(rodsegs) + 1) - c(dup)
        redsegs <- rodsegs[dup,,drop = FALSE]
      } else redsegs <- segments

    countslist <- lapply(seqlevels(redsegs), function(cc)
                         {                                 
#                           createIntervals <- function(inCluster = FALSE)
#                             {
#                               cummaxEnd <- cummax(end(dupTags))
#                               cumminStart <- cummax(start(dupTags))
#                               
#                               fIns <- cbind(findInterval(start(chrsegs), cummaxEnd) + 1, 
#                                             findInterval(end(chrsegs), cumminStart))
#                               if(inCluster) assign("fIns", fIns, envir = .GlobalEnv) else return(fIns)
#                               return(NULL)
#                             }
#
                           chrsegs <- IRanges(start = start(redsegs[seqnames(redsegs) == cc,]), end = end(redsegs[seqnames(redsegs) == cc,]))
                           if(length(chrsegs) == 0)
                             return(list(Cs = matrix(ncol = ncol(dataCs), nrow = 0), Ts = matrix(ncol = ncol(dataTs), nrow = 0)))
                           
                           
                           whchr <- which(as.character(seqnames(alignments)) == cc & start(alignments) <= max(end(chrsegs)) & end(alignments) >= min(start(chrsegs)))
                           
                           chrCs <- dataCs[whchr,,drop = FALSE]
                           chrTs <- dataTs[whchr,,drop = FALSE]
                           chralignments <- alignments[whchr,]
                           if(is.null(chralignments$multireads)) multireads <- 1 else multireads <- as.integer(chralignments$multireads)
                           
                           adjustCs <- do.call("cbind", lapply(1:ncol(chrCs), function(rr) as.numeric(chrCs[,rr]) / multireads))
                           adjustTs <- do.call("cbind", lapply(1:ncol(chrTs), function(rr) as.integer(chrTs[,rr]) / multireads))

                           if(length(chralignments) > 0)
                             {
                               ordTags <- order(start(chralignments), end(chralignments))
                               droTags <- order(end(chralignments), start(chralignments))
                               endsAbove <- findInterval(end(chrsegs) + 0.5, start(chralignments)[ordTags])
                               startsAbove <- findInterval(start(chrsegs) - 0.5, end(chralignments)[droTags])
                               
                               cenCs <- rbind(0L, apply(adjustCs[droTags,,drop = FALSE], 2, cumsum))
                               cstCs <- rbind(0L, apply(adjustCs[ordTags,,drop = FALSE], 2, cumsum))
                               chrUCs <- cstCs[endsAbove + 1L,,drop = FALSE] - cenCs[startsAbove + 1L,,drop = FALSE]
                               chrUCs[chrUCs < 0] <- 0L
                               
                               cenTs <- rbind(0L, apply(adjustTs[droTags,,drop = FALSE], 2, cumsum))
                               cstTs <- rbind(0L, apply(adjustTs[ordTags,,drop = FALSE], 2, cumsum))
                               chrUTs <- cstTs[endsAbove + 1L,,drop = FALSE] - cenTs[startsAbove + 1L,,drop = FALSE]
                               chrUTs[chrUTs < 0] <- 0L
                             } else {
                               chrUCs <- matrix(0L, ncol = ncol(chrCs), nrow = length(chrsegs))
                               chrUTs <- matrix(0L, ncol = ncol(chrTs), nrow = length(chrsegs))
                             }
                           
                           
#                           if((adjustMultireads & !is.null(mD@duplicates)) {
#                             
#
#                             
#                             dupTags <- ranges(chralignments)[rowSums(dupstr != "") > 0,]
#                             dupInfo <- dupstr[rowSums(dupstr != "") > 0,]
#                             
#                             if(length(dupTags) > 0)
#                               {
#                                 fIns <- createIntervals()
#                                 
#                                 splitDups <- function(withinCluster, Cid = FALSE) {
#                                   if(!Cid)
#                                     {
#                                       splitdup <- lapply(1:ncol(dupInfo), function(ii) strsplit(dupInfo[,ii], ","))
#                                     } else splitdup <- lapply(1:ncol(dupInfo), function(ii) strsplit(gsub("[0-9]*_[0-9]*_", "", dupInfo[,ii]), ","))
#                                   if(withinCluster) {
#                                     assign(c("splitdup", "splitdupC")[as.numeric(Cid) + 1], splitdup, envir = .GlobalEnv)
#                                     return(invisible(NULL))
#                                   } else return(splitdup)
#                                 }
#                                 
#                                 countNonUniques <- function(segii)
#                                     {                                   
#                                       if(fIns[segii,1L] > fIns[segii,2L]) {
#                                         return(list(Cs = rep(0L, ncol(dupInfo)), Ts = rep(0L, ncol(dupInfo))))
#                                       }
#                                       seltags <- fIns[segii,1L]:fIns[segii,2L]
#                                       subtags <- dupTags[seltags,, drop = FALSE]
#                                       seltags <- seltags[start(subtags) <= end(chrsegs)[segii] & end(subtags) >= start(chrsegs)[segii]]#                                                                           
#                                   
#                                   dupcolcts <- function(ii)
#                                     {
#                                       if(all(dupInfo[seltags,ii] == ""))
#                                         return(c(C = 0, T = 0))
#                                       spldup <- unlist(splitdup[[ii]][seltags])                                       
#                                       if(length(spldup) == 0) {
#                                         return(c(C = 0, T = 0))
#                                       } else {
#                                         Ccts <- sum(unlist(splitdupC[[ii]][seltags])[duplicated(spldup)] == "C")
#                                         Tcts <- sum(unlist(splitdupC[[ii]][seltags])[duplicated(spldup)] == "T")
#                                         return(c(C = Ccts, T = Tcts))
#                                       }
#                                     }
#                                   
#                                   cts <- sapply(1:ncol(dupInfo), dupcolcts)
#                                   
#                                   return(list(Cs = cts[1,], Ts = cts[2,]))
#                                }              
#                              
#                                whchk <- which(fIns[,1] < fIns[,2])
#                                rofin <- order(fIns[whchk,1], fIns[whchk,2])
#                                unqin <- .fastUniques(fIns[whchk,][rofin,])
#                                 
#                                 if(sum(unqin) > 1) {
#                                   if(!is.null(cl) & sum(unqin) > length(cl)) {
#                                     
#                                     environment(createIntervals) <- getCountEnv
#                                     environment(splitDups) <- getCountEnv
#                                     clusterCall(cl, clusterAssign,
#                                                 assignList = list(list(name = "chrsegs", data = chrsegs),
#                                        #                                                 list(name = "splitdup", data = splitdup),
#                                        #                                                 list(name = "splitdupC", data = splitdupC),
#                                                   list(name = "dupTags", data = dupTags),
#                                                   list(name = "dupInfo", data = dupInfo)))
#                                     clusterCall(cl, splitDups, withinCluster = TRUE, Cid = FALSE)
#                                     clusterCall(cl, splitDups, withinCluster = TRUE, Cid = TRUE)
#                                     clusterCall(cl, createIntervals, inCluster = TRUE)
#                                     environment(countNonUniques) <- getCountEnv
#                                     
#                                     dupCeeTee <- parLapply(cl, whchk[rofin[unqin]], countNonUniques)
#                                   } else {
#                                     splitdup <- splitDups(withinCluster = FALSE, Cid = FALSE)
#                                     splitdupC <- splitDups(withinCluster = FALSE, Cid = TRUE)
#                                     dupCeeTee <- lapply(whchk[rofin[unqin]], countNonUniques)
#                                   }
#                                     
#                                   repunq <- rep(1:sum(unqin), c(which(unqin)[-1] - 1, length(unqin)) - which(unqin) + 1)               #                    
#                                   chrDCs[whchk[rofin],] <- do.call("rbind", lapply(dupCeeTee, function(x) x$Cs))[repunq,]
#                                   chrDTs[whchk[rofin],] <- do.call("rbind", lapply(dupCeeTee, function(x) x$Ts))[repunq,]
#                                 }                                      
#                               }
#                           }
                           list(Cs = chrUCs, Ts = chrUTs)
                         })
                                        #countsmat <- matrix(countsmat, nrow = nrow(redsegs), ncol = length(mD@replicates), byrow = TRUE)
    Cs <- do.call("rbind", lapply(countslist, function(x) x$Cs))
    Ts <- do.call("rbind", lapply(countslist, function(x) x$Ts))
    
    if(!preFiltered) {
      Cs <- matrix(unlist(lapply(1:length(reps), function(x) rep(Cs[x,,drop = FALSE], reps[x]))), nrow = length(rodsegs), ncol = length(mD@replicates), byrow = TRUE)
      Cs <- Cs[order(rodering),, drop = FALSE]
      
      Cnas <- matrix(NA, nrow = length(segnas), ncol = ncol(Cs))
      Cnas[!segnas,] <- Cs
      
      Ts <- matrix(unlist(lapply(1:length(reps), function(x) rep(Ts[x,,drop = FALSE], reps[x]))), nrow = length(rodsegs), ncol = length(mD@replicates), byrow = TRUE)
      Ts <- Ts[order(rodering),, drop = FALSE]
      
      Tnas <- matrix(NA, nrow = length(segnas), ncol = ncol(Ts))
      Tnas[!segnas,] <- Ts

    } else {
      Cnas <- Cs
      Tnas <- Ts
    }

    if(!is.null(cl))
      clusterEvalQ(cl, rm(list = ls()))      

#    if(!as.matrix) {
#      Cdat <- do.call("DataFrame", lapply(1:ncol(Cnas), function(jj) Rle(Cnas[,jj])))
#      Tdat <- do.call("DataFrame", lapply(1:ncol(Tnas), function(jj) Rle(Tnas[,jj])))
#      return(list(Cs = Cdat, Ts = Tdat))
#    }
    return(list(Cs = Cnas, Ts = Tnas))
  }


.getCounts <- function(segments, aD, preFiltered = FALSE, useChunk = FALSE, as.matrix = FALSE, cl)
  {
    if(!is.null(cl))
      {
        clusterAssign <- function(assignList)
          {
            lapply(assignList, function(x) {
              assign(x$name, x$data, envir = .GlobalEnv)
              return(NULL)
            })
            return(NULL)
          }

        getCountEnv <- new.env(parent = .GlobalEnv)
        environment(clusterAssign) <- getCountEnv
      }

    alignments <- aD@alignments
    cdata <- aD@data
    gc()
    
    if(!preFiltered)
      {
        segnas <- as.vector(is.na(seqnames(segments)) | is.na(start(segments)) | is.na(end(segments)))
        segments <- segments[!segnas,,drop = FALSE]

        rodering <- order(as.integer(seqnames(segments)), start(segments), end(segments))
        rodsegs <- segments[rodering,, drop = FALSE]
        dup <- which(!duplicated(rodsegs))
        reps <- c(dup[-1], length(rodsegs) + 1) - c(dup)
        redsegs <- rodsegs[dup,,drop = FALSE]
      } else redsegs <- segments
    
    countsmat <- do.call("rbind", lapply(seqlevels(redsegs), function(cc)
                               {
                                 createIntervals <- function(inCluster = FALSE)
                                   {
                                     cummaxEnd <- cummax(end(dupTags))
                                     cumminStart <- cummax(start(dupTags))
                                     
                                     fIns <- cbind(findInterval(start(chrsegs), cummaxEnd) + 1, 
                                                   findInterval(end(chrsegs), cumminStart))
                                     if(inCluster) assign("fIns", fIns, envir = .GlobalEnv) else return(fIns)
                                     return(NULL)
                                   }
                                 
                                 chrsegs <- IRanges(start = start(redsegs[seqnames(redsegs) == cc,]), end = end(redsegs[seqnames(redsegs) == cc,]))
                                 if(length(chrsegs) == 0) return(matrix(ncol = ncol(cdata), nrow = 0))
                                 
                                 whchr <- which(as.character(seqnames(alignments)) == cc & start(alignments) <= max(end(chrsegs)) & end(alignments) >= min(start(chrsegs)))
                                 
                                 intData <- cdata[whchr,,drop = FALSE]
                                 
                                 chralignments <- alignments[whchr,]

                                 if("tags" %in% names(values(chralignments))) {
                                   if("chunkDup" %in% names(values(chralignments)) & useChunk) {
                                     nondupTags <- ranges(chralignments)[!chralignments$chunkDup,]
                                     nondupData <- intData[!chralignments$chunkDup,, drop = FALSE]
                                   } else {
                                     nondupTags <- ranges(chralignments)[chralignments$multireads == 1,]
                                     nondupData <- intData[chralignments$multireads == 1,, drop = FALSE]
                                   }
                                 } else {
                                   nondupTags <- ranges(chralignments)
                                   nondupData <- intData
                                 }
                                 
                                 if(length(nondupTags) > 0)
                                   {
                                     ordTags <- order(start(nondupTags), end(nondupTags))
                                     droTags <- order(end(nondupTags), start(nondupTags))
                                     
                                     cens <- rbind(0L, apply(nondupData[droTags,,drop = FALSE], 2, cumsum))
                                     csts <- rbind(0L, apply(nondupData[ordTags,,drop = FALSE], 2, cumsum))
                                     
                                     endsAbove <- findInterval(end(chrsegs) + 0.5, start(nondupTags)[ordTags])
                                     startsAbove <- findInterval(start(chrsegs) - 0.5, end(nondupTags)[droTags])
                                     
                                     chrUC <- csts[endsAbove + 1L,,drop = FALSE] - cens[startsAbove + 1L,,drop = FALSE]
                                     chrUC[chrUC < 0] <- 0L
                                   } else chrUC <- matrix(0L, ncol = ncol(intData), nrow = length(chrsegs))

                                 chrNC <- matrix(0L, ncol = ncol(intData), nrow = length(chrsegs))
                                 gc()

                                 if("tag" %in% names(values(alignments))) {
                                   if("chunkDup" %in% names(values(chralignments)) & useChunk) {
                                     dup <- which(!chralignments$chunkDup)
                                   } else dup <- which(chralignments$multireads > 1)

                                   selNC <- dup[getOverlaps(chralignments[dup,],
                                                            GRanges(seqnames = cc, chrsegs), whichOverlaps = FALSE, cl = NULL)]
                                   dupTags <- ranges(chralignments)[selNC,]
                                   dupTagID <- values(chralignments)$tag[selNC]
                                   dupTagID <- match(dupTagID, unique(dupTagID))
                                   dupData <- intData[selNC,, drop = FALSE]

                                   if(length(dupTags) > 0)
                                     {
                                       countNonUniques <- function(segii)
                                         {
                                           if(fIns[segii,1L] > fIns[segii,2L]) return(rep(0L, ncol(dupData)))
                                           seltags <- fIns[segii,1L]:fIns[segii,2L]
                                           tags <- dupTags[seltags,, drop = FALSE]
                                           seltags <- seltags[start(tags) <= end(chrsegs)[segii] & end(tags) >= start(chrsegs)[segii]]
                                           seltags <- seltags[!duplicated(dupTagID[seltags])]
                                           as.integer(colSums(dupData[seltags,,drop = FALSE]))
                                         }

                                       if(is.null(cl) | length(chrsegs) == 1) {
                                         fIns <- createIntervals()
                                         chrNC <- do.call("rbind", lapply(1:length(chrsegs), countNonUniques))
                                       } else {

                                         environment(createIntervals) <- getCountEnv
                                         environment(countNonUniques) <- getCountEnv

                                         splitdup <- max(1, round(prod(dim(dupData)) / 1e6))

                                         for(pp in 1:splitdup) {                                           
                                           if(splitdup > 1) spseg <- which(cut(1:length(chrsegs), breaks = splitdup, labels = FALSE) == pp) else spseg <- 1:length(chrsegs)
                                           sptag <- which(start(dupTags) <= max(end(chrsegs[spseg])) & end(dupTags) >= min(start(chrsegs[spseg])))
                                           if(length(sptag) > 0)
                                             {
                                               clusterCall(cl, clusterAssign,
                                                           assignList = list(list(name = "chrsegs", data = chrsegs[spseg]),
                                                             list(name = "dupTags", data = dupTags[sptag]),
                                                             list(name = "dupTagID", data = dupTagID[sptag]),
                                                             list(name = "dupData", data = dupData[sptag,,drop = FALSE])))
                                               clusterCall(cl, createIntervals, inCluster = TRUE)
                                               chrNC[spseg,] <- do.call("rbind", parLapplyLB(cl, 1:length(spseg), countNonUniques))
                                             }
                                         }
                                       }
                                       gc()
                                     }
                                 }
                                 chrUC + chrNC
                               }))
    #countsmat <- matrix(countsmat, nrow = nrow(redsegs), ncol = length(aD@replicates), byrow = TRUE)
    
    if(!preFiltered) {
      countsmat <- matrix(unlist(lapply(1:length(reps), function(x) rep(countsmat[x,,drop = FALSE], reps[x]))), nrow = length(rodsegs), ncol = length(aD@replicates), byrow = TRUE)
      countsmat <- countsmat[order(rodering),, drop = FALSE]
      
      countsnas <- matrix(NA, nrow = length(segnas), ncol = ncol(countsmat))
      countsnas[!segnas,] <- countsmat
    } else countsnas <- countsmat

    #if(!as.matrix) countData <- do.call("DataFrame", lapply(1:ncol(countsnas), function(jj) Rle(countsnas[,jj]))) else countData <- countsnas      

    #countData
    countsnas
  }

