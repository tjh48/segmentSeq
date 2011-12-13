getCounts <-
function(segments, aD, preFiltered = FALSE, as.matrix = FALSE, cl)
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
    if("tag" %in% names(values(alignments)) && class(values(alignments)$tag) != "integer") values(alignments)$tag <- as.integer(as.factor(values(alignments)$tag)) else if(!("tag" %in% names(values(alignments)))) values(alignments)$tag <- 1:length(alignments)
    cdata <- aD@data
    
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
                                 
                                 chrdata <- cdata[whchr,]
                                 intData <- sapply(1:ncol(chrdata), function(rr) as.integer(chrdata[,rr]))
                                 
                                 chralignments <- alignments[whchr,]

                                 nondupTags <- ranges(chralignments)[!values(chralignments)$chunkDup,]
                                 nondupData <- intData[!values(chralignments)$chunkDup,, drop = FALSE]
                                 
                                 if(length(nondupTags) > 0)
                                   {
                                     ordTags <- order(start(nondupTags), end(nondupTags))
                                     droTags <- order(end(nondupTags), start(nondupTags))
                                     
#                                     endsBelow <- findInterval(end(chrsegs), end(nondupTags)[droTags])
#                                     startsBelow <- findInterval(start(chrsegs) - 0.5, start(nondupTags)[ordTags])                                     
#                                     chrUC <- (cens[endsBelow + 1L,] - csts[startsBelow + 1L,])

                                     cens <- rbind(0L, apply(nondupData[droTags,,drop = FALSE], 2, cumsum))
                                     csts <- rbind(0L, apply(nondupData[ordTags,,drop = FALSE], 2, cumsum))
                                     
                                     endsAbove <- findInterval(end(chrsegs) + 0.5, start(nondupTags)[ordTags])
                                     startsAbove <- findInterval(start(chrsegs), end(nondupTags)[droTags])
                                     
                                     chrUC <- csts[endsAbove + 1L,] - cens[startsAbove + 1L,]
                                     chrUC[chrUC < 0] <- 0L
                                   } else chrUC <- matrix(0L, ncol = ncol(intData), nrow = nrow(chrsegs))
                                 
                                 dupTags <- ranges(chralignments)[values(chralignments)$chunkDup,]
                                 dupTagID <- values(chralignments)$tag[values(chralignments)$chunkDup]
                                 dupData <- intData[values(chralignments)$chunkDup,, drop = FALSE]
                                 
                                 if(length(dupTags) > 0)
                                   {
                                     if(!is.null(cl))
                                       {
                                         environment(createIntervals) <- getCountEnv
                                         clusterCall(cl, clusterAssign,
                                                     assignList = list(list(name = "chrsegs", data = chrsegs),
                                                       list(name = "dupTags", data = dupTags),
                                                       list(name = "dupTagID", data = dupTagID),
                                                       list(name = "dupData", data = dupData)))
                                         clusterCall(cl, createIntervals, inCluster = TRUE)
                                       } else fIns <- createIntervals()
                                     
                                     countNonUniques <- function(segii)
                                       {
                                         if(fIns[segii,1L] > fIns[segii,2L]) return(rep(0L, ncol(dupData)))
                                         seltags <- fIns[segii,1L]:fIns[segii,2L]
                                         tags <- dupTags[seltags,, drop = FALSE]
                                         seltags <- seltags[start(tags) <= end(chrsegs)[segii] & end(tags) >= start(chrsegs)[segii]]
                                         seltags <- seltags[!duplicated(dupTagID[seltags])]
                                         as.integer(colSums(dupData[seltags,,drop = FALSE]))
                                       }
                                     
                                     if(!is.null(cl))
                                       environment(countNonUniques) <- getCountEnv
                                     
                                     if(!is.null(cl)) chrNC <- parSapply(cl, 1:length(chrsegs), countNonUniques) else  chrNC <- sapply(1:length(chrsegs), countNonUniques)
                                                               
                                     
                                   } else chrNC <- matrix(0L, nrow = ncol(intData), ncol = length(chrsegs))
                                 
                                 chrUC + t(chrNC)
                               }))
    #countsmat <- matrix(countsmat, nrow = nrow(redsegs), ncol = length(aD@replicates), byrow = TRUE)
    
    if(!preFiltered) {
      countsmat <- matrix(unlist(lapply(1:length(reps), function(x) rep(countsmat[x,], reps[x]))), nrow = length(rodsegs), ncol = length(aD@replicates), byrow = TRUE)
      countsmat <- countsmat[order(rodering),, drop = FALSE]
      
      countsnas <- matrix(NA, nrow = length(segnas), ncol = ncol(countsmat))
      countsnas[!segnas,] <- countsmat
    } else countsnas <- countsmat

    if(!as.matrix) countData <- do.call("DataFrame", lapply(1:ncol(countsnas), function(jj) Rle(countsnas[,jj]))) else countData <- countsnas      

    countData
  }

