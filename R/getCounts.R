
getCounts <-
function(segments, aD, cl)
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
    if("tag" %in% names(alignments)) alignments$tag <- as.integer(as.factor(alignments$tag)) else alignments$tag <- 1:nrow(alignments)
    cdata <- aD@data
    
    segnas <- is.na(segments$chr) | is.na(segments$start) | is.na(segments$end)
    segments <- segments[!segnas,]

    rodering <- order(segments$chr, segments$start, segments$end)
    rodsegs <- segments[rodering,, drop = FALSE]
    dup <- which(!duplicated(rodsegs))
    reps <- c(dup[-1], nrow(rodsegs) + 1) - c(dup)
    redsegs <- rodsegs[dup,]
    

    countsmat <- unlist(lapply(unique(redsegs$chr), function(cc)
                               {
                                 createIntervals <- function(inCluster = FALSE)
                                   {
                                     cummaxEnd <- cummax(dupTags$end)
                                     cumminStart <- cummax(dupTags$start)
                                     
                                     fIns <- cbind(findInterval(chrsegs$start, cummaxEnd) + 1, 
                                                   findInterval(chrsegs$end, cumminStart))
                                     if(inCluster) assign("fIns", fIns, envir = .GlobalEnv) else return(fIns)
                                     return(NULL)
                                   }

                                 
                                 chrsegs <- data.frame(start = redsegs$start, end = redsegs$end)[redsegs$chr == cc,,drop = FALSE]
                                 chralignments <- subset(alignments, subset = alignments$chr == cc, select = c(start, end, tag))
                                 chrdata <- subset(cdata, subset = alignments$chr == cc)
                                 
                                 chralignments <- cbind(chralignments, duplicated = chralignments$tag %in% chralignments$tag[duplicated(chralignments$tag)])
                                 
                                 nondupTags <- subset(chralignments, subset = chralignments$duplicated == FALSE, select = c(start, end))
                                 nondupData <- chrdata[chralignments$duplicated == FALSE,,drop = FALSE]
                                 
                                 if(nrow(nondupTags) > 0)
                                   {
                                     ordTags <- order(nondupTags$start, nondupTags$end)
                                     droTags <- order(nondupTags$end, nondupTags$start)
                                     
                                     cens <- rbind(0L, apply(nondupData[droTags,,drop = FALSE], 2, cumsum))
                                     
                                     endsBelow <- findInterval(chrsegs$end, nondupTags$end[droTags])
                                     csts <- rbind(0L, apply(nondupData[ordTags,,drop = FALSE], 2, cumsum))
                                     startsBelow <- findInterval(chrsegs[,1] - 0.5, nondupTags$start[ordTags])
                                     
                                     chrUC <- (cens[endsBelow + 1L,] - csts[startsBelow + 1L,])
                                     chrUC[chrUC < 0] <- 0
                                   } else chrUC <- matrix(0L, ncol = ncol(chrdata), nrow = nrow(chrsegs))
                                 
                                 dupTags <- subset(chralignments, subset = chralignments$duplicated == TRUE, select = c(start, end, tag))
                                 dupData <- chrdata[chralignments$duplicated == TRUE,, drop = FALSE]
                                 
                                 if(nrow(dupTags) > 0)
                                   {
                                     if(!is.null(cl))
                                       {
                                         environment(createIntervals) <- getCountEnv
                                         clusterCall(cl, clusterAssign,
                                                     assignList = list(list(name = "chrsegs", data = chrsegs),
                                                       list(name = "dupTags", data = dupTags),
                                                       list(name = "dupData", data = dupData)))
                                         clusterCall(cl, createIntervals, inCluster = TRUE)
                                       } else fIns <- createIntervals()
                                     
                                     countNonUniques <- function(segii)
                                       {
                                         if(fIns[segii,1L] > fIns[segii,2L]) return(rep(0, ncol(dupData)))
                                         seltags <- fIns[segii,1L]:fIns[segii,2L]
                                         tags <- dupTags[seltags,, drop = FALSE]
                                         seltags <- seltags[tags$start <= chrsegs$end[segii] & tags$end >= chrsegs$start[segii]]
                                         seltags <- seltags[!duplicated(dupTags$tag[seltags])]
                                         as.integer(colSums(dupData[seltags,,drop = FALSE]))
                                       }
                                     
                                     if(!is.null(cl))
                                       environment(countNonUniques) <- getCountEnv
                                     
                                     if(!is.null(cl)) chrNC <- parSapply(cl, 1:nrow(chrsegs), countNonUniques) else  chrNC <- sapply(1:nrow(chrsegs), countNonUniques)
                                   } else chrNC <- matrix(0L, nrow = ncol(chrdata), ncol = nrow(chrsegs))
                                 
                                 t(chrUC + t(chrNC))
                               }))
    countsmat <- matrix(countsmat, nrow = nrow(redsegs), ncol = length(aD@replicates), byrow = TRUE)
    countsmat <- matrix(unlist(sapply(1:length(reps), function(x) rep(countsmat[x,], reps[x]))), nrow = nrow(rodsegs), ncol = length(aD@replicates), byrow = TRUE)
    countsmat[order(rodering),, drop = FALSE]

    countsnas <- matrix(NA, nrow = length(segnas), ncol = ncol(countsmat))
    countsnas[!segnas,] <- countsmat

    countsnas
  }

