clustSeg <- function(sD, aD, subRegion = NULL, getLikes = TRUE, cl)
  {
    if((missing(aD) | class(aD) != "alignmentData") & getLikes)
      stop("I can't assess the likelihoods of clustered data without an alignmentData object 'aD'")
     
    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)
    
    sD <- sD[order(as.factor(sD@segInfo$chr), sD@segInfo$start, -sD@segInfo$end),]


    dupStarts <- which(fastUniques(cbind(sD@segInfo$chr, sD@segInfo$start)))                                         
    emptyNulls <- with(sD@segInfo, data.frame(chr = chr[dupStarts], start = start[dupStarts] - leftSpace[dupStarts], end = start[dupStarts] - 1))
    emptyNulls <- emptyNulls[emptyNulls$start > 1, ]

    potnullD <- with(sD@segInfo, new("postSeg",
                                              data = rbind(matrix(0, ncol = ncol(sD), nrow = nrow(emptyNulls)),
                                                sD@data[leftSpace > 0 & rightSpace > 0,],
                                                sD@data[leftSpace > 0,],
                                                sD@data[rightSpace > 0,]),
                                     seglens = c(emptyNulls$end - emptyNulls$start + 1,
                                       (end - start + 1 + leftSpace + rightSpace)[leftSpace > 0 & rightSpace > 0],
                                       (end - start + 1 + leftSpace)[leftSpace > 0],
                                       (end - start + 1 + rightSpace)[rightSpace > 0]),
                                     libsizes = sD@libsizes,
                                     replicates = sD@replicates,
                                     annotation = data.frame(chr = c(emptyNulls$chr, chr[leftSpace > 0 & rightSpace > 0], chr[leftSpace > 0], chr[rightSpace > 0]),
                                       start = c(emptyNulls$start, (start - leftSpace)[leftSpace > 0 & rightSpace > 0], (start - leftSpace)[leftSpace > 0], start[rightSpace > 0]),
                                       end = c(emptyNulls$end, (end + rightSpace)[leftSpace > 0 & rightSpace > 0], end[leftSpace > 0], (end + rightSpace)[rightSpace > 0]),
                                       segType = "potnul",
                                       nullClass = c(rep("empty", nrow(emptyNulls)), rep("expBoth", sum(leftSpace > 0 & rightSpace > 0)), rep("expLeft", sum(leftSpace > 0)), rep("expRight", sum(rightSpace > 0))))
                                     ))

    potlociD <- with(sD@segInfo, new("postSeg",
                                     data = sD@data,
                                     seglens = end - start + 1,
                                     libsizes = sD@libsizes,
                                     replicates = sD@replicates,
                                     annotation = data.frame(chr = chr, start = start, end = end, segType = "potloc")))
    
    replicates <- potlociD@replicates

    message("Partitioning data...", appendLF = FALSE)
    
    clusterPosts <- lapply(unique(replicates), function(rep) {
      
      locDens <- log(colSums(t(potlociD@data[,replicates == rep,drop = FALSE]) / potlociD@libsizes[replicates == rep]) / sum(replicates == rep) / (potlociD@seglens))
      nulDens <- log(colSums(t(potnullD@data[,replicates == rep,drop = FALSE]) / potnullD@libsizes[replicates == rep]) / sum(replicates == rep) / (potnullD@seglens))

      nzlocs <- which(locDens != -Inf)
#      noOverlaps <- filterSegments(potlociD@annotation[nzlocs,], runif(length(nzlocs)), decreasing = TRUE)
      locCutoff <- bimodalSep(locDens[nzlocs], bQ = c(0.5, 1))

#      locCutoff <- bimodalSep(locDens[nzlocs], bQ = c(0.1, 0.9))

      locM <- as.numeric(locDens > locCutoff)
      nulM <- as.numeric(nulDens < locCutoff)
      
      znulls <- which(nulDens == -Inf)
      zerolens <- log(potnullD@seglens[znulls,])
      filZeros <- filterSegments(potnullD@annotation[znulls,], zerolens, decreasing = TRUE)

      lenCutoff <- bimodalSep(zerolens[filZeros])
      
      locM[log(potlociD@seglens) < lenCutoff & locM == 0] <- 0.5
      nulM[log(potnullD@seglens) < lenCutoff & nulM == 1] <- 0.5

      message(".", appendLF = FALSE)
      
      list(lociPost = locM, nullPost = nulM)
    })

    message("done.", appendLF = TRUE)

    potlociD@posteriors <- log(sapply(clusterPosts, function(x) x$lociPost))
    potnullD@posteriors <- log(sapply(clusterPosts, function(x) x$nullPost))

    if (is.null(subRegion)) {
      locSub <- 1:nrow(potlociD)
      nulSub <- 1:nrow(potnullD)
    } else {
      locSub <- sort(unique(c(unlist(apply(subRegion, 1, 
                                           function(sR) which(as.character(potlociD@annotation$chr) == as.character(sR[1]) &
                                                              potlociD@annotation$start >= as.numeric(sR[2]) &
                                                              potlociD@annotation$end <= as.numeric(sR[3])))))))
      nulSub <- sort(unique(c(unlist(apply(subRegion, 1, 
                                           function(sR) which(as.character(potnullD@annotation$chr) == as.character(sR[1]) &
                                                              potnullD@annotation$start >= as.numeric(sR[2]) &
                                                              potnullD@annotation$end <= as.numeric(sR[3])))))))

    }
    

    .processPosteriors(potlociD[locSub,], potnullD[nulSub,], aD = aD, lociCutoff = 1, nullCutoff = 1, getLikes = getLikes, cl = cl)
    
  }
