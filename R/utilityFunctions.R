
.filterSegments <- function(segs, orderOn, maxReport = Inf, ...)
  {    
    chrfilter <- function(chrsegs, suborderOn, ...)
      {
        message(".", appendLF = FALSE)
        cummaxEnd <- cummax(end(chrsegs))
        cumminStart <- cummax(start(chrsegs))
        
        fIns <- cbind(findInterval(start(chrsegs) - 0.5, cummaxEnd) + 1, 
                      findInterval(end(chrsegs), cumminStart))

        
        filtsegs <- order(suborderOn, ...)
    
        winsize <- 100
        ss <- 1
        
        filteredsegs <- c()
        
        while(length(filtsegs) > 0)
          {            
            chsam <- filtsegs[ss:min(ss + winsize -1, length(filtsegs))]
            oL <- lapply(chsam, function(x) setdiff((fIns[x,1L]:fIns[x,2L])[end(chrsegs)[fIns[x,1L]:fIns[x,2L]] >= start(chrsegs)[x] & start(chrsegs)[fIns[x,1L]:fIns[x,2L]] <= end(chrsegs)[x]], x))

            
            if(any(chsam %in% unlist(oL))) accsam <- !chsam %in% unlist(oL) else accsam <- rep(TRUE, length(chsam))
            accsam[1L] <- TRUE
            
            filteredsegs <- c(filteredsegs, chsam[accsam])

            if(length(filteredsegs) > maxReport) break()
            
            filtsegs <- setdiff(filtsegs, c(chsam[accsam], unlist(oL[accsam])))        
          }
        
        filteredsegs
      }



    newchrfilter <- function(chrsegs, suborderOn, ...)
      {
        message(".", appendLF = FALSE)
        
        filtsegs <- order(suborderOn, ...)
    
        winsize <- 10000
        ss <- 1
        
        filteredsegs <- c()
        
        while(length(filtsegs) > 0)
          {
            chsam <- filtsegs[ss:min(ss + winsize -1, length(filtsegs))]
            segsam <- chrsegs[chsam]
            chsam <- sort(chsam[findOverlaps(segsam, segsam, select = "first") == 1:length(chsam)])
            
            
            fIns <- cbind(findInterval(start(chrsegs) - 0.5, cummax(end(chrsegs[chsam]))) + 1,
                          findInterval(end(chrsegs) + 0.5, cummax(start(chrsegs[chsam]))))
            rejects <- which(fIns[,2] >= fIns[,1])            

            filteredsegs <- c(filteredsegs, chsam)

            if(length(filteredsegs) > maxReport) break()
            
            filtsegs <- filtsegs[!(filtsegs %in% rejects)]
          }
        
        filteredsegs
      }

    
    chrs <- levels(seqnames(segs))

    filtlist <- unlist(lapply(chrs, function(cc, orderOn, ...)
                              {
                                chrsamp <- which(seqnames(segs) == cc)
                                chrsegs <- ranges(segs)[chrsamp,]

                                rod <- order(start(chrsegs), end(chrsegs))
                                suborderOn <- (orderOn[chrsamp])[rod]
                                rodsegs <- chrsegs[rod,,drop = FALSE]
                                
                                chrfil <- newchrfilter(rodsegs, suborderOn, ...)
                                chrsamp[rod[chrfil]]
                              }, orderOn = orderOn, ...))
    
    filtlist
  }


.convertSegToLoci <- function(sD)
  {
    if(class(sD) == "segData") {
      lD <- new("lociData",
                data = matrix(ncol = ncol(sD), nrow = length(sD@coordinates)),
                libsizes = sD@libsizes,
                replicates = sD@replicates,
                coordinates = sD@coordinates,
                seglens = width(sD@coordinates))
      if(nrow(sD@data) > 0)
        lD@data <- sapply(1:ncol(sD), function(jj) as.integer(sD@data[,jj]))      
    } else if (class(sD) == "methSegs")
      {
        lD <- new("methData",
                  data = matrix(ncol = ncol(sD), nrow = length(sD@coordinates)),
                  pairData = matrix(ncol = ncol(sD), nrow = length(sD@coordinates)),
                  replicates = sD@replicates,
                  coordinates = sD@coordinates,
                  seglens = width(sD@coordinates))
        if(nrow(sD@Ts) > 0)
          lD@pairData <- sapply(1:ncol(sD), function(jj) as.integer(sD@Ts[,jj]))      
        if(nrow(sD@Cs) > 0)
          lD@data <- sapply(1:ncol(sD), function(jj) as.integer(sD@Cs[,jj]))      
      }
    if("seglens" %in% slotNames(sD) && nrow(sD@seglens) > 0)
      lD@seglens <- sapply(1:ncol(sD), function(jj) as.numeric(sD@seglens[,jj]))
    if(nrow(sD@locLikelihoods) > 0)
      lD@locLikelihoods <- sapply(1:ncol(sD@locLikelihoods), function(jj) as.double(sD@locLikelihoods[,jj]))
    lD
  }

.stratifySample <- function(stratSD, samplesize)
  {    
    sampData <- rowSums(sapply(1:ncol(stratSD), function(ii) as.integer(stratSD@data[,ii]) / stratSD@libsizes[ii]))
    sampData <- sampData[!duplicated(sampData)]
    
    sqnum <- length(sampData) / 1000
    squant <- quantile(sampData, 1:sqnum / sqnum)
    sqdup <- c(1, which(diff(squant) > min(1 / sD@libsizes) / 10))
    z <- cbind(as.numeric(squant[sqdup]), c(as.numeric(squant[sqdup[-1]]), max(sampData)))
    z[1,1] <- -Inf

    sy <- do.call("rbind",
                  lapply(1:nrow(z), function(ii) {
                    w <- z[ii,]
                    inbetweener <- which(sampData > w[1] & sampData <= w[2])
                    samplenum <- min(length(inbetweener), ceiling(samplesize / nrow(z)))
                    cbind(sample(inbetweener, size = samplenum, replace = FALSE), length(inbetweener) / samplenum)
                  })
                  )
    
    weights <- sy[,2]
    sy <- sy[,1]

    y <- sD@data[sy,,drop = FALSE]
    if(lensameFlag) seglensy <- seglens[sy] else seglensy <- seglens[sy,,drop = FALSE]
    
    ordData <- do.call(order, as.data.frame(cbind(y, seglensy)))
    y <- y[ordData,, drop = FALSE]
    weights <- weights[ordData]
    if(lensameFlag) seglensy <- seglensy[ordData] else seglensy <- seglensy[ordData,,drop = FALSE]
    
    dups <- c(1, which(rowSums((cbind(y,seglensy))[-1,] == (cbind(y, seglensy))[-nrow(y),]) != ncol(cbind(y, seglensy))) + 1)
    copies <- diff(c(dups, nrow(y) + 1))

    sy <- cbind(sampled = sy[ordData], representative = rep(1:length(copies), copies))

    y <- y[dups,, drop = FALSE]
    if(lensameFlag) seglensy <- seglensy[dups] else seglensy <- seglensy[dups,, drop = FALSE]
  }


.splitSD <- function(sD, largeness = 1e8)
  {
    winsize <- round(nrow(sD) / ceiling(prod(dim(sD)) / largeness))
    
    chunkSD <- values(findChunks(sD@coordinates, gap = 0, checkDuplication = FALSE))$chunk
    chunkWindows <- cbind(1:length(runLength(chunkSD)), cumsum(runLength(chunkSD)))    
    dupWin <- which(!duplicated(ceiling(chunkWindows[,2] / winsize)))
    dupWin <- cbind(dupWin, c(dupWin[-1] - 1, nrow(chunkWindows)))
    windowChunks <- lapply(1:nrow(dupWin), function(ii) unique(chunkSD)[chunkWindows[dupWin[ii,1]:dupWin[ii,2],1L]])
    
    sDsplit <- lapply(windowChunks, function(wc) {
      sD[which(chunkSD %in% wc),]
    })

    sDsplit
  }


.fastUniques <- function(x){
  if (nrow(x) > 1) {
    return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
  } else return(TRUE)
}

.mergeListLoci <- function(splitSeg) {
  mergeSeg <- new(class(splitSeg[[1]]),
                  locLikelihoods = do.call("rbind", lapply(splitSeg, function(x) x@locLikelihoods)),
                  coordinates = do.call("c", lapply(splitSeg, function(x) x@coordinates)),
                  data = do.call("rbind", lapply(splitSeg, function(x) x@data)),
                  replicates = splitSeg[[1]]@replicates,
                  libsizes = splitSeg[[1]]@libsizes,
                  groups = splitSeg[[1]]@groups,
                  annotation = do.call("rbind", lapply(splitSeg, function(x) x@annotation)),
                  seglens = do.call("rbind", lapply(splitSeg, function(x) x@seglens)))
  if(class(mergeSeg) == "methData")
    mergeSeg@pairData <- do.call("rbind", lapply(splitSeg, function(x) x@pairData))
  mergeSeg
}

.mergeSegData <- function(splitSeg) {

  message("Merging...", appendLF = FALSE)
  locLikelihoods <- do.call("DataFrame", lapply(1:ncol(splitSeg[[1]]@locLikelihoods), function(ii) {
    message(".", appendLF = FALSE)
    do.call("c", lapply(splitSeg, function(x) x@locLikelihoods[,ii]))
  }))
  data <- do.call("DataFrame", lapply(1:ncol(splitSeg[[1]]@data), function(ii) {
    message(".", appendLF = FALSE)
    do.call("c", lapply(splitSeg, function(x) x@data[,ii]))
  }))

  
  mergeSeg <- new("segData",                                                 
                  locLikelihoods = locLikelihoods,
                  coordinates = do.call("c", lapply(splitSeg, function(x) x@coordinates)),
                  data = data,
                  replicates = splitSeg[[1]]@replicates,
                  libsizes = splitSeg[[1]]@libsizes)
  message("done.")
  mergeSeg
}
