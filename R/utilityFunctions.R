% modification on git from copied files

thresholdFinder <- function(method, aM, subset, minprop = 0.05, bootstrap = 100, abstol = 1e-4, verbose = FALSE, cl = NULL, processAD.args = list(), heuristicSeg.args = list()) {
    if(!missing(subset)) aMS <- aM[subset,] else aMS <- aM
    nD <- normaliseNC(aMS)

    findProp <- function(method, nD, loci, minprop, verbose, cl, processAD.args, heuristicSeg.args) {
    
    repSep <- sapply(levels(nD@replicates), function(rep) {
        if(inherits(nD, "alignmentMeth")) {
            Cs <- rowSums(nD@Cs[,nD@replicates == rep, drop = FALSE])
            Ts <- rowSums(nD@Ts[,nD@replicates == rep, drop = FALSE])
            if(!is.null(loci))
                {
                    fO <- findOverlaps(loci@coordinates, nD@alignments)                    
                    ov <- split(subjectHits(fO), factor(queryHits(fO), levels = 1:nrow(loci)))
                if(any(loci@locLikelihoods < 0)) ll <- exp(loci@locLikelihoods) else ll <- loci@locLikelihoods
                sel <- ll[,rep] == 1 | rowSums(ll) == 0
                llsamp <- unlist(ov[sel])
                Cs <- Cs[llsamp]; Ts <- Ts[llsamp]                
                ovec <- c(); ovec[unlist(ov)] <- rep(1:length(ov), sapply(ov, length))
                ovec <- ovec[unlist(ov[sel])]
            } else ovec <- rep(1,length(Cs))            
        } else if(inherits(nD, "lociData")) {
            loci <- NULL
            Cs <- rowSums(nD@data[,nD@replicates == rep,1,drop = FALSE])
            Ts <- rowSums(nD@data[,nD@replicates == rep,2,drop = FALSE])
            if(nrow(mD@locLikelihoods) == length(Cs)) {
                if(any(mD@locLikelihoods < 0)) ll <- exp(mD@locLikelihoods) else ll <- mD@locLikelihoods
                sel <- ll[,rep] == 1 | rowSums(ll) == 0
                Cs <- Cs[sel]
                Ts <- Ts[sel]
            }
            ovec <- rep(1, length(Cs))
        }        
        
        if(method == "varsum") return(bimodalSeparator(Cs / (Cs + Ts)))
        if(method == "minden") {
            z <- Cs / (Cs + Ts); z <- z[!is.na(z)]
            dd <- density(z, from = 0, to = 1)
            xes <- dd$x[which(dd$y < c(dd$y[-1], Inf) & dd$y < c(Inf, dd$y[-length(dd$y)]))]
            return(min(xes[xes > quantile(z, minprop)]))
        }
        if(method == "beta" | method == "betaVarsum") {
            meanbeta <- function(x, Cs, Ts, ovec) {
                db <- dbeta(x, 0.5 + Cs, 0.5 + Ts)
                mean(sapply(split(db, ovec), mean))
            }
                                        #nCs <- Cs[Cs >0]; nTs <- Ts[Cs > 0]
            nCs <- Cs; nTs <- Ts
            
            divisions <- 1000; samp <- 100000
            basamp <- sample(1:length(nCs), min(samp, length(nCs)))
            
            if(!is.null(cl)) {
                mb <- parSapply(cl, seq_len(divisions - 1) / divisions, meanbeta, Cs = nCs[basamp], Ts = nTs[basamp], ovec = ovec[basamp])
            } else mb <- sapply(seq_len(divisions - 1) / divisions, meanbeta, Cs = nCs[basamp], Ts = nTs[basamp], ovec = ovec[basamp])
                
                                        #            lowdens <- which(mb / divisions < 1 / min(samp, length(nCs)))
#            suppressWarnings(lowdens <- min(lowdens[lowdens / divisions > minprop]) / divisions)
                                        #            if(lowdens < 1) return(lowdens)
            maxima <- c(which(mb > c(mb[-1], -Inf) & mb > c(-Inf, mb[-length(mb)])) / divisions, 1)
            minima <- c(which(mb < c(mb[-1], Inf) & mb < c(Inf, mb[-length(mb)])) / divisions)
            if(any(minima > minprop)) minmin <- min(minima[which(minima > minprop)]) else minmin <- max(minima)

            if(method == "beta")
                return(optimise(meanbeta, interval = c(max(maxima[maxima < minmin]), min(maxima[maxima > minmin])), Cs = nCs, Ts = nTs, ovec = ovec)$minimum)                        
            if(method == "betaVarsum") {
                fI <- findInterval(runif(1000000, min = 0, max = 1) , cumsum(mb) / divisions) / divisions
                return(bimodalSeparator(fI[fI > max(maxima[maxima < minmin]) & fI < min(maxima[maxima > minmin])]))
            }
        }
        if(method == "abc") {
            divisions <- 1000; samp <- 1000
            #nCs <- Cs[Cs >0]; nTs <- Ts[Cs > 0]
            nCs <- Cs; nTs <- Ts
            abcs <- do.call("rbind", lapply(sample(1:length(nCs), min(samp,length(nCs))), function(jj) {
                rs1 <- runif(1000, max(0, 0.5 + nCs[jj] - 5), 0.5+nCs[jj] + 5)
                rs2 <- runif(1000, max(0, 0.5 + nTs[jj] - 5), 0.5+nTs[jj] + 5)
                rnCs <- rbinom(1000, size = round(nCs[jj] + nTs[jj]), rbeta(1000, rs1, rs2))
                equal <- abs(rnCs - nCs[jj]) < 0.01 * (nCs[jj] + nTs[jj])                
                rowMeans(matrix(dbeta(rep(seq_len(divisions - 1) / divisions, sum(equal)), rep(rs1[equal], each = divisions - 1), rep(rs2[equal], each = divisions - 1)), ncol = sum(equal)))
            }))
            ma <- colMeans(abcs)
            maxima <- c(which(ma > c(ma[-1], -Inf) & ma > c(-Inf, ma[-length(ma)])) / divisions, 1)
            minima <- c(which(ma < c(ma[-1], -Inf) & ma < c(-Inf, ma[-length(ma)])) / divisions, 1)
            return(min(minima[minima > minprop]))
        }
    })
    return(repSep)
}


    
    if(bootstrap > 1) sD <- do.call("processAD", c(list(aD = aMS, verbose = FALSE, getCounts = TRUE, cl = cl), processAD.args))
    propVec <- c()
    for(pp in 1:bootstrap) {
        if(pp == 1) mD <- NULL else mD <- .inferMethNulls(hS, aMS, cl = cl)
        repSep <- findProp(method, nD, loci = mD, minprop = minprop, verbose = verbose, cl = cl)
        prop = mean(repSep)
        if(verbose) message("Currently estimated threshold: ", prop)
        if(any(abs(propVec - prop) < abstol) | pp == bootstrap) break()
        propVec[pp] <- prop
        hS <- do.call("heuristicSeg", c(list(sD = sD, aD = aMS, prop = prop, verbose = FALSE, cl = cl, getLikes = FALSE), heuristicSeg.args))
    }

    if(length(repSep) > 1 && sd(repSep) > 0.1) warning(paste("Standard deviation of automatically inferred methylation thresholds over replicates is greater than", sd(repSep), "; you may wish to examine the distribution of methylation for each replicate and specify a threshold manually."))    
    if(verbose) message("Automatically determined threshold for methylation locus; ", prop)
    prop
}

.printIRangesMatrix <- function(x)
  {
    cat("Matrix with ", nrow(x), " rows.\n")
    if(nrow(x) > 10)
      {      
        print(rbind(as.data.frame(x)[1:5,,drop = FALSE], matrix("...", nrow = 1, ncol = ncol(x), dimnames = list("...", colnames(x))), as.data.frame(x)[nrow(x) + (-4):0,,drop = FALSE]))
      } else print(x)  
  }

.printLocLikes <- function(locLike) {
  if(class(locLike) == "DataFrame") locLike <- .DF2matrix(locLike)
  
  if(any(exp(locLike) > 1, na.rm = TRUE))
    {
      cat('\nSlot "locLikelihoods":\n')
      modFunction <- identity
    } else {          
      cat('\nSlot "locLikelihoods" (stored on log scale):\n')        
      modFunction <- exp
    }
  .printIRangesMatrix(signif(modFunction(locLike), 5))
}


.filterSegments <- function(segs, orderOn, maxReport = Inf, ...)
  {
    chrfilter <- function(chrsegs, suborderOn, ...)
      {
        #message(".", appendLF = FALSE)
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
          #message(".", appendLF = FALSE)
        
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

    
    chrs <- unique(seqnames(segs))

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
        lD@data <- sD@data
    } else if (class(sD) == "segMeth")
      {
        lD <- new("methData",
                  data = list(matrix(ncol = ncol(sD), nrow = length(sD@coordinates)),
                    matrix(ncol = ncol(sD), nrow = length(sD@coordinates))),
                  replicates = sD@replicates,
                  coordinates = sD@coordinates,
                  seglens = width(sD@coordinates))                  
    
        if(nrow(sD@Ts) > 0 & nrow(sD@Cs) > 0)
          lD@data <- array(c(sD@Cs, sD@Ts), c(dim(sD@Cs), 2))
      }
    if("seglens" %in% slotNames(sD) && nrow(sD@seglens) > 0)
      lD@seglens <- sD@seglens
    lD@locLikelihoods <- .DF2matrix(sD@locLikelihoods)
      
    lD
  }

#.stratifySample <- function(stratSD, samplesize)
#  {    
#    sampData <- rowSums(sapply(1:ncol(stratSD), function(ii) as.integer(stratSD@data[,ii]) / stratSD@libsizes[ii]))
#    sampData <- sampData[!duplicated(sampData)]
#    
#    sqnum <- length(sampData) / 1000
#    squant <- quantile(sampData, 1:sqnum / sqnum)
#    sqdup <- c(1, which(diff(squant) > min(1 / sD@libsizes) / 10))
#    z <- cbind(as.numeric(squant[sqdup]), c(as.numeric(squant[sqdup[-1]]), max(sampData)))
#    z[1,1] <- -Inf

#    sy <- do.call("rbind",
#                  lapply(1:nrow(z), function(ii) {
#                    w <- z[ii,]
#                    inbetweener <- which(sampData > w[1] & sampData <= w[2])
#                    samplenum <- min(length(inbetweener), ceiling(samplesize / nrow(z)))
#                    cbind(sample(inbetweener, size = samplenum, replace = FALSE), length(inbetweener) / samplenum)
#                  })
#                  )
#    
#    weights <- sy[,2]
#    sy <- sy[,1]
#    y <- sD@data[sy,,drop = FALSE]
#    if(lensameFlag) seglensy <- seglens[sy] else seglensy <- seglens[sy,,drop = FALSE]
#    
#    ordData <- do.call(order, as.data.frame(cbind(y, seglensy)))
#    y <- y[ordData,, drop = FALSE]
#    weights <- weights[ordData]
#    if(lensameFlag) seglensy <- seglensy[ordData] else seglensy <- seglensy[ordData,,drop = FALSE]
#    
#    dups <- c(1, which(rowSums((cbind(y,seglensy))[-1,] == (cbind(y, seglensy))[-nrow(y),]) != ncol(cbind(y, seglensy))) + 1)
#    copies <- diff(c(dups, nrow(y) + 1))
#
#    sy <- cbind(sampled = sy[ordData], representative = rep(1:length(copies), copies))
#
#    y <- y[dups,, drop = FALSE]
#    if(lensameFlag) seglensy <- seglensy[dups] else seglensy <- seglensy[dups,, drop = FALSE]
#  }


.splitSD <- function(sD, chunkSD, largeness = 1e8)
  {
      winsize <- ceiling(nrow(sD) / ceiling(prod(dim(sD)) / largeness))
    
    if(missing(chunkSD)) chunkSD <- findChunks(sD@coordinates, gap = 0, checkDuplication = FALSE, justChunks = TRUE)
      chunkWindows <- cbind(1:length(runLength(chunkSD)), cumsum(runLength(chunkSD)))    
      winchunks <- ceiling(chunkWindows[,2] / winsize)
      dupWin <- sort(c(which(!duplicated(winchunks)), which(diff(winchunks) > 1) + 2))
      dupWin <- cbind(dupWin, c(dupWin[-1] - 1, nrow(chunkWindows)))
      dupWin <- dupWin[dupWin[,2] >= dupWin[,1],,drop = FALSE]
      
      windowChunks <- lapply(1:nrow(dupWin), function(ii) {
          unique(chunkSD)[chunkWindows[dupWin[ii,1]:dupWin[ii,2],1L]]
      })
      
      sDsplit <- lapply(windowChunks, function(wc) {
          which(chunkSD %in% wc)
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
                  data = do.call("abind", c(lapply(splitSeg, function(x) x@data), list(along = 1))),
                  replicates = splitSeg[[1]]@replicates,
                  sampleObservables = splitSeg[[1]]@sampleObservables,
                  groups = splitSeg[[1]]@groups,
                  annotation = do.call("rbind", lapply(splitSeg, function(x) x@annotation))
                  )
    if(length(splitSeg[[1]]@cellObservables) > 0) {
      mergeSeg@cellObservables <- lapply(1:length(splitSeg[[1]]@cellObservables), function(ii)
                                         do.call("abind", c(lapply(splitSeg, function(x) x@cellObservables[[ii]]), list(along = 1))))
      names(mergeSeg@cellObservables) <- names(splitSeg[[1]]@cellObservables)
    }
  if(length(splitSeg[[1]]@rowObservables) > 0) {
    mergeSeg@rowObservables <- lapply(1:length(splitSeg[[1]]@rowObservables), function(ii)
                                      do.call("abind", c(lapply(splitSeg, function(x) x@rowObservables[[ii]]), list(along = 1)))
                                              )
      names(mergeSeg@rowObservables) <- names(splitSeg[[1]]@rowObservables)
    }
  mergeSeg
}

.mergeSegData <- function(splitSeg) {

  message("Merging...", appendLF = FALSE)
  #locLikelihoods <- do.call("cbind", lapply(1:ncol(splitSeg[[1]]@locLikelihoods), function(ii) {
  #  message(".", appendLF = FALSE)
  #  do.call("c", lapply(splitSeg, function(x) x@locLikelihoods[,ii]))
  #}))

  #data <- do.call("cbind", lapply(1:ncol(splitSeg[[1]]@data), function(ii) {
  #  message(".", appendLF = FALSE)
  #  do.call("c", lapply(splitSeg, function(x) x@data[,ii]))
  #}))
  
  mergeSeg <- new(class(splitSeg[[1]]))
  mergeSeg@locLikelihoods = do.call("rbind", lapply(splitSeg, function(x) x@locLikelihoods))
  mergeSeg@coordinates = do.call("c", lapply(splitSeg, function(x) x@coordinates))
  mergeSeg@replicates = splitSeg[[1]]@replicates
  
  message("done.")
  mergeSeg
}


normaliseNC <- function(mD, nonconversion) {
    if(inherits(mD, "alignmentMeth")) {
        nonconversion <- mD@nonconversion
        Cs <- mD@Cs
        Ts <- mD@Ts
    } else if(inherits(mD, "countData")) {
       if(missing(nonconversion) & "nonconversion" %in% names(mD@sampleObservables))
           nonconversion <- mD@sampleObservables$nonconversion
       Cs <- mD@data[,,1]
       Ts <- mD@data[,,2]
    }
    
    ks <- pmin(Cs, t(t(Ts) * nonconversion / (1 - nonconversion)))

    Cs <- Cs - ks
    Ts <- Ts + ks
    Cs[Cs < 0] <- 0
    if(inherits(mD, "alignmentMeth")) {
        mD@Cs <- Cs
        mD@Ts <- Ts
    } else {
        mD@data[,,1] <- Cs
        mD@data[,,2] <- Ts
    }
    return(mD)
}

.methFunction <- function(methD, prop, locCutoff, nullP = FALSE, prior = c(0.5, 0.5))
  {
    if(nrow(methD) == 0) return(matrix(nrow = 0, ncol = length(levels(methD@replicates))))
    nD <- normaliseNC(methD)

    if(inherits(methD, "alignmentMeth")) {
        combCs <- as.vector(sapply(levels(methD@replicates), function(rep) rowSums(nD@Cs[,methD@replicates == rep,drop = FALSE])))
        combTs <- as.vector(sapply(levels(methD@replicates), function(rep) rowSums(nD@Ts[,methD@replicates == rep,drop = FALSE])))
    } else if(inherits(methD, "countData")) {
        combCs <- as.vector(sapply(levels(methD@replicates), function(rep) rowSums(nD@data[,methD@replicates == rep,1,drop = FALSE])))
        combTs <- as.vector(sapply(levels(methD@replicates), function(rep) rowSums(nD@data[,methD@replicates == rep,2,drop = FALSE])))
    }
    
    combDat <- cbind(combCs, combTs)
    rodDat <- order(combDat[,1], combDat[,2])

    unqDat <- .fastUniques(combDat[rodDat,])

    redDat <- combDat[rodDat[unqDat],]
    redTF <- rep(NA, nrow(redDat))
    plow <- redDat[,1] / rowSums(redDat) < prop

    if(!is.na(locCutoff)) {
                                        # Beta-distribution here is derived from being the conjugate of a binomial
    if(!nullP) {      
      redTF[plow] <- FALSE
      if(prior[2] == 0)
        redTF[redDat[,1] > 0 & redDat[,2] == 0] <- TRUE
      if(prior[1] == 0)
        redTF[redDat[,1] == 0 & redDat[,2] > 0] <- FALSE            
      pbetas <- pbeta(prop, prior[1] + redDat[which(!plow),1], prior[2] + redDat[which(!plow),2], lower.tail = FALSE)
      redTF[which(!plow)] <- pbetas > locCutoff      
    } else {
      redTF[!plow] <- FALSE
      if(prior[2] == 0)
        redTF[redDat[,1] > 0 & redDat[,2] == 0] <- FALSE
      if(prior[1] == 0)
        redTF[redDat[,1] == 0 & redDat[,2] > 0] <- TRUE                  
      pbetas <- pbeta(prop, prior[1] + redDat[which(plow),1], prior[2] + redDat[which(plow),2], lower.tail = TRUE)
      redTF[which(plow)] <- pbetas > locCutoff      
    }
  } else {
    redTF[plow] <- FALSE
    redTF[!plow] <- TRUE
  }

    TF <- rep(redTF, diff(c(which(unqDat), length(rodDat) + 1)))
    TF[rodDat] <- TF
    locM <- matrix(TF, ncol = length(levels(methD@replicates)))
    colnames(locM) <- levels(methD@replicates)
    #locM <- .matrix2Rle(locM)
    locM
  }

.matrix2Rle <- function(x) {
  dx <- do.call("DataFrame", apply(x, 2, Rle))
    if(!is.null(colnames(x))) colnames(dx) <- colnames(x) else colnames(dx) <- paste("X", 1:ncol(x), sep = "")
  dx
}

.DF2matrix <- function(x) {
  if(any(dim(x)) == 0) return(array(dim = dim(x)))
  dx <- do.call("cbind", lapply(as.list(x), as.vector))
  colnames(dx) <- colnames(x)
  dx
}

.rowSumDF <- function(x, na.rm = TRUE) {
  z <- as.list(x)
  if(na.rm) {
    nas <- Reduce("intersect", lapply(z, function(zz) which(is.na(zz))))
    z <- lapply(z, function(zz) {
      zz[which(is.na(zz))] == 0
      zz
    })
  }
  sumz <- Reduce("+", z)
  if(na.rm && length(nas) > 0) sumz[nas] <- NA
  sumz
}
