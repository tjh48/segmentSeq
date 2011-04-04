mergeSD <- function(..., aD, replicates, gap = 200, cl = NULL) {

  updateSD <- function(nsD, aD, cl = NULL)
    {
      if(is.unsorted(nsD@segInfo$chr)) nsD <- nsD[order(nsD@segInfo$chr),]
      updateSegs <- which(rowSums(is.na(nsD@data)) > 0)
      sI <- nsD@segInfo[updateSegs,]
      nsD@data[updateSegs,]    <- do.call("rbind",
                                          lapply(unique(sI$chr), function(chr) {
                                            chrAD <- aD[aD@alignments$chr == chr,]
                                            crsi <- sI[sI$chr == chr,]
                                            gap <- which(abs(diff(crsi$start)) > 1e4)
                                            splits <- cbind(c(1, gap + 1), c(gap, nrow(crsi)))
                                            do.call("rbind", lapply(1:nrow(splits), function(ii) {
                                              x <- as.numeric(splits[ii,])
                                              raD <- chrAD[chrAD@alignments$start >= crsi$start[x[1]] & chrAD@alignments$end <= crsi$end[x[2]],]
                                              getCounts(crsi[x[1]:x[2],], raD, cl = cl)
                                            }))
                                          }))
      nsD
    }

  
            binds <- list(...)
                                        #  gap <- 500
            maxloclen <- Inf
            
            fastUniques <- function(x)
              if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)
            
            binds <- lapply(binds, function(x) x[with(x@segInfo, order(as.factor(chr), start, end)),])
            subinds <- lapply(binds, function(x) subset(x@segInfo, subset = fastUniques(subset(x@segInfo, select = c("chr", "start"))), select = c("chr", "start", "end")))
            mergeSegs <- do.call("rbind", subinds)
            
            chrs <- unique(as.character(mergeSegs$chr))
            
            chrs <- do.call("rbind", lapply(binds, function(x) x@chrs))
            chrs <- chrs[!duplicated(chrs),]
            
            segs <- do.call("rbind", lapply(1:nrow(chrs), function(cc)
                                            {
                                        #      if(verbose){
                                        #        message("Chromosome: ", chrs$chr[cc])
                                        #        message("Finding start-stop co-ordinates...")
                                        #      }
                                              if(any(mergeSegs$chr == chrs$chr[cc]))
                                                {
                                                  chrTags <- mergeSegs[mergeSegs$chr == chrs$chr[cc],, drop = FALSE]          
                                                  
                                                  if(any(chrTags[,3L] > chrs$len[cc]))
                                                    warning(paste("Chromosome", chrs[cc], "has tags which extend over the given chromsome length."))
                                                  
                                                  chrSS <- chrTags[, 2:3, drop = FALSE]
                                                  if(nrow(chrSS) > 1)
                                                    {
                                                      chrSS <- chrSS[order(chrSS[,1L], chrSS[,2L]),]
                                                      chrSS <- chrSS[fastUniques(chrSS),]
                                                      chrmax <- matrix(c(cummax(chrSS[,1]), cummax(chrSS[,2])), ncol = 2)
                                                      ch <- which(chrmax[-1L,1L] > chrmax[-nrow(chrmax),2L])
                                                      startstop <- cbind(start = as.integer(c(min(chrmax[,1L]), chrmax[ch + 1L,1L])),
                                                                         end = as.integer(c(chrmax[ch,2L], max(chrmax[,2L]))))
                                                    } else if(nrow(chrSS) == 1) startstop <- cbind(starts = chrSS$start, ends = chrSS$end) else startstop <- matrix(nrow = 0, ncol = 0)   
                                                  gaps <- c(which(startstop[-nrow(startstop),2L] < startstop[-1L, 1L] - gap), nrow(startstop))
                                                  
                                        #          if(verbose)
                                        #            message("Defining potential subsegments...", appendLF = FALSE)
                                                  cgap <- cbind(1:nrow(startstop), rep(gaps, diff(c(0, gaps))))
                                                  csegs <- data.frame(start = rep(startstop[,1L], cgap[,2] - cgap[,1] + 1), end = startstop[unlist(lapply(1:nrow(cgap), function(ii) cgap[ii,1]:cgap[ii,2])), 2L])
                                        #          message("done!")
                                                  if(any(csegs[,2L] - csegs[,1L] + 1L > maxloclen))
                                                    {
                                                      csegs <- csegs[csegs[,2L] - csegs[,1L] + 1 <= maxloclen,]
                                                      csegs <- rbind(csegs, startstop[(!(startstop[,1L] %in% csegs[,1L])),])
                                                    }
                                                  
                                                  csegs <- csegs[order(csegs[,1L], csegs[,2L]),,drop = FALSE]
                                                  csegs <- cbind(csegs,
                                                                 leftSpace = (csegs[,1L] - c(0L, startstop[,2])[findInterval(csegs[,1L], startstop[,2L]) + 1L]) - 1L,
                                                                 rightSpace = c(startstop[,1L], chrs$len[cc] + 1L)[findInterval(csegs[,2L], startstop[,1L]) + 1L] - csegs[,2L] - 1L)
                                                  csegs <- cbind(chr = chrs$chr[cc], csegs)
                                                  
                                                } else csegs <- NULL
                                              csegs
                                            }))
            
            data <- do.call("rbind", lapply(chrs$chr, function(chr) {
              csegs <- segs[segs$chr == chr,]
              extractData <- function(x)
                {
                  chrx <- x[x@segInfo$chr == chr,]
                  chrx <- chrx[order(chrx@segInfo$end, chrx@segInfo$start),]
                  xends <- findInterval(csegs$end, chrx@segInfo$end)
                  xends[xends == 0] <- NA
                  xends[chrx@segInfo$start[xends] < csegs$start] <- NA
                  xstarts <- findInterval(chrx@segInfo$end[xends] - 0.5, chrx@segInfo$end) + 1
                  xregs <- cbind(xstarts, xends)
                  xdata <- matrix(NA, nrow = nrow(csegs), ncol = ncol(x))
                  xdata[is.na(xregs[,1]),] <- 0        
                  xending <- chrx@segInfo$end[xends]
                  
                  chrx <- chrx[order(chrx@segInfo$start, chrx@segInfo$end),]
                  xstarts <- findInterval(csegs$start - 0.5, chrx@segInfo$start) + 1
                  xstarting <- chrx@segInfo$start[xstarts]
                  xregs <- cbind(xstarting, xending, 1:nrow(csegs), NA, deparse.level = 0)[is.na(xdata[,1]),]
                  xlink <- cbind(as.matrix(subset(chrx@segInfo, select = c("start", "end"))), NA, 1:nrow(chrx), deparse.level = 0)
                  colnames(xlink) <- NULL
                  xlink <- rbind(xlink, xregs)
                  xlink <- xlink[order(xlink[,1], xlink[,2]),]
                  xuniques <- which(fastUniques(xlink[,1:2]))
                  xrow <- rep(xlink[xuniques,4], diff(c(xuniques, nrow(xlink))))
                  
                  xdata[xlink[which(!is.na(xlink[,3])),3],] <- chrx@data[xrow[which(!is.na(xlink[,3]))],]
                  xdata
                }
              data <- do.call("cbind", lapply(binds, extractData))
            }))
            
            z <- new("segData")
            z@data <- data
            z@segInfo <- segs
            z@libsizes <- unlist(lapply(binds, function(x) x@libsizes))
            z@chrs <- chrs
            z@replicates <- as.integer(replicates)
            
  updateSD(z, aD, cl = cl)
          }
          


setMethod("[", "segData", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      x@data <- x@data[i,, drop = FALSE]
      x@segInfo <- x@segInfo[i,,drop = FALSE]
    }

  if(!missing(j))
    {
      x@replicates <- as.integer(x@replicates[j])
      x@data <- x@data[,j,drop = FALSE]
      x@libsizes <- x@libsizes[j]
    }

  x
})

setMethod("dim", "segData", function(x) {
  dim(x@data)
})


setValidity("segData", function(object) {
  validmess <- c()
  valid <- TRUE
  if(length(object@libsizes) != ncol(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@libsizes' slot must equal number of columns of '@data' slot.")
    }
  if(nrow(object@rightData) > 0 & any(dim(object@rightData) != dim(object@data)))
    {
      valid <- FALSE
      validmess <- c(validmess, "If '@rightData' slot is non-empty, the dimensions of '@rightData' slot must equal those of the '@data' slot.")
    }
  if(nrow(object@segInfo) != nrow(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Number of rows of '@segInfo' slot not same as '@data' slot.")
    }
  if(length(object@replicates) != ncol(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@replicates' slot must equal number of columns of '@data' slot.")
    }
  if(!all(apply(object@data, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(!all(apply(object@data, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(valid) return(valid) else validmess
})


setMethod("show", "segData", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  cat('\nSlot "data":\n')
  if(nrow(object) > 5)
    {
      print(object@data[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@data)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
  cat('\nSlot "replicates":\n')
  print(object@replicates)
  cat('\nSlot "segInfo":\n')
  if(nrow(object@segInfo) > 5)
    {
      print(object@segInfo[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@segInfo)
  if(length(object@priorType) > 1)
    {
      cat('Slot "priors":\n')
      cat(paste('Priors are of type:', object@priorType), '\n')
    }
})


