filterSegments <- function(segs, orderOn, ...)
  {    
    chrfilter <- function(segs, suborderOn, ...)
      {
        cummaxEnd <- cummax(segs$end)
        cumminStart <- cummax(segs$start)
        
        fIns <- cbind(findInterval(segs$start, cummaxEnd) + 1, 
                      findInterval(segs$end, cumminStart))
        
        filtsegs <- order(suborderOn, ...)
    
        winsize <- 10
        ss <- 1
        
        filteredsegs <- c()
        
        while(length(filtsegs) > 0)
          {
            chsam <- filtsegs[ss:min(ss + winsize -1, length(filtsegs))]
            oL <- sapply(chsam, function(x) setdiff((fIns[x,1L]:fIns[x,2L])[segs$end[fIns[x,1L]:fIns[x,2L]] >= segs$start[x] & segs$start[fIns[x,1L]:fIns[x,2L]] <= segs$end[x]], x))
            if(any(chsam %in% unlist(oL))) accsam <- !chsam %in% unlist(oL) else accsam <- rep(TRUE, length(chsam))
            accsam[1L] <- TRUE
            
            filteredsegs <- c(filteredsegs, chsam[accsam])
            
            filtsegs <- setdiff(filtsegs, c(chsam[accsam], unlist(oL[accsam])))        
          }
        
        filteredsegs
      }

    chrs <- unique(segs$chr)

    filtlist <- unlist(lapply(chrs, function(cc, orderOn, ...)
                              {
                                chrsamp <- segs$chr == cc
                                chrsegs <- subset(segs, chrsamp, select = c(start, end))

                                rod <- order(chrsegs$start, chrsegs$end)
                                suborderOn <- (orderOn[chrsamp])[rod]
                                rodsegs <- chrsegs[rod,,drop = FALSE]
                                
                                chrfil <- chrfilter(rodsegs, suborderOn, ...)
                                which(chrsamp)[rod[chrfil]]
                              }, orderOn = orderOn))
    
    filtlist
  }
