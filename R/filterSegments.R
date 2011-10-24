.filterSegments <- function(segs, orderOn, ...)
  {    
    chrfilter <- function(chrsegs, suborderOn, ...)
      {
        cummaxEnd <- cummax(end(chrsegs))
        cumminStart <- cummax(start(chrsegs))
        
        fIns <- cbind(findInterval(start(chrsegs) - 0.5, cummaxEnd) + 1, 
                      findInterval(end(chrsegs), cumminStart))
        
        filtsegs <- order(suborderOn, ...)
    
        winsize <- 10
        ss <- 1
        
        filteredsegs <- c()
        
        while(length(filtsegs) > 0)
          {
            chsam <- filtsegs[ss:min(ss + winsize -1, length(filtsegs))]
            oL <- lapply(chsam, function(x) setdiff((fIns[x,1L]:fIns[x,2L])[end(chrsegs)[fIns[x,1L]:fIns[x,2L]] >= start(chrsegs)[x] & start(chrsegs)[fIns[x,1L]:fIns[x,2L]] <= end(chrsegs)[x]], x))
            if(any(chsam %in% unlist(oL))) accsam <- !chsam %in% unlist(oL) else accsam <- rep(TRUE, length(chsam))
            accsam[1L] <- TRUE
            
            filteredsegs <- c(filteredsegs, chsam[accsam])
            
            filtsegs <- setdiff(filtsegs, c(chsam[accsam], unlist(oL[accsam])))        
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
                                
                                chrfil <- chrfilter(rodsegs, suborderOn, ...)
                                chrsamp[rod[chrfil]]
                              }, orderOn = orderOn, ...))
    
    filtlist
  }
