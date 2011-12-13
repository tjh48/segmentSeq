getOverlaps <- function(coordinates, segments, overlapType = "overlapping", whichOverlaps = TRUE, cl)
  {
    if(missing(cl)) cl <- NULL
    
    if(overlapType == "overlapping") {
      segord <- order(as.factor(seqnames(segments)), start(segments), end(segments))
    } else if (overlapType == "contains") {
      segord <- order(as.factor(seqnames(segments)), end(segments), start(segments))
    } else if(overlapType == "within") {
      segord <- order(as.factor(seqnames(segments)), start(segments), end(segments))
    }
    
    coordord <- order(as.factor(seqnames(coordinates)), start(coordinates), end(coordinates))
      
    seg <- segments[segord,, drop = FALSE]
    coord <- coordinates[coordord,, drop = FALSE]

    chrOverlaps <- lapply(seqlevels(coordinates), function(chr) {
      whchr <- which(seqnames(coord) == chr)
      chrcoord <- coord[whchr,, drop = FALSE]
      whseg <- which(seqnames(seg) == chr)
      chrseg <- seg[whseg,,drop = FALSE]

      if(length(chrcoord) == 0) return()
      if(length(chrseg) == 0 & !whichOverlaps) return(rep(FALSE, length(chrcoord)))
      if(length(chrseg) == 0 & whichOverlaps) return(rep(NA, length(chrcoord)))

      if(overlapType == "overlapping") {
        fIns <- cbind(findInterval(start(chrcoord), cummax(end(chrseg))) + 1,
                      findInterval(end(chrcoord), cummax(start(chrseg))))
      } else if (overlapType == "contains") {
        fIns <- cbind(findInterval(start(chrcoord) - 0.5, cummax(start(chrseg))) + 1,
                      findInterval(end(chrcoord), end(chrseg)))
      } else if(overlapType == "within") {
        fIns <- cbind(findInterval(end(chrcoord) - 0.5, cummax(end(chrseg))) + 1,
                      findInterval(start(chrcoord), start(chrseg)))
      }
        
      if(!whichOverlaps & overlapType %in% c("overlapping", "contains"))
        return(chrOverlaps <- as.list(fIns[,2] >= fIns[,1]))          

      chrOverlaps <- list()
      chrOverlaps[1:length(chrcoord)] <- NA
      
      if(!whichOverlaps & overlapType == "within")
        {
          coordCheck <- fIns[,1] <= fIns[,2]
          chrOverlaps[!coordCheck] <- FALSE
          if(any(coordCheck))
            {
              rodseg <- chrseg[order(start(chrseg), -end(chrseg)),,drop = FALSE]
              chrOverlaps[which(coordCheck)[which(end(rodseg)[match(start(chrseg)[fIns[coordCheck,1]], start(rodseg))] >= end(chrcoord)[coordCheck])]] <- TRUE
            }
          coordCheck <- is.na(chrOverlaps)
          if(any(coordCheck))
            {
              rodseg <- chrseg[order(end(chrseg), start(chrseg)),, drop = FALSE]
              chrOverlaps[which(coordCheck)[which(start(rodseg)[match(end(chrseg)[fIns[coordCheck,2]], end(rodseg))] <= start(chrcoord)[coordCheck])]] <- TRUE
            }
          coordCheck <- which(is.na(chrOverlaps))
          if(length(coordCheck) == 0) return(chrOverlaps)
        } else if(whichOverlaps & overlapType == "within") {
          coordCheck <- which(fIns[,1] <= fIns[,2])
        } else if(overlapType %in% c("overlapping", "contains"))
        {
          chrOverlaps[fIns[,1] == fIns[,2]] <- segord[whseg[fIns[fIns[,1] == fIns[,2],1]]]
          coordCheck <- which(fIns[,1] < fIns[,2])
        }
        
      if(length(coordCheck) > 0)
        {
          checkOverlaps <- function(co, whichOverlaps = TRUE)
            {
              start <- co[1]
              end <- co[2]
              rfIStart <- co[3]
              rfIEnd <- co[4]
              
              selseg <- rseg[rfIStart:rfIEnd,,drop = FALSE]
              if(overlapType == "overlapping") {
                return(segord[rwhseg[(rfIStart:rfIEnd)[end(selseg) >= start & start(selseg) <= end]]])
              } else if(overlapType == "within") {
                whichWithin <- start(selseg) <= start & end(selseg) >= end
                if(whichOverlaps) {
                  if(any(whichWithin)) return(segord[rwhseg[(rfIStart:rfIEnd)[whichWithin]]]) else return(NA)
                } else {
                  if(any(whichWithin)) return(TRUE) else return(FALSE)
                }
              } else if(overlapType == "contains") {
                return(segord[rwhseg[(rfIStart:rfIEnd)[start(selseg) >= start & end(selseg) <= end]]])
              }
            }
          
          if(!is.null(cl))
            {
              clustAssign <- function(object, name)
                {
                  assign(name, object, envir = .GlobalEnv)
                  NULL
                }
              overlapsEnv <- new.env(parent = .GlobalEnv)
              environment(clustAssign) <- overlapsEnv
              environment(checkOverlaps) <- overlapsEnv
              clusterCall(cl, clustAssign, fIns, "fIns")
              clusterCall(cl, clustAssign, chrseg, "chrseg")
              clusterCall(cl, clustAssign, chrcoord, "chrcoord")
              clusterCall(cl, clustAssign, segord, "segord")
              clusterCall(cl, clustAssign, whseg, "whseg")
              clusterCall(cl, clustAssign, overlapType, "overlapType")
            }          

          chkCoord <- chrcoord[coordCheck,,drop = FALSE]
          fIns <- fIns[coordCheck,, drop = FALSE]
          
          numwin <- ceiling(length(chkCoord) / min(50000, length(chkCoord)))
          ends <- round(quantile(end(chkCoord), probs = 1:numwin / numwin))
          windows <- data.frame(start = c(1, ends[-length(ends)] + 1), end = ends)

          for(ii in 1:nrow(windows))
            {
              if(ii %% 10 == 0) message(".", appendLF = FALSE)
              wrcoord <- which(end(chkCoord) >= windows$start[ii] & end(chkCoord) <= windows$end[ii])
              rcoord <- chkCoord[wrcoord,, drop = FALSE]
              rfIns <- fIns[wrcoord,, drop = FALSE]
              adj <- min(rfIns) - 1

              if(!is.null(cl))
                {
                  clustRedSeg <- function(minseg, maxseg) {
                    library(GenomicRanges)
                    rwhseg <- whseg[minseg:maxseg]
                    rseg <- chrseg[minseg:maxseg,,drop = FALSE]
                    assign("rseg", rseg, envir = .GlobalEnv)
                    assign("rwhseg", rwhseg, envir = .GlobalEnv)
                    NULL
                  }
                  
                  environment(clustRedSeg) <- overlapsEnv

                  clusterCall(cl, clustRedSeg, min(rfIns), max(rfIns))

                } else{
                  rseg <- chrseg[min(rfIns):max(rfIns),, drop = FALSE]
                  rwhseg <- whseg[min(rfIns):max(rfIns)]
                }

              rfIns <- rfIns - adj

              if(!is.null(cl)) {
                apResult <- parApply(cl, cbind(start(rcoord), end(rcoord), rfIns), 1, checkOverlaps, whichOverlaps = whichOverlaps)
              } else apResult <- apply(cbind(start(rcoord), end(rcoord), rfIns), 1, checkOverlaps, whichOverlaps = whichOverlaps)

              if(whichOverlaps) {
                if(is.matrix(apResult)) apResult <- as.list(as.data.frame(apResult))
                chrOverlaps[coordCheck[wrcoord]] <- as.list(apResult)
              } else chrOverlaps[coordCheck[wrcoord]] <- apResult
            }
          
        }

      chrOverlaps
      
    })

    if(!is.null(cl))
      clusterEvalQ(cl, rm(list = ls()))

    overlaps <- list()
    for(chr in seqlevels(coord))
      {
        whchr <- which(seqnames(coord) == as.character(chr))
        if(length(whchr) > 0)
          overlaps[coordord[whchr]] <- chrOverlaps[[which(levels(seqnames(coord)) == chr)]]        
      }
    if(whichOverlaps)
      {
        overlaps <- lapply(overlaps, function(x) if(all(is.na(x))) return(integer(0)) else return(x))        
        return(overlaps)
      } else return(unlist(overlaps))
  }

