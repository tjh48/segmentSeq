getOverlaps <- function(coordinates, segments, overlapType = "overlapping", whichOverlaps = TRUE, cl)
  {
    if(overlapType == "overlapping") {
      segord <- order(as.factor(segments$chr), segments$start, segments$end)
    } else if (overlapType == "contains") {
      segord <- order(as.factor(segments$chr), segments$end, segments$start)
    } else if(overlapType == "within") {
      segord <- order(as.factor(segments$chr), segments$start, segments$end)
    }
    
    coordord <- order(as.factor(coordinates$chr), coordinates$start, coordinates$end)
      
    seg <- segments[segord,, drop = FALSE]
    coord <- coordinates[coordord,, drop = FALSE]

    chrOverlaps <- lapply(unique(coord$chr), function(chr) {
      whchr <- which(coord$chr == chr)
      chrcoord <- coord[whchr,, drop = FALSE]
      whseg <- which(seg$chr == chr)
      chrseg <- seg[whseg,,drop = FALSE]

      if(overlapType == "overlapping") {
        fIns <- cbind(findInterval(chrcoord$start, cummax(chrseg$end)) + 1,
                      findInterval(chrcoord$end, cummax(chrseg$start)))
      } else if (overlapType == "contains") {
        fIns <- cbind(findInterval(chrcoord$start - 0.5, cummax(chrseg$start)) + 1,
                      findInterval(chrcoord$end, chrseg$end))
      } else if(overlapType == "within") {
        fIns <- cbind(findInterval(chrcoord$end - 0.5, cummax(chrseg$end)) + 1,
                      findInterval(chrcoord$start, chrseg$start))
      }
        
      if(!whichOverlaps & overlapType %in% c("overlapping", "contains"))
        return(chrOverlaps <- as.list(fIns[,2] >= fIns[,1]))          

      chrOverlaps <- list()
      chrOverlaps[1:nrow(chrcoord)] <- NA
      
      if(!whichOverlaps & overlapType == "within")
        {
          coordCheck <- fIns[,1] <= fIns[,2]
          chrOverlaps[!coordCheck] <- FALSE
          if(any(coordCheck))
            {
              rodseg <- chrseg[with(chrseg, order(start, -end)),,drop = FALSE]
              chrOverlaps[which(coordCheck)[which(rodseg$end[match(chrseg$start[fIns[coordCheck,1]], rodseg$start)] >= chrcoord$end[coordCheck])]] <- TRUE
            }
          coordCheck <- is.na(chrOverlaps)
          if(any(coordCheck))
            {
              rodseg <- chrseg[with(chrseg, order(end, start)),, drop = FALSE]
              chrOverlaps[which(coordCheck)[which(rodseg$start[match(chrseg$end[fIns[coordCheck,2]], rodseg$end)] <= chrcoord$start[coordCheck])]] <- TRUE
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
                return(segord[rwhseg[(rfIStart:rfIEnd)[selseg$end >= start & selseg$start <= end]]])
              } else if(overlapType == "within") {
                whichWithin <- selseg$start <= start & selseg$end >= end
                if(whichOverlaps) {
                  if(any(whichWithin)) return(segord[rwhseg[(rfIStart:rfIEnd)[whichWithin]]]) else return(NA)
                } else {
                  if(any(whichWithin)) return(TRUE) else return(FALSE)
                }
              } else if(overlapType == "contains") {
                return(segord[rwhseg[(rfIStart:rfIEnd)[selseg$start >= start & selseg$end <= end]]])
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
          
          numwin <- ceiling(nrow(chkCoord) / min(5000, nrow(chkCoord)))
          ends <- round(quantile(chkCoord$end, probs = 1:numwin / numwin))
          windows <- data.frame(start = c(1, ends[-length(ends)] + 1), end = ends)
          
          for(ii in 1:nrow(windows))
            {
              if(ii %% 10 == 0) message(".", appendLF = FALSE)
              wrcoord <- which(chkCoord$end >= windows$start[ii] & chkCoord$end <= windows$end[ii])
              rcoord <- chkCoord[wrcoord,, drop = FALSE]
              rfIns <- fIns[wrcoord,, drop = FALSE]
              adj <- min(rfIns) - 1
              
              if(!is.null(cl))
                {
                  clustRedSeg <- function(minseg, maxseg) {
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
                apResult <- parApply(cl, cbind(rcoord$start, rcoord$end, rfIns), 1, checkOverlaps, whichOverlaps = whichOverlaps)
              } else apResult <- apply(cbind(rcoord$start, rcoord$end, rfIns), 1, checkOverlaps, whichOverlaps = whichOverlaps)

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
    for(chr in unique(coord$chr))
      {
        whchr <- which(coord$chr == as.character(chr))
        overlaps[coordord[whchr]] <- chrOverlaps[[which(unique(coord$chr) == chr)]]        
      }
    if(whichOverlaps)
      {
        overlaps <- lapply(overlaps, function(x) if(all(is.na(x))) return(integer(0)) else return(x))        
        return(overlaps)
      } else return(unlist(overlaps))
  }

