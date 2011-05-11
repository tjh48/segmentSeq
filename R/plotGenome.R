plotGenome <-
function(aD, sD, chr = 1, limits = c(0, 1e4), samples = NULL, plotType = "pileup", plotDuplicated = FALSE, density = 0, showNumber = TRUE, logScale = FALSE, ...)
{
  if(!(chr %in% aD@chrs$chr)) stop(paste("Chromosome", chr, "is not defined in the 'aD@chrs$chr' slot."))
  
  reducedAD <- aD[aD@alignments$chr == chr & aD@alignments$start <= limits[2] & aD@alignments$end >= limits[1],]
  reducedAD@chrs <- reducedAD@chrs[reducedAD@chrs$chr == chr,]
  
  libsizes = reducedAD@libsizes
  cTags <- reducedAD@alignments
  cdata <- reducedAD@data

  if(is.null(samples)) samples <- 1:ncol(reducedAD@data)

  libnames = reducedAD@libnames[samples]

  plot(NA, NA, xlim = limits, ylim = c(0L, length(samples) + 1L), xlab = "", ylab = "", axes = FALSE,...)
  axis(side = 1, at = 0:10 * round(diff(limits) / 10)  + limits[1])
  axis(side = 2, at = 1:length(libnames), labels = libnames, las = 1)
  
  if(plotDuplicated & plotType == "pileup") segments(limits[1], samples, limits[2], samples, col = "cyan", lwd = 1) else segments(limits[1], samples - 0.45, limits[2], samples - 0.45, col = "cyan", lwd = 1)    
  segments(limits[1], samples - 0.5, limits[2], samples - 0.5, col = "black", lwd = 1)
  segments(limits[1], samples + 0.5, limits[2], samples + 0.5, col = "black", lwd = 1)

  message("Computing plot...", appendLF = FALSE)
  
  pAD <- processAD(reducedAD, gap = 0, cl = NULL, verbose = FALSE)
  
  if(plotType == "pileup") {
    rectPileup <- do.call("rbind", lapply(1:length(samples), function(uu)
           {
             sTags <- cbind(subset(cTags, select = c("chr", "start", "end", "matches")), count = cdata[,samples[uu]])
             if(plotDuplicated)
               {
                 
                 uTags <- sTags[sTags$matches == 1,]
                 uTags <- uTags[order(uTags$start, uTags$end),]
                 uTags <- uTags[uTags$count > 0,, drop = FALSE]
                 
                 dTags <- sTags[sTags$matches > 1,]
                 dTags <- dTags[order(dTags$start, dTags$end),]
                 dTags <- dTags[dTags$count > 0,, drop = FALSE]
                 
                 if(nrow(uTags) > 0)
                   {
                     cpu <- coverage(IRanges(start=rep(uTags$start, uTags$count), end=rep(uTags$end, uTags$count)))
                     rectcpu <- cbind(start(cpu), uu, end(cpu), runValue(cpu))                     
                   } else rectcpu <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                 if(nrow(dTags) > 0)
                   {
                     cpd <- coverage(IRanges(start=rep(dTags$start, dTags$count), end=rep(dTags$end, dTags$count)))
                     rectcpd <- cbind(start(cpd), uu, end(cpd), -runValue(cpd))
                   } else rectcpd <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                 rectcpu <- rbind(rectcpu, rectcpd)
                                        #                 rect(rectcpu[rectcpu[,4] > 0,1], uu, rectcpu[rectcpu[,4] > 0,3], uu + (rectcpu[rectcpu[,4] > 0,4]) / maxscale, col = "black")
                                        #                 rect(rectcpd[rectcpd[,4] > 0,1], uu, rectcpd[rectcpd[,4] > 0,3], uu - (rectcpd[rectcpd[,4] > 0,4]) / maxscale, col = "black")
               } else {
                 tagOverlaps <- getOverlaps(pAD@segInfo, sTags, overlapType = "contains", cl = NULL)

                 coverageOverlaps <- do.call("rbind", lapply(1:nrow(pAD), function(ii) {
                   tags <- tagOverlaps[[ii]]
                   rectcpu <- matrix(ncol = 4, nrow = 0)
                   if(length(tags) > 0)
                     {
                       cpu <- coverage(IRanges(start=rep(sTags$start[tags], sTags$count[tags]), end=rep(sTags$end[tags], sTags$count[tags])))
                       if(length(cpu) > 0)
                         {
                           rectcpu <- matrix(cbind(start(cpu), uu - 0.45, end(cpu), runValue(cpu)), ncol = 4)
                           rectcpu <- rectcpu[rectcpu[,4] > 0,,drop = FALSE]
                           rectcpu[,4] <- rectcpu[,4] / sum(sTags$count[tags]) * pAD@data[ii,samples[uu]]
                         }
                     }
                   rectcpu
                 }))
                   
#                 scaleChunk <- sapply(getOverlaps(pAD@segInfo, data.frame(chr = pAD@segInfo$chr[1], start = coverageOverlaps[,1], end = coverageOverlaps[,3]), overlapType = "contains", cl = NULL), function(x) if(length(x) == 0) return(0) else return(sum(coverageOverlaps[x,4] * (coverageOverlaps[x,3] - coverageOverlaps[x,1] + 1))))
                 
#                 rectChunk <- unlist(getOverlaps(
#                                                 data.frame(chr = pAD@segInfo$chr[1], start = coverageOverlaps[,1], end = coverageOverlaps[,3])
#                                                 , pAD@segInfo, overlapType = "within", cl = NULL))
                 
#                 coverageOverlaps[,4] <- coverageOverlaps[,4] / scaleChunk[rectChunk] * pAD@data[rectChunk,samples[uu]] * (pAD@segInfo$end[rectChunk] - pAD@segInfo$start[rectChunk] + 1)
                 
                 
               }
               
             coverageOverlaps[,4] <- coverageOverlaps[,4] / reducedAD@libsizes[samples[uu]]
             
#             message(".", appendLF = FALSE)
             coverageOverlaps
           }))

    rectPileup[,4] <- rectPileup[,4] / min(rectPileup[,4])

    if(logScale)
      {
        rectPileup[rectPileup[,4] > 0,4] <- log2(rectPileup[rectPileup[,4] > 0, 4])
        rectPileup[rectPileup[,4] < 0,4] <- -log2(-rectPileup[rectPileup[,4] < 0, 4])
      }

    rectPileup[,4] <- rectPileup[,4] / max(rectPileup[,4])

    rect(rectPileup[,1], rectPileup[,2], rectPileup[,3], rectPileup[,2] + rectPileup[,4], col = "black", border = "black")
    
  } else if(plotType == "chunk") {
    rectcpu <- t(t(pAD@data / (pAD@segInfo$end - pAD@segInfo$start + 1)) / pAD@libsizes)
    if(logScale) {
      rectcpu <- rectcpu / min(rectcpu[rectcpu > 0])
      rectcpu <- log2(rectcpu)
    }

    maxscale <- max(rectcpu) / 0.9
    sapply(1:length(samples), function(uu) {
      samp <- samples[uu]
      if(any(rectcpu[,samp] > 0))
        rect(pAD@segInfo$start[rectcpu[,samp] > 0], uu - 0.45, pAD@segInfo$end[rectcpu[,samp] > 0], uu - 0.45 + (rectcpu[rectcpu[,samp] > 0, samp]) / maxscale, col = "black")
    })
  }

  message(".done!")

  if(!missing(sD))
    {
      breakpoints <- sD@annotation
      brlims <- which(breakpoints$chr == chr & breakpoints$start <= limits[2] & breakpoints$end >= limits[1])
      if(length(brlims) > 0)
        {
          bps <- breakpoints[brlims,]
          bps$start <- sapply(bps$start, max, limits[1])
          bps$end <- sapply(bps$end, min, limits[2])

          if(!is.null(sD@posteriors) & all(dim(sD@posteriors) != 0)) {
            sapply(1:ncol(sD), function(ss) {
              alpha = (exp(sD@posteriors[brlims,sD@replicates[ss]]))
              alpha[is.na(alpha)] <- 0
              if(any(alpha > 0.1))
                {
                  brcols <- rgb((rep(c(1,0,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                                (rep(c(0,1,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                                (rep(c(0,0,1), ceiling(length(brlims) / 3)))[1:length(brlims)], alpha = alpha)
                  angle <- rep(c(45, 115, 165), ceiling(length(brlims) / 3))
                  rect(bps[alpha > 0.1,2], -.5 + ss, bps[alpha>0.1,3], ss + .5, density = density, col = brcols[alpha > 0.1], border = brcols[alpha > 0.1], angle = angle[alpha > 0.1])
                  segments(bps[,2], -.5 + ss, bps[,3], ss + .5, col = brcols)
                }
            })
          } else {
            brcols <- rgb((rep(c(1,0,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                          (rep(c(0,1,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                          (rep(c(0,0,1), ceiling(length(brlims) / 3)))[1:length(brlims)], alpha = 1)
            rect(bps[,2], 0 + .5, bps[,3], length(samples) + 0.5, density = density, col = brcols, border = brcols, angle = 45)
            segments(bps[,2], -.5 + ss, bps[,3], ss + .5, col = brcols)
          }

          if(showNumber) text((bps$start + bps$end) / 2, y = 0, labels = brlims, col = "black")            
        }
    }

  
  invisible(NULL)
}

