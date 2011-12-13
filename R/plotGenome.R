plotGenome <-
function(aD, locData, chr = 1, limits = c(0, 1e4), samples = NULL, plotType = "pileup", plotDuplicated = FALSE, density = 0, showNumber = TRUE, logScale = FALSE, cap = Inf, ...)
{
  if(!(chr %in% seqlevels(aD@alignments))) stop(paste("Chromosome", chr, "does not seem to exist within the alignment data of the 'aD' object."))
  
  reducedAD <- aD[seqnames(aD@alignments) == chr & start(aD@alignments) <= limits[2] & end(aD@alignments) >= limits[1],]

  seqlevels(reducedAD@alignments) <- chr
  
  libsizes = reducedAD@libsizes
  cTags <- reducedAD@alignments
  cdata <- reducedAD@data

  if(is.null(samples)) samples <- 1:ncol(reducedAD@data)

  libnames = reducedAD@libnames[samples]

  plot(NA, NA, xlim = limits, ylim = c(0L, length(samples) + 1L), xlab = "", ylab = "", axes = FALSE, ...)
  axis(side = 1, at = 0:10 * round(diff(limits) / 10)  + limits[1], ...)
  axis(side = 2, at = 1:length(libnames), labels = libnames, las = 1, ...)
  
#  if(plotDuplicated & plotType == "pileup") segments(limits[1], 1:length(samples), limits[2], 1:length(samples), col = "cyan", lwd = 1) else segments(limits[1], 1:length(samples) - 0.45, limits[2], 1:length(samples)- 0.45, col = "cyan", lwd = 1)    
  segments(limits[1], 1:length(samples) - 0.5, limits[2], 1:length(samples) - 0.5, col = "black", lwd = 1)
  segments(limits[1], 1:length(samples) + 0.5, limits[2], 1:length(samples) + 0.5, col = "black", lwd = 1)

  message("Computing plot...", appendLF = FALSE)
  
  pAD <- processAD(reducedAD, gap = 0, cl = NULL, verbose = FALSE)
  
  if(plotType == "pileup") {
    rectPileup <- do.call("rbind", lapply(1:length(samples), function(uu)
           {
             sTags <- cTags
             values(sTags)$count <- cdata[,samples[uu]]
             
             if(plotDuplicated)
               {                 
                 uTags <- sTags[values(sTags)$matches == 1,]
                 uTags <- uTags[order(start(uTags), end(uTags)),]
                 uTags <- uTags[values(uTags)$count > 0,, drop = FALSE]
                 
                 dTags <- sTags[values(sTags)$matches > 1,]
                 dTags <- dTags[order(start(dTags), end(dTags)),]
                 dTags <- dTags[values(dTags)$count > 0,, drop = FALSE]
                 
                 if(nrow(uTags) > 0)
                   {
                     cpu <- coverage(IRanges(start=rep(start(uTags), as.integer(values(uTags)$count)),
                                             end=rep(end(uTags), as.integer(values(uTags)$count))))
                     rectcpu <- cbind(start(cpu), uu, end(cpu), runValue(cpu))                     
                   } else rectcpu <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                 if(nrow(dTags) > 0)
                   {
                     cpd <- coverage(IRanges(start=rep(start(dTags), as.integer(values(dTags)$count)),
                                             end=rep(end(dTags), as.integer(values(dTags)$count))))
                     rectcpd <- cbind(start(cpd), uu, end(cpd), -runValue(cpd))
                   } else rectcpd <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                 coverageOverlaps <- rbind(rectcpu, rectcpd)                 
               } else {
                 tagOverlaps <- getOverlaps(pAD@coordinates, sTags, overlapType = "contains", cl = NULL)

                 coverageOverlaps <- do.call("rbind", lapply(1:nrow(pAD), function(ii) {
                   tags <- tagOverlaps[[ii]]
                   rectcpu <- matrix(ncol = 4, nrow = 0)
                   if(length(tags) > 0)
                     {
                       cpu <- coverage(IRanges(start=rep(start(sTags)[tags], as.integer(values(sTags)$count[tags])),
                                               end=rep(end(sTags)[tags], as.integer(values(sTags)$count[tags]))))
                       if(length(cpu) > 0)
                         {
                           rectcpu <- matrix(cbind(start(cpu), uu - 0.45, end(cpu), runValue(cpu)), ncol = 4)
                           rectcpu[,4] <- sapply(rectcpu[,4], min, cap)
                           rectcpu <- rectcpu[rectcpu[,4] > 0,,drop = FALSE]                           
                           rectcpu[,4] <- rectcpu[,4] / sum(values(sTags)$count[tags]) * as.integer(pAD@data[ii,samples[uu]])
                           
                         }
                     }
                   rectcpu
                 }))
               }
               
             coverageOverlaps[,4] <- coverageOverlaps[,4] / reducedAD@libsizes[samples[uu]]             
             return(coverageOverlaps)
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
    rectcpu <- sapply(1:ncol(pAD), function(ii) as.integer(pAD@data[,ii]) / width(pAD@coordinates) / pAD@libsizes[ii])

    if(logScale) {
      rectcpu <- rectcpu / min(rectcpu[rectcpu > 0])
      rectcpu <- log2(rectcpu)
    }

    maxscale <- max(rectcpu) / 0.9
    sapply(1:length(samples), function(uu) {
      samp <- samples[uu]
      if(any(rectcpu[,samp] > 0))
        rect(start(pAD@coordinates)[rectcpu[,samp] > 0], uu - 0.45, end(pAD@coordinates)[rectcpu[,samp] > 0], uu - 0.45 + (rectcpu[rectcpu[,samp] > 0, samp]) / maxscale, col = "black")
    })
  }

  message(".done!")

  if(!missing(locData))
    {
      breakpoints <- locData@coordinates
      brlims <- which(seqnames(breakpoints) == chr & start(breakpoints) <= limits[2] & end(breakpoints) >= limits[1])
      if(length(brlims) > 0)
        {
          bps <- breakpoints[brlims,]
          start(bps)[start(bps) < limits[1]] <- limits[1]
          end(bps)[end(bps) > limits[2]] <- limits[2]

          if(!is.null(locData@locLikelihoods) & all(dim(locData@locLikelihoods) != 0)) {
            sapply(1:length(samples), function(ss) {
              alpha = (exp(locData@locLikelihoods[brlims,locData@replicates[samples[ss]]]))
              alpha[is.na(alpha)] <- 0
              if(any(alpha > 0.1))
                {
                  brcols <- rgb((rep(c(1,0,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                                (rep(c(0,0,1), ceiling(length(brlims) / 3)))[1:length(brlims)],
                                (rep(c(0,1,0), ceiling(length(brlims) / 3)))[1:length(brlims)], alpha = alpha)
                  angle <- rep(c(45, 115, 165), ceiling(length(brlims) / 3))
                  rect(start(bps)[alpha > 0.1], -.5 + ss, end(bps)[alpha>0.1], ss + .5, density = density, col = brcols[alpha > 0.1], border = brcols[alpha > 0.1], angle = angle[alpha > 0.1])
                  segments(start(bps)[alpha > 0.1], -.5 + ss, end(bps)[alpha > 0.1], ss + .5, col = brcols[alpha > 0.1])
                }
            })
          } else {
            brcols <- rgb((rep(c(1,0,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                          (rep(c(0,1,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                          (rep(c(0,0,1), ceiling(length(brlims) / 3)))[1:length(brlims)], alpha = 1)
            rect(start(bps), 0 + .5, end(bps), length(samples) + 0.5, density = density, col = brcols, border = brcols, angle = 45)            
          }

          if(showNumber) text((start(bps) + end(bps)) / 2, y = 0, labels = brlims, col = "black")            
        }
    }

  
  invisible(NULL)
}

