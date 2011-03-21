plotGenome <-
function(aD, sD, chr = 1, limits = c(0, 1e4), samples = NULL, plotType = "pileup", plotDuplicated = FALSE, ...)
{  
  libsizes = aD@libsizes
  cTags <- aD@alignments
  cdata <- aD@data
  

  if(!(chr %in% aD@chrs$chr)) stop(paste("Chromosome", chr, "is not defined in the 'aD@chrs$chr' slot."))
  
  samples <- 1:ncol(aD@data)

  libnames = aD@libnames[samples]

  plot(NA, NA, xlim = limits, ylim = c(0L, length(samples) + 1L), xlab = "", ylab = "", axes = FALSE, ...)
  axis(side = 1, at = 0:10 * round(diff(limits) / 10)  + limits[1])
  axis(side = 2, at = 1:length(libnames), labels = libnames, las = 1)
  
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
              alpha = (exp(sD@posteriors[brlims,sD@replicates[ss]]))^4
              alpha[is.na(alpha)] <- 0
              brcols <- rgb((rep(c(1,0,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                            (rep(c(0,1,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                            (rep(c(0,0,1), ceiling(length(brlims) / 3)))[1:length(brlims)], alpha = alpha)
              rect(bps[,2], -.5 + ss, bps[,3], ss + .5, density = 2, brcols, angle = 0:5 * 180 / 7)
            })
          } else {
            brcols <- rgb((rep(c(1,0,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                          (rep(c(0,1,0), ceiling(length(brlims) / 3)))[1:length(brlims)],
                          (rep(c(0,0,1), ceiling(length(brlims) / 3)))[1:length(brlims)], alpha = 1)
            rect(bps[,2], 0 + .5, bps[,3], length(samples) + 0.5, density = 2, brcols, angle = 0:5 * 180/7)
          }
          
          text((bps$start + bps$end) / 2, y = 0, labels = brlims, col = "black")
        }
    }

  if(plotDuplicated) segments(limits[1], samples, limits[2], samples, col = "cyan", lwd = 1) else segments(limits[1], samples - 0.45, limits[2], samples - 0.45, col = "cyan", lwd = 1)
    
  segments(limits[1], samples - 0.5, limits[2], samples - 0.5, col = "black", lwd = 1)
  segments(limits[1], samples + 0.5, limits[2], samples + 0.5, col = "black", lwd = 1)

  selTags <- which(cTags$chr == chr & cTags$start <= limits[2] & cTags$end >= limits[1])

  sapply(1:length(samples), function(uu)
         {
           sTags <- cbind(cTags, count = cdata[,uu])[selTags,]
           
           if(plotDuplicated)
             {
               uTags <- sTags[sTags$matches <= 1,]
               uTags <- uTags[order(uTags$start, uTags$end),]
               uTags <- uTags[uTags$count > 0,, drop = FALSE]
               
               dTags <- sTags[sTags$matches > 1,]
               dTags <- dTags[order(dTags$start, dTags$end),]
               dTags <- dTags[dTags$count > 0,, drop = FALSE]
               
               if(plotType == "pileup")
                 {
                   if(nrow(uTags) > 0)
                     {
                       cpu <- coverage(IRanges(start=rep(uTags$start, uTags$count), end=rep(uTags$end, uTags$count)))
                       rectcpu <- cbind(start(cpu), 0, end(cpu), runValue(cpu))
                     } else rectcpu <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                   if(nrow(dTags) > 0)
                     {
                       cpd <- coverage(IRanges(start=rep(dTags$start, dTags$count), end=rep(dTags$end, dTags$count)))
                       rectcpd <- cbind(start(cpd), 0, end(cpd), runValue(cpd))
                     } else rectcpd <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                   maxscale <- max(c(log(rectcpu[,4]), log(rectcpd[,4]), log(2))) * 2
                   
                   rect(rectcpu[rectcpu[,4] > 0,1], uu, rectcpu[rectcpu[,4] > 0,3], uu + log(rectcpu[rectcpu[,4] > 0,4]) / maxscale, col = "black")
                   rect(rectcpd[rectcpd[,4] > 0,1], uu, rectcpd[rectcpd[,4] > 0,3], uu - log(rectcpd[rectcpd[,4] > 0,4]) / maxscale, col = "black")
                 }
             } else {
               cpu <- coverage(IRanges(start=rep(sTags$start, sTags$count), end=rep(sTags$end, sTags$count)))
               if(length(cpu) > 0)
                 {
                   rectcpu <- cbind(start(cpu), 0, end(cpu), runValue(cpu))
                   maxscale <- max(log(rectcpu[,4]), log(2)) / 0.9
                   rect(rectcpu[rectcpu[,4] > 0,1], uu - 0.45, rectcpu[rectcpu[,4] > 0,3], uu - 0.45 + log(rectcpu[rectcpu[,4] > 0,4]) / maxscale, col = "black")
                 }
             }
           NULL
         })
  invisible(NULL)
}

