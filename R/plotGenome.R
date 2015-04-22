plotGenome <-
function(aD, loci, chr = 1, limits = c(0, 1e4), samples = NULL, plotType = "pileup", plotDuplicated = FALSE, density = 0, showNumber = TRUE, logScale = FALSE, cap = Inf, ...)
{
  if(!(chr %in% seqlevels(aD@alignments))) stop(paste("Chromosome", chr, "does not seem to exist within the alignment data of the 'aD' object."))
  
  reducedAD <- aD[which(seqnames(aD@alignments) == chr & start(aD@alignments) <= limits[2] & end(aD@alignments) >= limits[1]),]

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
  
  pAD <- processAD(reducedAD, gap = 0, cl = NULL, verbose = FALSE, getCounts = TRUE)
  
  if(plotType == "pileup") {
    rectPileup <- do.call("rbind", lapply(1:length(samples), function(uu)
           {
             sTags <- cTags
             values(sTags)$count <- cdata[,samples[uu]]
             
             if(plotDuplicated)
               {                 
                 uTags <- sTags[values(sTags)$multireads == 1,]
                 uTags <- uTags[order(start(uTags), end(uTags)),]
                 uTags <- uTags[values(uTags)$count > 0,, drop = FALSE]
                 
                 dTags <- sTags[values(sTags)$multireads > 1,]
                 dTags <- dTags[order(start(dTags), end(dTags)),]
                 dTags <- dTags[values(dTags)$count > 0,, drop = FALSE]
                 
                 if(length(uTags) > 0)
                   {
                     cpu <- coverage(IRanges(start=rep(start(uTags), as.integer(values(uTags)$count)),
                                             end=rep(end(uTags), as.integer(values(uTags)$count))))
                     rectcpu <- cbind(start(cpu), uu, end(cpu), runValue(cpu))                     
                   } else rectcpu <- matrix(c(NA, NA, NA, 1), nrow = 1, ncol = 4)
                 if(length(dTags) > 0)
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
                           if(any(rectcpu[,4] <= 0)) {
                             warning("There are non-zero count tags whose summed count is zero - something's gone wrong with the input data.")
                             rectcpu <- rectcpu[rectcpu[,4] > 0,,drop = FALSE]
                           }
                         }
                     }
                   if(any(rectcpu[,4] == 0)) message(ii)
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

    rect(rectPileup[,1] - 0.5, rectPileup[,2], rectPileup[,3] + 0.5, rectPileup[,2] + rectPileup[,4], col = "black", border = "black")
    
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
        rect(start(pAD@coordinates)[rectcpu[,samp] > 0] - 0.5, uu - 0.45, end(pAD@coordinates)[rectcpu[,samp] > 0] + 0.5, uu - 0.45 + (rectcpu[rectcpu[,samp] > 0, samp]) / maxscale, col = "black")
    })
  }

  message(".done!")

  if(!missing(loci))
    {
      breakpoints <- loci@coordinates
      brlims <- which(seqnames(breakpoints) == chr & start(breakpoints) <= limits[2] & end(breakpoints) >= limits[1])
      if(length(brlims) > 0)
        {
          bps <- breakpoints[brlims,]
          start(bps)[start(bps) < limits[1]] <- limits[1]
          end(bps)[end(bps) > limits[2]] <- limits[2]

          if(!is.null(loci@locLikelihoods) & all(dim(loci@locLikelihoods) != 0)) {
            sapply(1:length(samples), function(ss) {
              alpha = (exp(loci@locLikelihoods[brlims,loci@replicates[samples[ss]]]))
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



plotMeth <- function(aM, loci, chr, limits, samples, showNumber = TRUE, rgb = c(1,0,0), angle = 45, cap, add = FALSE)
  {
    normalise = FALSE
    if(missing(samples)) samples <- 1:ncol(aM)
    
    redADPlus <- aM[which(seqnames(aM@alignments) == chr & start(aM@alignments) >= limits[1] & end(aM@alignments) <= limits[2] & strand(aM@alignments) == "+"),]
    redADMinus <- aM[which(seqnames(aM@alignments) == chr & start(aM@alignments) >= limits[1] & end(aM@alignments) <= limits[2] & strand(aM@alignments) == "-"),]
    
    plusCs <- redADPlus@Cs / as.integer(redADPlus@alignments$multireads)
    plusTs <- redADPlus@Ts / as.integer(redADPlus@alignments$multireads)
    
    minusCs <- redADMinus@Cs / as.integer(redADMinus@alignments$multireads)
    minusTs <- redADMinus@Ts / as.integer(redADMinus@alignments$multireads)
        
    if(!missing(cap))
     {
       plusTotal <- plusCs + plusTs
       plusCs[plusTotal > cap] <- (plusCs / (plusTotal) * cap)[plusTotal > cap]
       plusTs[plusTotal > cap] <- (plusTs / (plusTotal) * cap)[plusTotal > cap]
       
       minusTotal <- minusCs + minusTs
       minusCs[minusTotal > cap] <- (minusCs / (minusTotal) * cap)[minusTotal > cap]
       minusTs[minusTotal > cap] <- (minusTs / (minusTotal) * cap)[minusTotal > cap]
     }
   
    maxVal <- 2.5 * max(plusCs + plusTs, minusCs + minusTs)

    par(mar = c(2, 5, 6, 2))
    if(!add) plot(NA, NA, xlim = limits, ylim = c(!showNumber - 0.5, length(samples) + 0.5), axes = FALSE, xlab = "Position", ylab = "")

    if(!missing(loci))
      {
        brlims <- which(seqnames(loci@coordinates) == chr & end(loci@coordinates) >= limits[1] & start(loci@coordinates) <= limits[2])
        
        redloci <- loci[brlims,]
                                        #       if(nrow(redloci@locLikelihoods) == 0)
                                        #         redloci@locLikelihoods <- matrix(0, ncol = length(redloci@replicates), nrow = nrow(redloci))
        
        if(nrow(redloci) > 0) {
          plotLoci <- TRUE
          xtext <- (start(redloci@coordinates) + end(redloci@coordinates)) / 2          
          dup <- duplicated(round(xtext / diff(limits), 1))          
          ytext <- unlist(apply(cbind(which(!dup), c(which(!dup)[-1] - 1, length(dup))), 1, function(x) x[1]:x[2] - x[1] + 1)) %% 5 * 0.2
          if(showNumber) text(xtext, ytext - 1, labels = brlims, col = "black", cex = 1)
        } else plotLoci <- FALSE         
      } else plotLoci <- FALSE

    sapply(samples, function(sN) {
      message("Plotting sample: ", sN)
      abline(h = which(sN == samples))
      
      if(plotLoci) {
        
        whrep <- which(levels(redloci@replicates) == redloci@replicates[which(samples == sN)])

        plotSelLoci <- function(selloc, loweradj, upperadj)
          {        
            if(nrow(selloc) > 0)
              {
                alpha = exp(selloc@locLikelihoods[,whrep]) * 0.75
                alpha[is.na(alpha)] <- 0            
                brcols <- rgb(rgb[1] * 0.7, rgb[2] * 0.7, rgb[3] * 0.7, alpha = alpha^2)            
                rect(start(selloc@coordinates) - 0.5, which(samples == sN) - loweradj, end(selloc@coordinates) + 0.5, which(samples == sN) + upperadj, col = brcols, border = brcols, density = 2, angle = angle)
              }
          }
        
      
        plotSelLoci(redloci[as.vector(redloci@locLikelihoods[,whrep] > -Inf & strand(redloci@coordinates) == "+"),], 0, 0.45)
        plotSelLoci(redloci[as.vector(redloci@locLikelihoods[,whrep] > -Inf & strand(redloci@coordinates) == "-"),], 0.45, 0)
        plotSelLoci(redloci[as.vector(redloci@locLikelihoods[,whrep] > -Inf & strand(redloci@coordinates) == "*"),], 0.45, 0.45)

      }
     
      if(!normalise)
       {         
         if(nrow(redADPlus) > 0) {
           rect(xleft = start(redADPlus@alignments) - 0.5, xright = end(redADPlus@alignments) + 0.5,
                ybottom = which(samples == sN) + (plusCs[,samples == sN]) / maxVal,
                ytop = which(samples == sN) + (plusCs[,samples == sN] + plusTs[,samples == sN]) / maxVal,
                col = "black", border = "black")
           
           rect(xleft = start(redADPlus@alignments) - 0.5, xright = end(redADPlus@alignments) + 0.5,
                ybottom = which(samples == sN),
                ytop = which(samples == sN) + (plusCs[,samples == sN]) / maxVal,
                col = rgb(rgb[1], rgb[2], rgb[3]), border = rgb(rgb[1], rgb[2], rgb[3]))
         }
         
         if(nrow(redADMinus) > 0) {
           rect(xleft = start(redADMinus@alignments) - 0.5, xright = end(redADMinus@alignments) + 0.5,
                ybottom = which(samples == sN) - (minusCs[,samples == sN]) / maxVal,
                ytop = which(samples == sN) - (minusCs[,samples == sN] + minusTs[,samples == sN]) / maxVal,
                col = "black", border = "black")
           
           rect(xleft = start(redADMinus@alignments) - 0.5, xright = end(redADMinus@alignments) + 0.5,
                ybottom = which(samples == sN),
                ytop = which(samples == sN) - (minusCs[,samples == sN]) / maxVal,
                col = rgb(rgb[1], rgb[2], rgb[3]), border = rgb(rgb[1], rgb[2], rgb[3]))
         }
         
       } #else if(normalise) {
                                        #rect(xleft = start(redADPlus@alignments), xright = end(redADPlus@alignments), ybottom = which(samples == sN), ytop = which(samples == sN) + (redADPlus@Cs[,ii]) / ((redADPlus@Cs[,ii]) + (redADPlus@Ts[,ii])) / 2.5, col = "red", border = "red")
          #rect(xleft = start(redADPlus@alignments), xright = end(redADPlus@alignments), ybottom = which(samples == sN), ytop = which(samples == sN) - (redADPlus@Cs[,ii]) / ((redADPlus@Cs[,ii]) + (redADPlus@Ts[,ii])) / 2.5, col = "red", border = "red")
        #}                 
    })
    
   axis(1)
    axis(2, labels = aM@libnames[samples], at = 1:length(samples), las = 2)
   
  }

