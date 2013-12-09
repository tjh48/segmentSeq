.averageMethylationRegions <- function (mD, modcod, redOvers, position, samples, coordinates, lenMod, cuts, surrounding = 0)
{
  redMD <- mD  
  uus <- unique(unlist(samples))
  redsamp <- lapply(samples, match, uus)

  modC <- redMD@Cs[unlist(redOvers), uus, drop = FALSE]/as.integer(redMD@alignments$multireads[unlist(redOvers)])
  modT <- redMD@Ts[unlist(redOvers), uus, drop = FALSE]/as.integer(redMD@alignments$multireads[unlist(redOvers)])
  lenOvers <- sapply(redOvers, length)
  startData <- cbind(start = unlist(redOvers), codstart = rep(start(coordinates), 
                                                 lenOvers), codwidth = rep(width(coordinates), lenOvers))
  startData[, 1] <- start(redMD@alignments)[startData[, 1]] - startData[, 2]
  minus <- (rep(as.character(strand(coordinates)), lenOvers) == "-")
  if (surrounding > 0) {
    rights <- (startData[, 1] >= startData[, 3])
    insides <- (startData[, 1] >= 0 & startData[, 1] < startData[, 
                           3])
    lefts <- (startData[, 1] < 0)
    startData[insides, 1] <- startData[insides, 1]/startData[insides, 
                                                             3] * surrounding
    startData[rights, 1] <- startData[rights, 1] - startData[rights, 
                                                             3] + surrounding
    startData[minus & insides, 1] <- surrounding - startData[minus & 
                                                             insides, 1]
    startData[minus & lefts] <- startData[minus & lefts] * 
      -1 + surrounding
    startData[minus & rights] <- (startData[minus & rights] - 
                                  surrounding) * -1
    } else {
    startData[, 1] <- startData[, 1]/startData[, 3]
    startData[minus] <- 1 - startData[minus]
    startData[, 1] <- startData[, 1] * lenMod
  }
  if (surrounding > 0) {
    splits <- cut(startData[, 1], breaks = (-cuts):(2 * cuts)/cuts * surrounding, include.lowest = TRUE, labels = FALSE)
  } else splits <- cut(startData[, 1], breaks = 0:cuts/cuts * lenMod, include.lowest = TRUE, labels = FALSE)
                       
  methProps <- lapply(redsamp, function(samp) {
    Cs <- modC[, samp, drop = FALSE]
    Ts <- modT[, samp, drop = FALSE]
    sampProps <- rowMeans(Cs/(Cs + Ts), na.rm = TRUE)
    sapply(split(sampProps, factor(splits, levels = 1:c(3 * 
                                             cuts, cuts)[as.numeric(surrounding == 0) + 1])), 
           mean, na.rm = TRUE)
  })

  message(".", appendLF = FALSE)
  return(methProps)
}

.plotProfile <- function(position, profiles, samples, col, surrounding, ylim, lenMod, add, ...)
  {
    lapply(1:length(profiles), function(ii) {
      if(ii == 1 & add == FALSE) {
        coverage <- profiles[[ii]]
        plot(x = position, y = coverage, type= "l", col = col[ii], ylim = ylim, ...)
      } else lines(x = position, y = profiles[[ii]], col = col[ii], ...)
    })
    if(!add) {
      if(surrounding > 0) abline(v = c(surrounding, surrounding + lenMod), col = "orange", lty = 3, lwd = 2)
      legend(x = "topright", lty = 1, col = col, legend = names(samples))
    }
    
  }

.subProfile <- function(coordinates, position, lenMod, mD, samples, cuts, surrounding = 0)
  {
    modcod <- coordinates
    if(surrounding > 0) {
      start(modcod) <- pmax(start(modcod) - surrounding, 1)
      end(modcod) <- end(modcod) + surrounding
    }
            
    methOvers <- findOverlaps(modcod, mD@alignments, select = "all")
    overMD <- sort(unique(unlist(methOvers@subjectHits)))
    redOvers <- split(match(methOvers@subjectHits, overMD), as.factor(methOvers@queryHits))
    redMD <- mD[overMD,]
    
    if(class(mD) == "alignmentMeth") {      
      return(.averageMethylationRegions(mD = redMD, modcod = modcod, redOvers = redOvers, position = position, samples = samples, coordinates = coordinates, lenMod = lenMod, cuts = cuts, surrounding = surrounding))
    } else {
      uus <- unique(unlist(samples))
      redsamp <- lapply(samples, match, uus)
      
      chrwidths <- sapply(seqlevels(mD@alignments), function(chr) max(c(0, end(modcod[seqnames(modcod) == chr]))))
      
      if(class(mD) == "alignmentData") {
        RKPM <- t(t(getCounts(modcod, redMD, cl = NULL) / width(modcod)) / mD@libsizes) * 1e9
        sampRKPM <- sapply(redsamp, function(samp) {
          rowMeans(RKPM[,samp,drop = FALSE])
        })
        minRKPM <- max(quantile(sampRKPM[sampRKPM != 0], 0.1), 10)
        
                                        #      RKPM <- t(t(getCounts(coordinates, redMD, cl = NULL) / width(coordinates)) / mD@libsizes) * 1e9
                                        #      end(leftMod[strand(leftMod) != "-"]) <- start(coordinates[strand(leftMod) != "-"]) - 1
                                        #      start(leftMod[strand(leftMod) == "-"]) <- end(coordinates[strand(leftMod) == "-"]) + 1
                                        #      leftRKPM <- t(t(getCounts(leftMod, redMD, cl = NULL) / width(leftMod)) / mD@libsizes) * 1e9
                                        #      start(rightMod[strand(rightMod) != "-"]) <- end(coordinates[strand(rightMod) != "-"]) + 1
                                        #      end(rightMod[strand(rightMod) == "-"]) <- start(coordinates[strand(rightMod) == "-"]) - 1
                                        #      rightRKPM <- t(t(getCounts(rightMod, redMD, cl = NULL) / width(rightMod)) / mD@libsizes) * 1e9
        
        
        codbase <- lapply(seqlevels(mD@alignments), function(chr) {
          message(".", appendLF = FALSE)
          bases <- as.integer(ranges(modcod[seqnames(modcod) == chr]))
          id <- rep(which(seqnames(modcod) == chr), width(modcod[seqnames(modcod) == chr]))
          adjPos <- bases - start(modcod)[as.integer(id)] + 1                
          centre <- Rle(bases >= start(coordinates)[id] & bases <= end(coordinates)[id])
          leftBound <- start(coordinates)[id] - start(modcod)[id]
          
          adjPos[which(centre)] <- (adjPos[which(centre)] - leftBound[which(centre)]) / (width(coordinates)[id])[which(centre)] * lenMod        
          minus <- Rle(as.character(strand(coordinates))[id] == "-")
          adjPos[which(centre & minus)] <- lenMod - adjPos[which(centre & minus)] + 1
          
          if(surrounding > 0) {
            left <- Rle(bases < start(coordinates)[id])
            right <-Rle(bases > end(coordinates)[id])
            adjPos[which(left)] <- adjPos[which(left)] + surrounding - leftBound[which(left)]
            adjPos[which(right)] <- adjPos[which(right)] - (width(coordinates)[id])[which(right)] - leftBound[which(right)]
            adjPos[which((left | right) & minus)] <- surrounding - adjPos[which((left | right) & minus)] + 1
            right[which((right | left) & minus)] <- !right[which((right | left) & minus)]
            left[which((right | left) & minus)] <- !left[which((right | left) & minus)]        
            return(DataFrame(bases = Rle(bases), adjPos = Rle(adjPos), id = Rle(id), left = (left), centre = (centre), right = (right)))
          } else return(DataFrame(bases = Rle(bases), adjPos = Rle(adjPos), id = Rle(id), centre = (centre)))
        })
        
        mergeBase <- do.call("rbind", codbase)[,-1]      
        splitBases <- function(wBase) if(any(wBase)) split(which(wBase), factor(cut(as.integer(mergeBase$adjPos[which(wBase)]), breaks = cuts, labels = FALSE, include.lowest = TRUE), levels = 1:cuts))
        leftSplit <- splitBases(mergeBase$left)
        rightSplit <- splitBases(mergeBase$right)
        centreSplit <- splitBases(mergeBase$centre)
        
        weightSplit <- function(split) {       
          numberWithin <- table(mergeBase$id[split])
          w = 1 / Rle(as.integer(numberWithin[match(as.integer(mergeBase$id[split]), names(numberWithin))]))
        }
        leftWeight <- lapply(leftSplit, weightSplit)
        rightWeight <- lapply(rightSplit, weightSplit)
        centreWeight <- lapply(centreSplit, weightSplit)
        
        overCod <- unique(unlist(redOvers))
        
        repDat <- lapply(redsamp, function(samp) {
          message(".", appendLF = FALSE)
          
          sampRKPM <- rowMeans(RKPM[,samp,drop = FALSE])
          sampRKPM[sampRKPM < minRKPM] <- minRKPM
          
          weights <- rowSums(t(t(redMD@data[overCod,samp] / redMD@alignments$matches[overCod]) / redMD@libsizes[samp]) * 1e6) / length(samp)
          covRep <- coverage(redMD@alignments[overCod[weights > 0]], weight = weights[weights > 0], width = chrwidths)
          
          nonzeros <- which(rowSums(redMD@data[,samp,drop = FALSE]) > 0)
          covNZ <- coverage(redMD@alignments[intersect(nonzeros, overCod)], width = chrwidths) > 0
          covRep <- covRep * covNZ
          
          coverage <- do.call("c",
                              lapply(1:length(covRep), function(chr) {
                                if(nrow(codbase[[chr]]) > 0) {            
                                  coverage <- covRep[[chr]][codbase[[chr]]$bases]
                                  coverage[coverage < 0] <- 0
                                  return(as.numeric(coverage))
                                } else return(c())
                              }))
          
                                        #        coverBase <- coverBase[which((coverBase$id %in% which(sampRKPM > 0) & coverBase$centre) |
                                        #                               (coverBase$id %in% which(sampLRKPM > 0) & coverBase$left) |
                                        #                               (coverBase$id %in% which(sampRRKPM > 0) & coverBase$right)),]
          
                                        #        if(any(coverBase$left)) leftSummary <- sapply(split(coverBase$coverage[coverBase$left], factor(cut(coverBase$adjPos[coverBase$left], breaks = cuts, labels = FALSE, include.lowest = TRUE), levels = 1:cuts)), mean, na.rm = TRUE) else leftSummary <- rep(0, cuts)
                                        #        if(any(coverBase$right)) rightSummary <- sapply(split(coverBase$coverage[coverBase$right], factor(cut(coverBase$adjPos[coverBase$right], breaks = cuts, labels = FALSE, include.lowest = TRUE), levels = 1:cuts)), mean, na.rm = TRUE) else rightSummary <- rep(0, cuts)
          
          summariseCoverage <- function(split, splitWeight, codRKPM) {
            sapply(1:length(split), function(ii) {            
              weighted.mean(coverage[split[[ii]]] / codRKPM[as.integer(mergeBase$id[split[[ii]]])], w = splitWeight[[ii]], na.rm = TRUE) * mean(codRKPM, trim = 0.1)
            })
          }
          
          if(surrounding > 0) {
            sumCov <- c(summariseCoverage(leftSplit, leftWeight, sampRKPM),
                        summariseCoverage(centreSplit, centreWeight, sampRKPM),
                        summariseCoverage(rightSplit, rightWeight, sampRKPM))
          } else sumCov <- summariseCoverage(centreSplit, centreWeight, sampRKPM)
          sumCov
        })
      }
            
      message(".", appendLF = FALSE)
      return(repDat)
    }
  }

averageProfiles <- function(mD, samples, coordinates, cuts, maxcuts = 200, bw = 5000, surrounding = 0, add = FALSE, col, ylim, ...)
  {
    message("Plotting...", appendLF = FALSE)
    
    if(missing(cuts)) cuts <- max(ceiling(median(width(coordinates)) / bw * length(coordinates)), 5)
#    cuts <- ceiling(cuts)
    cuts <- min(cuts, maxcuts)
    if(missing(samples)) samples <- mD@replicates
    if(is.factor(samples)) {
      sampNames <- levels(samples)
      samples <- lapply(levels(samples), function(rep) which(samples== rep))
      names(samples) <- sampNames
    }
    if(!is.list(samples)) samples <- list(samples)
    if(missing(col)) col <- rainbow(length(samples))

    if(length(unique(width(coordinates))) == 1) lenMod <- unique(width(coordinates)) else lenMod <- c(1000, surrounding)[as.integer(surrounding != 0) + 1]
    if(surrounding > 0) {
      position <- c(0.5 + 0:(cuts - 1) / cuts * surrounding,
                    surrounding + (0.5 + 0:(cuts - 1) / cuts * lenMod),
                    surrounding + lenMod + 0.5 + 0:(cuts - 1) / cuts * surrounding)
    } else position <- (0.5 + 0:(cuts - 1)) / cuts * lenMod      

    if(ceiling(length(coordinates) / 5000) > 1) {
      splitcod <- split(1:length(coordinates), cut(1:length(coordinates), breaks = ceiling(length(coordinates) / 5000), labels = FALSE))
    } else splitcod <- list(1:length(coordinates))
    
    splitProf <- lapply(splitcod, function(x) .subProfile(coordinates[x], position = position, lenMod = lenMod, mD = mD, samples = samples, cuts = cuts, surrounding = surrounding))
    profiles <- lapply(1:length(samples), function(ii) colMeans(do.call("rbind", lapply(splitProf, function(x) x[[ii]]))))
    message(".done!")
    if(is.null(ylim) || missing(ylim)) ylim <- c(0, max(sapply(profiles, max, na.rm = TRUE), na.rm = TRUE) * 1.1)
    .plotProfile(position = position, profiles = profiles, samples = samples, col = col, surrounding = surrounding, ylim = ylim, lenMod = lenMod,add = add, ...)
      
    invisible(list(position, profiles))
  }

plotMethDistribution <- function(meth, samples, bw = 1e-3, subtract, chrs, centromeres, add = FALSE, col, legend = TRUE, ...)
  {
    if(missing(samples)) samples <- meth@replicates
    if(missing(chrs)) chrs <- NULL
    
    if(is.factor(samples)) {
      namSamp <- levels(samples)
      samples <- lapply(levels(samples), function(rep) which(samples== rep))
      names(samples) <- namSamp
    }
    if(!is.list(samples)) samples <- list(samples) else samples <- samples
    if(missing(col)) col <- rainbow(length(samples))

    if(missing(subtract)) subtract <- NULL
    if(!is.list(subtract)) subtract <- lapply(1:length(samples), function(rep) subtract)
    
    if(missing(centromeres)) centromeres <- NULL
    
    methDist <- lapply(1:length(samples), function(samps)
                       .plotSampleMeth(meth[,samples[[samps]]], bw = bw, subtract = subtract[[samps]], col = col[samps], centromeres = centromeres, add = c(add, TRUE)[as.integer(samps > 1) + 1], chrs = chrs, ...))

    if(length(samples) > 1 & legend) legend(x = "topright", legend = names(samples), lwd = 2, lty = 1, col = col)
    invisible(methDist)
  }
    
.plotSampleMeth <- function(meth, bw = 1e-3, subtract, centromeres, chrs, add = FALSE, col = col, ...) #adjustByCs = FALSE, ...) {
  {
    chrlens <- seqlengths(meth@alignments)
    if(any(is.na(chrlens)))
      chrlens[is.na(chrlens)] <- sapply(names(chrlens)[is.na(chrlens)], function(nana) max(end(meth@alignments[seqnames(meth@alignments) == nana])))
    if(!is.null(chrs)) {
      chrlens <- chrlens[names(chrlens) %in% chrs]
    } else chrs <- names(chrlens)

    sumchr <- cumsum(c(0, chrlens))
    CHGpos <- do.call("c", lapply(1:length(chrs), function(chr) {
      (end(meth@alignments[seqnames(meth@alignments) == chrs[chr]]) + start(meth@alignments[seqnames(meth@alignments) == chrs[chr]])) / 2 + sumchr[chr]
    }))

    adjCs <- rowSums(meth@Cs / as.integer(meth@alignments$multireads))
    adjTs <- rowSums(meth@Ts / as.integer(meth@alignments$multireads))
    
    breaks <- do.call("c", lapply(1:length(chrlens), function(ii) {
      clen <- chrlens[ii]
      cuts <- round(clen * bw)
      acl <- clen / round(clen / cuts)
      1:(clen / acl) * acl + cumsum(c(0, chrlens))[ii]
    }))
    
                                        #  breaks <- do.call("c", lapply(1:length(seqlevels(meth@alignments)), function(ii) {
                                        #    chr <- seqlevels(meth@alignments)[ii]
                                        #    quantile(unique(CHGpos), probs = seq(0, 1, 0.1)) + c(0, cumsum(chrlens))[ii]
                                        #  }))
    
    
    ints <- findInterval(breaks, c(CHGpos), all.inside = TRUE)
    sumC <- cumsum(adjCs)[ints]; sumT <- cumsum(adjTs)[ints]
    countC <- sumC - c(0, sumC[-length(sumC)]); countT <- sumT - c(0, sumT[-length(sumT)])
    methdiv <- (countC / (countC + countT))
    
                                        #  if(adjustByCs) {
                                        #    siteNums <- findInterval(breaks, unique(CHGpos))
                                        #    Cadj <- siteNums - c(0, siteNums[-length(breaks)])
                                        #    methdiv <- methdiv * Cadj / max(Cadj)
                                        #  }
    
    if(!is.null(subtract)) {
      methdiv <- methdiv - subtract
      ylim <- c(-1.1, 1)
    } else ylim = c(-0.1, 1)

    if(missing(col)) col = "black"
    
    methylation <- methdiv[!is.na(methdiv)]
    position <- breaks[!is.na(methdiv)] #- cuts / 2
    if(!add) {
      plot(x = position, y = methylation, type = "l", axes = FALSE, ylim = ylim, col = col, ylab = "Proportion of methylation", ...)
      axis(2, at = pretty(0:1, n = 5))
      if(length(chrlens) > 1)
        segments(x0 = cumsum(chrlens)[-length(chrlens)], y0 = 0, y1 = 1, col = "red", lty = 2, lwd = 3)    
      if(!missing(centromeres) && !is.null(centromeres)) segments(x0 = c(centromeres[,1] + c(0, cumsum(chrlens)[-length(chrlens)]), centromeres[,2] + c(0, cumsum(chrlens)[-length(chrlens)])), y0 = 0, y1 = 1, lty = 2, lwd = 2, col = "blue")
      text(names(chrlens), srt = 10, adj = 1, y = ylim[1] + 0.05, x = cumsum(c(0, chrlens[-length(chrlens)])) + chrlens / 2, cex = 2.5)
    } else lines(x = position, y = methylation, col = col, ...)
    invisible(data.frame(position = breaks, methylation = methdiv))
  }   
