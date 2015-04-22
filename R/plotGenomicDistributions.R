.averageMethylationRegions <- function (redMD, subcoord, redOvers, position, samples, lenMod, cuts, surrounding)
{
  uus <- unique(unlist(samples))
#  redsamp <- lapply(samples, match, uus)

  nonconversion <- redMD@nonconversion
  modC <- redMD@Cs[unlist(redOvers), uus, drop = FALSE]/as.integer(redMD@alignments$multireads[unlist(redOvers)])
  modT <- redMD@Ts[unlist(redOvers), uus, drop = FALSE]/as.integer(redMD@alignments$multireads[unlist(redOvers)])

  ks <- t(t(modT) * nonconversion / (1 - nonconversion))
  modC <- modC - ks
  modT <- modT + ks
  modC[modC < 0] <- 0    
  
  lenOvers <- sapply(redOvers, length)
  startData <- cbind(start = unlist(redOvers),
                     codstart = rep(start(subcoord), lenOvers),
                     codwidth = rep(width(subcoord), lenOvers))                                                 
  startData[, 1] <- start(redMD@alignments)[startData[, 1]] - startData[, 2]
  minus <- (rep(as.character(strand(subcoord)), lenOvers) == "-")
  if (surrounding > 0) {
    rights <- (startData[, 1] >= startData[, 3])
    insides <- (startData[, 1] >= 0 & startData[, 1] < startData[, 3])                           
    lefts <- (startData[, 1] < 0)
    startData[insides, 1] <- startData[insides, 1]/startData[insides, 3] * surrounding
    startData[rights, 1] <- startData[rights, 1] - startData[rights, 3] + surrounding
    startData[minus & insides, 1] <- surrounding - startData[minus & insides, 1]
    startData[minus & lefts] <- startData[minus & lefts] * -1 + surrounding
    startData[minus & rights] <- (startData[minus & rights] - surrounding) * -1                                  
  } else {
    startData[, 1] <- startData[, 1]/startData[, 3]
    startData[minus] <- 1 - startData[minus]
    startData[, 1] <- startData[, 1] * lenMod
  }        
  if (surrounding > 0) {
    splits <- cut(startData[, 1], breaks = (-cuts):(2 * cuts)/cuts * surrounding, include.lowest = TRUE, labels = FALSE)
  } else splits <- cut(startData[, 1], breaks = 0:cuts/cuts * lenMod, include.lowest = TRUE, labels = FALSE)
  splits <- factor(splits, levels = 1:c(3 * cuts, cuts)[as.numeric(surrounding == 0) + 1])

  
  splitcod <- split(rep(1:length(redOvers), sapply(redOvers, length)), splits)
  splitdat <- split(1:nrow(modC), splits)
  sumstats <- mapply(function(x, y) {
    geneProps <- do.call("rbind", lapply(split(x, y), function(sp) {
      sumC <- colSums(modC[sp,,drop = FALSE])
      sumT <- colSums(modT[sp,,drop = FALSE])
      sumC / (sumC + sumT)
    }))
    sgProps <- lapply(samples, function(samps) geneProps[,samps,drop = FALSE])
    sumstat <- lapply(sgProps, function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), quantile(x, 1:19 / 20, na.rm = TRUE)))
    sumstat
  }, splitdat, splitcod, SIMPLIFY = FALSE)

  sampstats <- lapply(names(samples), function(samnam) do.call("rbind", lapply(sumstats, function(x) x[[samnam]])))

  message(".", appendLF = FALSE)
  return(sampstats)
}

plotAverageProfile <- function(position, profiles, col, surrounding, ylim, add = FALSE, meanOnly = TRUE, legend = TRUE, titles, ...)
  {
    if(meanOnly) {
      if(missing(col)) col <- rainbow(length(profiles))
      if(is.null(col)) col <- rainbow(length(profiles))
      lapply(1:length(profiles), function(ii) {
        coverage <- profiles[[ii]]
        if(is.matrix(coverage)) coverage <- coverage[,1]
        if(ii == 1 & !add) {
          plot(x = position, y = coverage, type= "l", col = col[ii], ylim = ylim, xlab = "position", ylab = "", ...)
        } else lines(x = position, y = coverage, col = col[ii], ...)
      })
      if(!add) {
        if(surrounding > 0) abline(v = c(surrounding, max(position) + diff(position[1:2]) - 0.5 - surrounding), col = "dark grey", lty = 3, lwd = 2)
        if("cex" %in% names(list(...))) legcex <- list(...)$cex else legcex <- 1
        if("lwd" %in% names(list(...))) leglwd <- list(...)$lwd else leglwd <- 1
        if(legend) legend(x = "topleft", lty = 1, lwd = leglwd, col = col, legend = names(profiles), cex = legcex, bty = "n")
      }
    } else {      
      if(names(dev.cur()) == "null device" || all(par()$mfrow == c(1,1))) par(mfrow = c(1,length(profiles)))

      if(missing(col)) col <- do.call("rbind", lapply(1:10, function(s) rainbow(length(profiles), s = s / 10)))
      if(is.null(col)) col <- do.call("rbind", lapply(1:10, function(s) rainbow(length(profiles), s = s / 10)))
      
      for(mm in 1:length(profiles)) {
        plot(x = NA, y = NA, ylim = ylim, xlim = range(position), ylab = "", xlab = "position", ...)

        if(!missing(titles))
          title(main = titles[mm], ...)
        
        for(ii in 1:9)
          rect(position - diff(position[1:2]) / 2, profiles[[mm]][,2 + ii], position + diff(position[1:2]) / 2, profiles[[mm]][,22 - ii], col = col[ii,mm], border = col[ii,mm])
        lines(x = position, y = profiles[[mm]][,1], ylim = c(0,1), col = col[10,mm])        
        if(surrounding > 0) {
          seglen <- length(position) / 3
          abline(v = c(weighted.mean(position[0:1 + seglen], w = c(0.1, 0.9)),
                   weighted.mean(position[seglen * 2 + 1:2], w = c(0.9, 0.1))),
                   col = "dark grey", lty = 3, lwd = 2)
        }
      }
    }
  }

.subProfile <- function(subcoord, position, lenMod, mD, samples, cuts, surrounding = 0)
  {
    modcod <- subcoord
    if(surrounding > 0) {
      start(modcod) <- pmax(start(modcod) - surrounding, 1)
      end(modcod) <- end(modcod) + surrounding
    }
            
    methOvers <- findOverlaps(modcod, mD@alignments, select = "all")
    overMD <- sort(unique(unlist(methOvers@subjectHits)))
    redOvers <- split(match(methOvers@subjectHits, overMD), factor(methOvers@queryHits, levels = 1:length(modcod)))

    uus <- unique(unlist(samples))
    redsamp <- lapply(samples, match, uus)
    redMD <- mD[overMD,uus]
    
    if(class(mD) == "alignmentMeth") {      
      return(.averageMethylationRegions(redMD = redMD, subcoord = subcoord, redOvers = redOvers, position = position, samples = redsamp, lenMod = lenMod, cuts = cuts, surrounding = surrounding))
    } else {            
      chrwidths <- sapply(seqlevels(mD@alignments), function(chr) max(c(0, end(modcod[seqnames(modcod) == chr]))))
      
      if(class(mD) == "alignmentData") {
        RPKM <- t(t(getCounts(modcod, redMD, cl = NULL) / width(modcod)) / redMD@libsizes) * 1e9
        sampRPKM <- sapply(redsamp, function(samp) {
          rowMeans(RPKM[,samp,drop = FALSE])
        })
        minRPKM <- max(quantile(sampRPKM[sampRPKM != 0], 0.1), 10)
        
                                        #      RPKM <- t(t(getCounts(subcoord, redMD, cl = NULL) / width(subcoord)) / mD@libsizes) * 1e9
                                        #      end(leftMod[strand(leftMod) != "-"]) <- start(subcoord[strand(leftMod) != "-"]) - 1
                                        #      start(leftMod[strand(leftMod) == "-"]) <- end(subcoord[strand(leftMod) == "-"]) + 1
                                        #      leftRPKM <- t(t(getCounts(leftMod, redMD, cl = NULL) / width(leftMod)) / mD@libsizes) * 1e9
                                        #      start(rightMod[strand(rightMod) != "-"]) <- end(subcoord[strand(rightMod) != "-"]) + 1
                                        #      end(rightMod[strand(rightMod) == "-"]) <- start(subcoord[strand(rightMod) == "-"]) - 1
                                        #      rightRPKM <- t(t(getCounts(rightMod, redMD, cl = NULL) / width(rightMod)) / mD@libsizes) * 1e9
        
        
        codbase <- lapply(seqlevels(redMD@alignments), function(chr) {
          message(".", appendLF = FALSE)
          bases <- as.integer(ranges(modcod[seqnames(modcod) == chr]))
          id <- rep(which(seqnames(modcod) == chr), width(modcod[seqnames(modcod) == chr]))
          adjPos <- bases - start(modcod)[as.integer(id)] + 1                
          centre <- Rle(bases >= start(subcoord)[id] & bases <= end(subcoord)[id])
          leftBound <- start(subcoord)[id] - start(modcod)[id]
          
          adjPos[which(centre)] <- (adjPos[which(centre)] - leftBound[which(centre)]) / (width(subcoord)[id])[which(centre)] * lenMod        
          minus <- Rle(as.character(strand(subcoord))[id] == "-")
          adjPos[which(centre & minus)] <- lenMod - adjPos[which(centre & minus)] + 1
          
          if(surrounding > 0) {
            left <- Rle(bases < start(subcoord)[id])
            right <-Rle(bases > end(subcoord)[id])
            adjPos[which(left)] <- adjPos[which(left)] + surrounding - leftBound[which(left)]
            adjPos[which(right)] <- adjPos[which(right)] - (width(subcoord)[id])[which(right)] - leftBound[which(right)]
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
          
          sampRPKM <- rowMeans(RPKM[,samp,drop = FALSE])
          sampRPKM[sampRPKM < minRPKM] <- minRPKM
          
          weights <- rowSums(t(t(redMD@data[overCod,samp] / redMD@alignments$multireads[overCod]) / redMD@libsizes[samp]) * 1e6) / length(samp)
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
          
                                        #        coverBase <- coverBase[which((coverBase$id %in% which(sampRPKM > 0) & coverBase$centre) |
                                        #                               (coverBase$id %in% which(sampLRPKM > 0) & coverBase$left) |
                                        #                               (coverBase$id %in% which(sampRRPKM > 0) & coverBase$right)),]
          
                                        #        if(any(coverBase$left)) leftSummary <- sapply(split(coverBase$coverage[coverBase$left], factor(cut(coverBase$adjPos[coverBase$left], breaks = cuts, labels = FALSE, include.lowest = TRUE), levels = 1:cuts)), mean, na.rm = TRUE) else leftSummary <- rep(0, cuts)
                                        #        if(any(coverBase$right)) rightSummary <- sapply(split(coverBase$coverage[coverBase$right], factor(cut(coverBase$adjPos[coverBase$right], breaks = cuts, labels = FALSE, include.lowest = TRUE), levels = 1:cuts)), mean, na.rm = TRUE) else rightSummary <- rep(0, cuts)
          
          summariseCoverage <- function(split, splitWeight, codRPKM) {
            sapply(1:length(split), function(ii) {            
              weighted.mean(coverage[split[[ii]]] / codRPKM[as.integer(mergeBase$id[split[[ii]]])], w = splitWeight[[ii]], na.rm = TRUE) * mean(codRPKM, trim = 0.1)
            })
          }
          
          if(surrounding > 0) {
            sumCov <- c(summariseCoverage(leftSplit, leftWeight, sampRPKM),
                        summariseCoverage(centreSplit, centreWeight, sampRPKM),
                        summariseCoverage(rightSplit, rightWeight, sampRPKM))
          } else sumCov <- summariseCoverage(centreSplit, centreWeight, sampRPKM)
          sumCov
        })
      }
            
      message(".", appendLF = FALSE)
      return(repDat)
    }
  }

averageProfiles <- function(mD, samples, coordinates, cuts, maxcuts = 200, bw = 5000, surrounding = 0, add = FALSE, col, ylim, meanOnly = TRUE, plot = TRUE, ...)
  {
    message("Calculating...", appendLF = FALSE)
    
    if(missing(cuts))
      cuts <- max(ceiling(median(width(coordinates)) / bw * length(coordinates)), 5)
#    cuts <- ceiling(cuts)
    cuts <- min(cuts, maxcuts)
    if(missing(samples))
      samples <- mD@replicates
    if(is.factor(samples)) {
      sampNames <- levels(samples)
      samples <- lapply(levels(samples), function(rep) which(samples== rep))
      names(samples) <- sampNames
    }
    if(!is.list(samples)) samples <- list(samples)
    if(missing(col))
      col <- NULL

    if(length(unique(width(coordinates))) == 1) lenMod <- unique(width(coordinates)) else lenMod <- c(1000, surrounding)[as.integer(surrounding != 0) + 1]
    
    if(surrounding > 0) {
      position <- c(0.5 + 0:(cuts - 1) / cuts * surrounding,
                    surrounding + (0.5 + 0:(cuts - 1) / cuts * lenMod),
                    surrounding + lenMod + 0.5 + 0:(cuts - 1) / cuts * surrounding)
    } else position <- (0.5 + 0:(cuts - 1)) / cuts * lenMod      

    if(ceiling(length(coordinates) / 5000) > 1) {
      splitcod <- split(1:length(coordinates), cut(1:length(coordinates), breaks = ceiling(length(coordinates) / 5000), labels = FALSE))
    } else splitcod <- list(1:length(coordinates))
    
    splitProf <- lapply(splitcod, function(x) .subProfile(subcoord = coordinates[x], position = position, lenMod = lenMod, mD = mD, samples = samples, cuts = cuts, surrounding = surrounding))

    profiles <- lapply(1:length(samples), function(ii) Reduce("+", lapply(splitProf, function(x) x[[ii]])) / length(splitProf))
    names(profiles) <- names(samples)
    message(".done!")    
    if(missing(ylim) || is.null(ylim)) ylim <- c(0, max(sapply(profiles, max, na.rm = TRUE), na.rm = TRUE) * 1.1)
    
    if(plot) plotAverageProfile(position = position, profiles = profiles, surrounding = surrounding, ylim = ylim,add = add, meanOnly = meanOnly, ...)
      
    invisible(list(position = position, profiles = profiles))
  }

plotMethDistribution <- function(meth, samples, bw = 1e-3, subtract, chrs, centromeres, add = FALSE, col, legend = TRUE, ...)
  {
    if(missing(samples)) samples <- meth@replicates
    if(missing(chrs)) chrs <- NULL
    if(missing(col)) col <- rainbow(length(samples))
    
    if(is.factor(samples)) {
      namSamp <- levels(samples)
      samples <- lapply(levels(samples), function(rep) which(samples== rep))
      names(samples) <- namSamp
    }
    if(!is.list(samples)) samples <- list(samples) else samples <- samples

    if(missing(subtract)) subtract <- NULL
    if(!is.list(subtract)) subtract <- lapply(1:length(samples), function(rep) subtract)
    
    if(missing(centromeres)) centromeres <- NULL
    
    methDist <- lapply(1:length(samples), function(samps)
                       .plotSampleMeth(meth[,samples[[samps]]], bw = bw, subtract = subtract[[samps]], col = col[samps], centromeres = centromeres, add = c(add, TRUE)[as.integer(samps > 1) + 1], chrs = chrs, ...))

    if(length(samples) > 1 & legend) legend(x = "topright", legend = names(samples), lwd = 2, lty = 1, col = col)
    invisible(methDist)
  }
    
.plotSampleMeth <- function(meth, bw = 1e-2, subtract, centromeres, chrs, add = FALSE, col = col, ...) #adjustByCs = FALSE, ...) {
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
      plot(x = position, y = methylation, type = "l", axes = FALSE, xlim = c(0, max(position) * 1.1), ylim = ylim, col = col, ylab = "", ...)
      axis(2, at = pretty(0:1, n = 5))
      if(length(chrlens) > 1)
        segments(x0 = cumsum(chrlens)[-length(chrlens)], y0 = 0, y1 = 1, col = "red", lty = 2, lwd = 3)    
      if(!missing(centromeres) && !is.null(centromeres)) segments(x0 = c(centromeres[,1] + c(0, cumsum(chrlens)[-length(chrlens)]), centromeres[,2] + c(0, cumsum(chrlens)[-length(chrlens)])), y0 = 0, y1 = 1, lty = 2, lwd = 2, col = "blue")
      text(names(chrlens), srt = 10, adj = 1, y = ylim[1] + 0.05, x = cumsum(c(0, chrlens[-length(chrlens)])) + chrlens / 2, cex = 1.5)
    } else lines(x = position, y = methylation, col = col, ...)
    invisible(data.frame(position = breaks, methylation = methdiv))
  }
