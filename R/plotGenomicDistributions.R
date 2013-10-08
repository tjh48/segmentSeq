averageMethylationRegions <- function(mD, samples, coordinates, cuts, bw = 0.01, surrounding = 0, add = FALSE, col, ...)
  {
    if(missing(cuts)) cuts <- max(length(coordinates) * bw, 20)
    cuts <- ceiling(cuts)
    if(missing(samples)) samples <- mD@replicates
    if(is.factor(samples)) samples <- sapply(levels(samples), function(rep) which(samples== rep))
    if(!is.list(samples)) samples <- list(samples)
    if(missing(col)) col <- rainbow(length(samples))
    modcod <- coordinates
    start(modcod) <- pmax(start(modcod) - surrounding, 0)
    end(modcod) <- end(modcod) + surrounding
    aligns <- mD@alignments
    values(aligns) <- NULL
    methOvers <- findOverlaps(modcod, aligns, select = "all")
    overMD <- sort(unique(unlist(methOvers@subjectHits)))
    redOvers <- split(match(methOvers@subjectHits, overMD), as.factor(methOvers@queryHits))
    redMD <- mD[overMD,]

    uus <- unique(unlist(samples))
    redsamp <- lapply(samples, match, uus)    
    modC <- redMD@Cs[unlist(redOvers),uus,drop = FALSE] / as.integer(redMD@alignments$multireads[unlist(redOvers)])
    modT <- redMD@Ts[unlist(redOvers),uus,drop = FALSE] / as.integer(redMD@alignments$multireads[unlist(redOvers)])
    
    lenOvers <- sapply(redOvers, length)

    redcoord <- coordinates[unique(methOvers@queryHits),]    
    startData <- cbind(start = unlist(redOvers), codstart = rep(start(redcoord), lenOvers), codwidth = rep(width(redcoord), lenOvers))

    startData[,1] <- start(redMD@alignments)[startData[,1]] - startData[,2]
    minus <- (rep(as.character(strand(redcoord)), lenOvers) == "-")
    if(surrounding > 0) {
      rights <- (startData[,1] >= startData[,3])
      insides <- (startData[,1] >= 0 & startData[,1] < startData[,3])
      lefts <- (startData[,1] < 0)    
      startData[insides,1] <- startData[insides,1] / startData[insides,3] * surrounding
      startData[rights,1] <- startData[rights,1] - startData[rights,3] + surrounding

      startData[minus & insides,1] <- surrounding - startData[minus & insides,1]
      startData[minus & lefts] <- startData[minus & lefts] * -1 + surrounding
      startData[minus & rights] <- (startData[minus & rights] - surrounding) * - 1
    } else {
      startData[,1] <- startData[,1] / startData[,3]
      startData[minus] <- 1 - startData[minus]
      if(length(unique(startData[,3])) == 1) lenMod <- unique(startData[,3]) else lenMod <- 1
      startData[,1] <- startData[,1] * lenMod
    }

    if(surrounding > 0) {
      splits <- cut(startData[,1], breaks = (-cuts):(2 * cuts) / cuts * surrounding, include.lowest = TRUE, labels = FALSE)
      position <- (0.5 + (-cuts):(2 * cuts - 1)) /cuts * surrounding
      weightings <- rep(1, nrow(startData))
      weightings[startData[,1] >= 0 & startData[,1] < surrounding] <- 1 / startData[startData[,1] >= 0 & startData[,1] < surrounding,3]
    } else {      
      splits <- cut(startData[,1], breaks = 0:cuts / cuts * lenMod, include.lowest = TRUE, labels = FALSE)
      position <- (0.5 + 1:(cuts - 1)) / cuts
      weightings <- 1 / startData[,3]
    }

    methProps <- lapply(redsamp, function(samp) {
      Cs <- modC[, samp, drop = FALSE]
      Ts <- modT[, samp, drop = FALSE]
      sampProps <- rowMeans(Cs / (Cs + Ts), na.rm = TRUE)
      sampWeights <- weightings
      sampWeights[is.na(sampProps)] <- 0
      z <- sapply(
             split(1:length(sampProps), factor(splits, levels = 1:c(3 * cuts, cuts)[as.numeric(surrounding == 0) + 1])),
             function(zzz) weighted.mean(sampProps[zzz], w = sampWeights[zzz]))
    })

    lapply(1:length(methProps), function(ii) {
      if(ii == 1 & add == FALSE) {
        methylation <- methProps[[ii]]
        plot(x = position, y = methylation, type= "l", col = col[ii], ...)
      } else lines(x = position, y = methProps[[ii]], col = col[ii], ...)
    })
    if(!add) {
      if(surrounding > 0) abline(v = c(0, surrounding), col = "orange", lty = 3)
      legend(x = "topright", lty = 1, col = col, legend = names(samples))
    }

    invisible(list(position, methProps))
  }
    
plotMethDistribution <- function(meth, samples, bw = 1e-3, subtract, centromeres, add = FALSE, ...) #adjustByCs = FALSE, ...) {
  {
    chrlens <- seqlengths(meth@alignments)
    if(any(is.na(chrlens)))
      chrlens[is.na(chrlens)] <- sapply(names(chrlens)[is.na(chrlens)], function(nana) max(end(meth@alignments[seqnames(meth@alignments) == nana])))
    
    CHGpos <- (end(meth@alignments) + start(meth@alignments)) / 2
    CHGpos <- CHGpos + c(0, cumsum(chrlens))[match(as.character(seqnames(meth@alignments)), seqlevels(meth@alignments))]
    adjCs <- rowSums(meth@Cs[,samples,drop = FALSE] / as.integer(meth@alignments$multireads))
    adjTs <- rowSums(meth@Ts[,samples,drop = FALSE] / as.integer(meth@alignments$multireads))
    
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
    
    
    ints <- findInterval(breaks, c(CHGpos))
    sumC <- cumsum(adjCs)[ints]; sumT <- cumsum(adjTs)[ints]
    countC <- sumC - c(0, sumC[-length(sumC)]); countT <- sumT - c(0, sumT[-length(sumT)])
    methdiv <- (countC / (countC + countT))
    
                                        #  if(adjustByCs) {
                                        #    siteNums <- findInterval(breaks, unique(CHGpos))
                                        #    Cadj <- siteNums - c(0, siteNums[-length(breaks)])
                                        #    methdiv <- methdiv * Cadj / max(Cadj)
                                        #  }
    
    if(!missing(subtract)) {
      methdiv <- methdiv - subtract
      ylim <- c(-1.1, 1)
    } else ylim = c(-0.1, 1)
    
    methylation <- methdiv[!is.na(methdiv)]
    position <- breaks[!is.na(methdiv)] #- cuts / 2
    if(!add) {
      plot(x = position, y = methylation, type = "l", axes = FALSE, ylim = ylim, ...)
      axis(2)
      abline(v = cumsum(chrlens)[-length(chrlens)], col = "red", lty = 2, lwd = 3)    
      if(!missing(centromeres)) abline(v = c(centromeres[,1] + c(0, cumsum(chrlens)[-length(chrlens)]), centromeres[,2] + c(0, cumsum(chrlens)[-length(chrlens)])), lty = 2, lwd = 2, col = "blue")
      text(seqlevels(meth@alignments), srt = 30, adj = 1, y = ylim[1] + 0.05, x = cumsum(c(0, chrlens[-length(chrlens)])) + chrlens / 2)
    } else lines(x = position, y = methylation, ...)
    invisible(data.frame(position = breaks, methylation = methdiv))
  }   
