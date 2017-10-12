% modification on git from copied files
.zeroInMeth <- function(aD, smallSegs)
    {

        
    zeroCs <- !is.na(findOverlaps(aD@alignments, smallSegs, select = "arbitrary")) #.getOverlaps(aD@alignments, smallSegs, whichOverlaps = FALSE)
#    zeroCs <- rowSums(sapply(1:ncol(aD), function(ii) as.integer(aD@Cs[,ii]))) == 0
    
    chrBreaks <- which(seqnames(aD@alignments)[-nrow(aD)] != seqnames(aD@alignments)[-1])
    chrBreaks <- cbind(c(1, chrBreaks + 1), c(chrBreaks, nrow(aD)))
    
    whichZero <- which(zeroCs)

    if(length(whichZero) > 0) {
      zeroBlocks <- do.call("rbind", lapply(1:nrow(chrBreaks), function(ii) {
        chrZeros <- whichZero[whichZero >= chrBreaks[ii,1] & whichZero <= chrBreaks[ii,2]]      
        adjZeros <- cbind(chrZeros[c(1, which(diff(chrZeros) > 1) + 1)], chrZeros[c(which(diff(chrZeros) > 1), length(chrZeros))])
      }))
      
      zeroCoords <- GRanges(seqnames(aD@alignments)[zeroBlocks[,1]], IRanges(start = start(aD@alignments)[zeroBlocks[,1]], end = end(aD@alignments)[zeroBlocks[,2]]))
    } else {
      zeroCoords <- GRanges()
      seqinfo(zeroCoords) <- seqinfo(aD@alignments)
    }
      
    zeroCoords
  }      


.constructMethNulls <- function(emptyNulls, sDP, locDef, minlen)
    {
        
        overLoc <- sDP@coordinates[!is.na(findOverlaps(sDP@coordinates, locDef, type = "within", select = "arbitrary")),]
                                        #.getOverlaps(sDP@coordinates, locDef, overlapType = "within", whichOverlaps = FALSE, ignoreStrand = FALSE, cl = NULL),]
                     
    
        leftRight <- do.call("rbind", lapply(seqlevels(overLoc), function(chr) {
            leftids <- Rle(findInterval(start(overLoc[seqnames(overLoc) == chr,]), end(emptyNulls[seqnames(emptyNulls) == chr,])))
            leftids[leftids <= 0] <- NA      
            left <- Rle(start(overLoc[which(seqnames(overLoc) == chr),]) - rep(start(emptyNulls[seqnames(emptyNulls) == chr])[runValue(leftids)], runLength(leftids)))
            
            rightids <- Rle(findInterval(end(overLoc[which(seqnames(overLoc) == chr)]), end(emptyNulls[seqnames(emptyNulls) == chr,])) + 1)
            rightids[rightids > sum(seqnames(emptyNulls) == chr)] <- NA
            right = Rle(rep(end(emptyNulls[seqnames(emptyNulls) == chr])[runValue(rightids)], runLength(rightids)) - end(overLoc[which(seqnames(overLoc) == chr),]))
            
            DataFrame(left, right)
        }))
        
        leftGood <- (!is.na(leftRight[,'left']))
        rightGood <- (!is.na(leftRight[,'right']))
        
        nullCoords = GRanges(
            seqnames = c(seqnames(emptyNulls), seqnames(overLoc)[c(which(leftGood & rightGood), which(leftGood), which(rightGood))]),
            IRanges(
                start = as.integer(c(Rle(start(emptyNulls)),
                    (start(overLoc) - leftRight[,"left"])[which(leftGood & rightGood)],
                    (start(overLoc) - leftRight[,"left"])[which(leftGood)],
                    start(overLoc)[which(rightGood)])),
                end = as.integer(c(Rle(end(emptyNulls)),
                    (end(overLoc) + leftRight[,"right"])[which(leftGood & rightGood)],
                    (end(overLoc))[which(leftGood)],
                    (end(overLoc) + leftRight[,'right'])[which(rightGood)]))
            )
        )
        rm(leftRight)
        if(!missing(minlen)) nullCoords <- nullCoords[width(nullCoords) >= minlen,]
        if(length(nullCoords) > 0) {
            nullCoords <- unique(nullCoords)
            nullCoords <- nullCoords[order(as.integer(seqnames(nullCoords)), start(nullCoords), end(nullCoords))]
            nullCoords <- nullCoords[which(.getOverlaps(coordinates = nullCoords, segments = locDef, overlapType = "within", whichOverlaps = FALSE, cl = NULL))]
            # for some reason this is incredibly slow using findOverlaps
            #nullCoords <- nullCoords[which(!is.na(findOverlaps(nullCoords, locDef, type = "within", select = "first")))]
        }
        potnullD <- new("lociData",
                        replicates = sDP@replicates,
                        coordinates = nullCoords, locLikelihoods = matrix(NA, nrow = 0, ncol = nlevels(sDP@replicates)))
        potnullD
    }

.constructNulls <- function(emptyNulls, sDWithin, locDef, forPriors = FALSE, samplesize, aD = aD, cl = NULL)
  {
    # find the gap to the left and right of each element of sD within a locus (sDWithin)
    leftRight <- matrix(NA, ncol = 2, nrow = nrow(sDWithin))
    colnames(leftRight) <- c("left", "right")

    for(ss in unique(strand(sDWithin@coordinates)))      
        leftRight[which(strand(sDWithin@coordinates) == ss),] <- do.call("rbind",
                                  lapply(unique(seqnames(sDWithin@coordinates)), function(chr) {
                                      sDsel <- which(seqnames(sDWithin@coordinates) == chr & strand(sDWithin@coordinates) == ss)
                                      empsel <- which(seqnames(emptyNulls) == chr & strand(emptyNulls) == ss)
                                      left <- start(sDWithin@coordinates[sDsel,]) -
                                          start(emptyNulls[empsel])[match(start(sDWithin@coordinates[sDsel,]), (end(emptyNulls[empsel,]) + 1))]
                                      right <- end(emptyNulls[empsel])[match(end(sDWithin@coordinates[sDsel,]), (start(emptyNulls[empsel,]) - 1))] -
                                          end(sDWithin@coordinates[sDsel,])
                                      cbind(left, right)
                                  }))

    # select from the empty regions those within or adjacent (unless looking to construct priors) to a potential locus
    
    empties <- which(!is.na(findOverlaps(emptyNulls, sDWithin@coordinates, type = "within", select = "first")))

    if(is.na(forPriors) || !forPriors) empties <- union(empties,
                                                         unlist(lapply(levels(seqnames(sDWithin@coordinates)), function(chr) {
                                                             leftEmpties <- match(start(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (end(emptyNulls[seqnames(emptyNulls) == chr,]) + 1))
                                                             rightEmpties <- match(end(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (start(emptyNulls[seqnames(emptyNulls) == chr,]) - 1))
                                                             which(seqnames(emptyNulls) == chr)[union(leftEmpties, rightEmpties)]
                                                         }))
                                                         )    
    empties <- empties[!is.na(empties)]

    leftGood <- (!is.na(leftRight[,'left']))
    rightGood <- (!is.na(leftRight[,'right']))

    # combine regions with gaps to the left, gaps to the right, gaps on both side, plus empty regions within loci
    
    nullCoords = GRanges(
      seqnames = c(
        seqnames(emptyNulls[empties,]),
        seqnames(sDWithin@coordinates)[c(which(leftGood & rightGood), which(leftGood), which(rightGood))]),
      IRanges(
              start = c(start(emptyNulls[empties,]),
                (start(sDWithin@coordinates) - leftRight[,"left"])[which(leftGood & rightGood)],
                (start(sDWithin@coordinates) - leftRight[,"left"])[which(leftGood)],
                (start(sDWithin@coordinates))[which(rightGood)]),
              end = c(end(emptyNulls[empties,]),
                (end(sDWithin@coordinates) + leftRight[,"right"])[which(leftGood & rightGood)],
                (end(sDWithin@coordinates))[which(leftGood)],
                (end(sDWithin@coordinates) + leftRight[,'right'])[which(rightGood)])),
      sDID = c(rep(-1, length(empties)), which(leftGood & rightGood), which(leftGood), which(rightGood)),
      strand = c(
        strand(emptyNulls[empties,]),
        strand(sDWithin@coordinates)[c(which(leftGood & rightGood), which(leftGood), which(rightGood))])
      )
      

                                        # which constructed coordinates lie within a locus?
    nullCoords <- nullCoords[which(.getOverlaps(coordinates = nullCoords, segments = locDef, overlapType = "within", whichOverlaps = FALSE, cl = NULL))]
    #nullCoords <- nullCoords[!is.na(findOverlaps(nullCoords, locDef, type = "within", select = "first")),]


    if(is.na(forPriors)) {
        if(length(nullCoords) == 0) return(nullCoords) else return(nullCoords[sample(1:length(nullCoords), 1)])        
    }
    if(forPriors) {
        nullCoords <- nullCoords[.filterSegments(nullCoords, runif(length(nullCoords)), maxReport = samplesize),]
    }
                                        # empty regions carry no data    
    
    if(class(aD) == "alignmentData") {
      nullData <- matrix(0L, nrow = length(nullCoords), ncol = ncol(sDWithin))
      colnames(nullData) <- colnames(sDWithin@data)
      if(nrow(sDWithin@data) > 0) nullData[nullCoords$sDID > 0,] <- sDWithin@data[nullCoords$sDID[nullCoords$sDID > 0],]
      values(nullCoords) <- NULL
      potnullD <- new("lociData", data = nullData,
                      libsizes = libsizes(sDWithin),
                      replicates = sDWithin@replicates,
                      coordinates = nullCoords)
    } else if(class(aD) == "alignmentMeth") {
      values(nullCoords) <- NULL
      potnullD <- new("lociData",
                      data = list(matrix(NA, ncol = ncol(sDWithin), nrow = length(nullCoords)), 
                          matrix(NA, ncol = ncol(sDWithin), nrow = length(nullCoords))),        
                      replicates = sDWithin@replicates,
                      coordinates = nullCoords)      
      nullCounts <- getCounts(potnullD@coordinates, aD, cl = cl)
      potnullD@Cs <- nullCounts$Cs
      potnullD@Ts <- nullCounts$Ts
      potnullD <- potnullD[rowSums(potnullD@Cs) + rowSums(potnullD@Ts) > 0,]
    }
      
    potnullD
  }


.constructNullPriors <- function(emptyNulls, sDWithin, locDef, samplesize)
  {
    
    leftRight <- do.call("rbind", lapply(seqlevels(sDWithin@coordinates), function(chr) {
      left <- start(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]) -
        start(emptyNulls[seqnames(emptyNulls) == chr])[match(start(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (end(emptyNulls[seqnames(emptyNulls) == chr,]) + 1))]
      right <- end(emptyNulls[seqnames(emptyNulls) == chr])[match(end(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),]), (start(emptyNulls[seqnames(emptyNulls) == chr,]) - 1))] - end(sDWithin@coordinates[which(seqnames(sDWithin@coordinates) == chr),])
      cbind(left, right)
    }))

    empover <- !is.na(findOverlaps(emptyNulls, sDWithin@coordinates, type = "within", select = "first"))
    if(any(empover)) emptyNulls <- emptyNulls[which(empover),] else emptyNulls <- emptyNulls[FALSE,]
    
    leftGood <- !is.na(leftRight[,'left'])
    rightGood <- !is.na(leftRight[,'right'])    
    
    nullCoords <- GRanges(
      seqnames = c(
        seqnames(emptyNulls),
        seqnames(sDWithin@coordinates)[leftGood & rightGood],
        seqnames(sDWithin@coordinates)[leftGood],
        seqnames(sDWithin@coordinates)[rightGood]),
                    IRanges(
                      start = c(
                        start(emptyNulls),
                        (start(sDWithin@coordinates) - leftRight[,"left"])[leftGood & rightGood],
                        (start(sDWithin@coordinates) - leftRight[,"left"])[leftGood],
                        (start(sDWithin@coordinates))[rightGood]),              
                      end = c(end(emptyNulls),
                        (end(sDWithin@coordinates) + leftRight[,"right"])[leftGood & rightGood],
                        (end(sDWithin@coordinates))[leftGood],
                        (end(sDWithin@coordinates) + leftRight[,'right'])[rightGood])
                      ))
    
    overLoci <- which(!is.na(findOverlaps(coordinates = nullCoords, segments = locDef, type = "within", select = "first")))
    filNulls <- sort(overLoci[.filterSegments(nullCoords[overLoci], runif(length(overLoci)), maxReport = samplesize)])
    splitNulls <- c(length(emptyNulls), sum(leftGood & rightGood), sum(leftGood), sum(rightGood))
    emptyData <- do.call("cbind", lapply(1:ncol(sDWithin), function(x) (rep(0L, sum(filNulls <= splitNulls[1])))))    

    lrChoose <- which(leftGood & rightGood)[filNulls[filNulls > cumsum(splitNulls)[1] & filNulls <= cumsum(splitNulls)[2]] - cumsum(splitNulls)[1]]
    leftChoose <- which(leftGood)[filNulls[filNulls > cumsum(splitNulls)[2] & filNulls <= cumsum(splitNulls)[3]] - cumsum(splitNulls)[2]]
    rightChoose <- which(rightGood)[filNulls[filNulls > cumsum(splitNulls)[3] & filNulls <= cumsum(splitNulls)[4]] - cumsum(splitNulls)[3]]
    
    if(class(sDWithin) == "segData") {
      colnames(emptyData) <- colnames(sDWithin@data)
      potnullD <- new("segData", data = rbind(emptyData, sDWithin@data[c(lrChoose, leftChoose, rightChoose),]),
                      libsizes = libsizes(sDWithin),
                      replicates = sDWithin@replicates,
                      coordinates = nullCoords[filNulls])
    } else if(class(sDWithin) == "segMeth") {
      colnames(emptyData) <- colnames(sDWithin@Cs)
      potnullD <- new("segMeth",
                      Cs = rbind(emptyData, sDWithin@Cs[c(lrChoose, leftChoose, rightChoose),]),
                      Ts = rbind(emptyData, sDWithin@Ts[c(lrChoose, leftChoose, rightChoose),]),
                      replicates = sDWithin@replicates,
                      coordinates = nullCoords[filNulls])
    }
    
    potnullD <- potnullD[order(as.integer(seqnames(potnullD@coordinates)), start(potnullD@coordinates), end(potnullD@coordinates)),]
    
    potnullD
  }
