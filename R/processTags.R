.tagDuplicated <- function(x)
    if(!"tag" %in% colnames(values(x))) return(duplicated(x)) else duplicatedIntegerQuads(seqnames(x), start(x), strand(x), match(x$tag, unique(x$tag)))

readMeths <- function(files, dir = ".", libnames, replicates, nonconversion, chrs)#, splitStrands = TRUE)#, duplicateFiles)
  {
    if(missing(nonconversion))
      nonconversion = rep(0, length(files))
    if(any(nonconversion < 0) | any(nonconversion > 1)) stop("Non-conversion rates must be between zero and one.")
    
    replicates <- as.factor(replicates)
    if(dir != "") files <- paste(dir, files, sep = "/")

    if(missing(chrs)) chrs <- NULL
    
    message("Reading files...", appendLF = FALSE)
    methReads <- lapply(files, function(file, chrs) {
        message(".", appendLF = FALSE)
        mreads <- read.delim(file, as.is = TRUE, header = FALSE)
        mreads <- mreads[mreads[,2] >= 0,]
        if(!is.null(chrs)) mreads <- mreads[mreads[,1] %in% chrs,]
        mreads <- mreads[mreads[,3] %in% c("+", "-", "*"),]
      mr <- GRanges(seqnames = mreads[,1], IRanges(mreads[,2], width = 1), strand = mreads[,3], #multireads = as.integer(mreads[,6]),
                    Cs = as.integer(mreads[,4]), Ts = as.integer(mreads[,5]))
      if(ncol(mreads) == 6) mr$multireads <- Rle(mreads[,6]) else mr$multireads <- Rle(1)
      mr
    }, chrs = chrs)
    message("done!", appendLF = TRUE)

    chrs <- sort(unique(do.call("c", lapply(methReads, seqlevels))))
    methReads <- lapply(methReads, function(gr) {
      chrnames <- factor(as.character(seqnames(gr)), levels = chrs)
      seqlevels(gr) <- chrs
      seqnames(gr) <- chrnames
      gr
    })
    
    chrTag <- methReads[[1]]

    if(length(methReads) > 1) {
      message("Finding unique cytosines...", appendLF = FALSE)      
      for(ii in 2:length(methReads)) {
        message(".", appendLF = FALSE)
        chrTag <- append(chrTag, methReads[[ii]])
        chrTag <- chrTag[order(as.integer(seqnames(chrTag)), as.integer(start(chrTag)), as.integer(strand(chrTag)), as.integer(values(chrTag)$multireads)),]
        chrTag <- chrTag[!duplicated(chrTag) | c(TRUE, as.integer(values(chrTag)$multireads[-length(chrTag)]) != as.integer(values(chrTag)$multireads[-1])),]
      }    
      message("done!", appendLF = TRUE)
    }

    values(chrTag)$Cs <- values(chrTag)$Ts <- 0      
    message("Processing samples...", appendLF = FALSE)      
    mReads <- lapply(methReads, function(mreads) {
      message(".", appendLF = FALSE)
      colnames(chrTag) <- colnames(mreads)
      mreads <- append(mreads, chrTag)
      mreads <- mreads[order(as.integer(seqnames(mreads)), as.integer(start(mreads)), as.integer(strand(mreads)), as.integer(values(mreads)$multireads)),]
      mreads <- mreads[!duplicated(mreads) | c(TRUE, as.integer(values(mreads)$multireads[-length(mreads)]) != as.integer(values(mreads)$multireads[-1])),]
      list(Cs = values(mreads)$Cs, Ts = values(mreads)$Ts)
    })
      
    message("done!", appendLF = TRUE)
    
    Cs <- do.call("cbind", lapply(mReads, function(x) x$Cs))
    Ts <- do.call("cbind", lapply(mReads, function(x) x$Ts))
    
    colnames(Cs) <- libnames
    colnames(Ts) <- libnames
    values(chrTag) <- values(chrTag)$multireads
    colnames(values(chrTag)) <- "multireads"
    
    mD = new("alignmentMeth", alignments = chrTag, Cs = Cs, Ts = Ts, replicates = replicates, libnames = libnames, nonconversion = nonconversion)
    
                                        #    if(!missing(duplicateFiles))
                                        #      {
                                        #        message("Reading duplicate files...", appendLF = FALSE)
                                        #        dupReads <- lapply(duplicateFiles, function(file) {
                                        #          message(".", appendLF = FALSE)
                                        #          dupreads <- read.delim(file, as.is = TRUE, header = FALSE)
                                        #          GRanges(seqnames = dupreads[,1], IRanges(dupreads[,2], width = 1), strand = dupreads[,3], multiread = Rle(dupreads[,6]), Cs = Rle(dupreads[,4]), Ts = Rle(#dupreads[,5]), id = Rle(dupreads[,7]))
                                        #        })
                                        #        message("done!", appendLF = TRUE)
    ##        mD@duplication <- dupReads
                                        #     } else mD@duplication <- NULL
    
                                        #            if(splitStrands) mD <- strandSplitter(mD)
    
    mD
  }

.fastUniques <- function(x, na = FALSE){
    if (nrow(x) > 1) {
          if(na) {
                  return(c(TRUE, (rowSums(x[-1L,,drop = FALSE] != x[-nrow(x),,drop = FALSE], na.rm = TRUE) > 0) | rowSums((is.na(x[-1L,,drop= FALSE])== is.na(x[-nrow(x),,drop = FALSE]))) != ncol(x)))
                } else return(c(TRUE, rowSums(x[-1L, , drop = FALSE] == x[-nrow(x),, drop = FALSE]) != ncol(x)))
        } else return(TRUE)
  }


readBAM <-
function(files, dir = ".", replicates, libnames, chrs, chrlens, countID = NULL, minlen = 15, maxlen = 1000, multireads = 1000, polyLength, estimationType = "quantile", discardTags = FALSE, verbose = TRUE, filterReport = NULL)
  {
    if(missing(polyLength)) polyLength <- NULL
    if(!is.null(polyLength)) polyBase <- do.call("rbind", lapply(c("A", "C", "G", "T"), function(polybase) {
      suppressWarnings(apply(matrix(c(".", rep(polybase, polyLength + 1)),
                                    nrow = polyLength + 1, ncol = polyLength + 1), 1, paste, sep = "", collapse = ""))
    }))
    if(length(files) != length(replicates)) stop("The 'replicates' vector needs to be the same length as the 'files' vector")

    if(missing(chrs)) chrs <- NULL else {
        chrs <- as.character(chrs)
        seqinf <- Seqinfo(seqnames = chrs)#, seqlengths = chrlens)
    }
    
    replicates <- as.factor(replicates)

    

    if(!all(levels(replicates) %in% replicates))
      stop("There appear to be additional levels in your (factor) replicates specification which are not present in the replicates vector.")
    
#    if(any(chrlens != as.integer(chrlens)))
#      stop("The 'chrlens' vector must be castable as an integer")
#    chrlens <- as.integer(chrlens)    

    if(missing(libnames))
      libnames <- sub(".*/", "", files)

    if(dir != "") files <- paste(dir, files, sep = "/")

#    if(class(chrs) != "character")
#      stop("'chrs' must be of type 'character'.")

    if(verbose)
      message("Reading files...", appendLF = FALSE)

    sampleNumbers <- 1:length(files)

    Tags <- lapply(sampleNumbers, function(ii) {
      tags <- scanBam(files[ii])[[1]]
      if(!is.null(countID)) {
          counts <- as.integer(scanBam(files[ii], param = ScanBamParam(tag=countID))[[1]][[1]][[1]])
      } else counts <- (rep(1L, length(tags$seq)))
           
      keepReads <- which(!is.na(tags$pos))      
      ir <- IRanges(start = as.integer(tags$pos[keepReads]), width = tags$qwidth[keepReads])
      aln <- GRanges(seqnames = tags$rname[keepReads], ir, strand = tags$strand[keepReads], count = counts[keepReads])
      
      if(!discardTags) {
          revcomp <- sapply(tags$flag, function(flag) intToBits(flag)[5] == 1)
          if(any(revcomp)) tags$seq[revcomp] <- reverseComplement(tags$seq[revcomp])
          aln$tag = (as.character(tags$seq[keepReads]))
          aln <- aln[order(as.factor(seqnames(aln)), as.integer(start(aln)), as.factor(values(aln)$tag)),]
      } else aln <- sort(aln)

      dupTags <- .tagDuplicated(aln)
      if(any(dupTags)) {
          count <- diff(c(which(!dupTags), length(dupTags) + 1))
          if(!is.null(countID)) count <- sapply(split(aln$count, rep(1:length(count), count)), sum)
          aln <- aln[!dupTags,]
          values(aln)$count <- count
      }

      
      filterTags <- rep(NA, 4)     
      filterInfo <- function(seltags, filname) {
          if(!is.null(filterReport) & !discardTags) {
              write(unique(as.character(aln$tag[!seltags])), file = paste(gsub(".*/", "", files[ii]), filterReport, filname, sep = "_"))
              return(paste(length(unique(aln$tag[!seltags])), sum(!seltags), sep = ":"))
          } else return(NA)
      }
      
      if(!is.null(chrs)) {
          chrtags <- as.integer(seqnames(aln)) %in% which(seqlevels(aln) %in% chrs)
          filterTags[1] <- filterInfo(chrtags, "chrs")
          aln <- aln[chrtags,]
      } else filterTags[1] <- NA
      if(length(aln) == 0) warning(paste("There were no tags left in sample", ii, "after selection by chromosome. Are you sure you've got the right chromosome names?"))

      widths <- width(aln)
      goodwidths <- widths >= minlen & widths <= maxlen
      filterTags[2] <- filterInfo(goodwidths, "widths")
      aln <- aln[goodwidths,]      
      
      if(!is.null(multireads) & !discardTags) {
        tabtags <- table(as.character(aln$tag))
        goodmult <- as.character(aln$tag) %in% names(tabtags[tabtags < multireads])
        filterTags[3] <- filterInfo(goodmult, "multireads")
        aln <- aln[goodmult,]
      }

      if(!is.null(polyLength) & !discardTags) {        
        polyBaseRemove <- unique(unlist(lapply(polyBase, grep, as.character(values(aln)$tag))))
        goodpolyBase <- rep(TRUE, length(aln))
        if(length(polyBaseRemove) > 0) goodpolyBase[polyBaseRemove] <- FALSE
        filterTags[4] <- filterInfo(goodpolyBase, "poly")
        aln <- aln[goodpolyBase,]
      }

      if(!is.null(filterReport)) {
        filDat <- matrix(c(file = files[ii], filterTags), nrow = 1)
        colnames(filDat) <- c("file", "chrs", "length", "multireads", "polyLength")
        write.table(filDat, file = paste(filterReport, ".txt", sep = ""), append = (ii > 1), row.names = FALSE, col.names = ii == 1, sep = "\t", quote = FALSE)
      }
      
      message(".", appendLF = FALSE)
      if(!is.null(chrs)) seqinfo(aln, new2old = match(seqlevels(seqinf), seqlevels(aln))) <- seqinf
      aln
  })
    
    message(".done!")
                                        #if(!missing(tempFile)) save(Tags, file = paste("tags_", tempFile, sep = ""))
    if(is.null(chrs)) seqinf <- seqinfo(Tags[[1]])
    
    aD <- .processTags(Tags, verbose = verbose, estimationType = estimationType, seqinf = seqinf, libnames = libnames, replicates = replicates, discardTags = discardTags)
    aD
    
}


readGeneric <-
    function(files, dir = ".", replicates, libnames, chrs, chrlens,
             cols, header = TRUE, minlen = 15, maxlen = 1000, multireads = 1000, polyLength, estimationType = "quantile", discardTags = FALSE, verbose = TRUE, filterReport = NULL, ...)
        {
            if(discardTags) stop("readGeneric does not currently support discarding tag information. Either include this data, or switch to BAM files (recommended).")
            if(missing(polyLength)) polyLength <- NULL
            if(!is.null(polyLength)) polyBase <- do.call("rbind", lapply(c("A", "C", "G", "T"), function(polybase) {
                suppressWarnings(apply(matrix(c(".", rep(polybase, polyLength + 1)),
                                              nrow = polyLength + 1, ncol = polyLength + 1), 1, paste, sep = "", collapse = ""))
            }))
            
            if(missing(chrs)) chrs <- NULL else {
                chrs <- as.character(chrs)
                seqinf <- Seqinfo(seqnames = chrs)#, seqlengths = chrlens)
            }
            
            replicates <- as.factor(replicates)
            
            if(!all(levels(replicates) %in% replicates))
                stop("There appear to be additional levels in your (factor) replicates specification which are not present in the replicates vector.")


            if(!missing(chrlens) & !missing(chrs)) {
                if(any(chrlens != as.integer(chrlens)))
                    stop("The 'chrlens' vector must be castable as an integer")
                chrlens <- as.integer(chrlens)
            }
            
            countPresent <- TRUE
            tagPresent <- TRUE
            strandPresent <- TRUE
            
    if(missing(cols)) cols <- NULL
            
    if(!is.null(cols))
      {
        if(!(all(c("chr", "start", "end") %in% names(cols))))
          stop("'cols' argument must contain named values for 'chr', 'start', 'end' or be NULL")
        if(any(c(!("count" %in% names(cols)), is.na(cols[names(cols) == "count"]))))
          {
            countPresent <- FALSE
            warning("No 'count' column specified in 'cols' argument; I'll assume that the file contains non-redundant reads")
          }
        if(any(c(!("tag" %in% names(cols)), is.na(cols[names(cols) == "tag"]))))
          {
            tagPresent <- FALSE
            warning("No 'tag' column specified in 'cols' argument; the 'alignData' object will omit sequence information.")
          }
        if(any(c(!("strand" %in% names(cols)), is.na(cols[names(cols) == "strand"]))))
          {
            strandPresent <- FALSE
          }        
      } else if(is.null(cols) & header == FALSE) warning("No 'cols' argument supplied and 'header = FALSE'. Using default values for columns")

    if(missing(libnames))
        libnames <- sub(".*/", "", files)

    files <- paste(dir, files, sep = "/")

    if(verbose)
      message("Reading files...", appendLF = FALSE)

    sampleNumbers <- 1:length(files)
    
    Tags <- lapply(sampleNumbers, function(ii, cols, header, ...) {
      filetags <- read.table(files[ii], header = header, as.is = TRUE)
      if(header & is.null(cols))
        {
          if(all(c("chr", "start", "end") %in% names(filetags)))
            {
              chrcol <- which(names(filetags) == "chr")
              startcol <- which(names(filetags) == "start")
              endcol <- which(names(filetags) == "end")
            } else stop(paste("Couldn't find appropriate column names (and columns were not specified) in file:", files[ii]))
          if("tag" %in% names(filetags))
            {
              tagcol <- which(names(filetags) == "tag")
            } else {
              tagPresent <- FALSE
              warning("No 'tag' column found in file; the 'alignData' object will omit sequence information.")
            }
          if("count" %in% names(filetags))
            {
              countcol <- which(names(filetags) == "count")
            } else {
              countPresent <- FALSE
              warning("No 'count' column found in file; I'll assume that the file contains non-redundant reads")
            }
          if("strand" %in% names(filetags)) strandcol <- which(names(filetags) == "strand") else strandPresent <- FALSE
            
        } else if(!header & is.null(cols))
          {
            chrcol <- 1L
            tagcol <- 2L
            countcol <- 3L
            startcol <- 4L
            endcol <- 5L
            strandcol <- 6L            
            tagPresent <- TRUE
            countPresent <- TRUE
            strandPresent <- TRUE
          } else if(!is.null(cols)) {
            chrcol <- cols[names(cols) == "chr"]
            if(tagPresent) tagcol <- cols[names(cols) == "tag"] else tagcol <- NA          
            if(countPresent) countcol <- cols[names(cols) == "count"] else countcol <- NA
            if(strandPresent) strandcol <- cols[names(cols) == "strand"] else strandcol <- NA
            startcol <- cols[names(cols) == "start"]
            endcol <- cols[names(cols) == "end"]
          }

#      if(strandPresent & !is.na(strand)) filetags <- filetags[filetags[,strandcol] %in% strand,]

      filterTags <- rep(NA, 4)     
      filterInfo <- function(seltags, filname) {
        if(!is.null(filterReport)) {
          write(unique(filetags[!seltags, tagcol]), file = paste(gsub(".*/", "", files[ii]), filterReport, filname, sep = "_"))
          if(tagPresent) return(paste(length(unique(filetags[!seltags,tagcol])), sum(!seltags), sep = ":")) else return(sum(!seltags))
        } else return(NA)
      }

      if(!is.null(chrs)) {
          chrtags <- filetags[,chrcol] %in% chrs        
          filterTags[1] <- filterInfo(chrtags, "chrs")
          filetags <- filetags[chrtags,]  
      } else filterTags[1] <- NA
      
      widths <- as.integer(filetags[, endcol]) - as.integer(filetags[,startcol]) + 1      
      goodwidths <- widths >= minlen & widths <= maxlen
      filterTags[2] <- filterInfo(goodwidths, "widths")
      filetags <- filetags[goodwidths,]
      
      if(tagPresent & !is.null(multireads)) {
        tabtags <- table(filetags[,tagcol])
        goodmult <- filetags[,tagcol] %in% names(tabtags[tabtags < multireads])
        filterTags[3] <- filterInfo(goodmult, "multireads")
        filetags <- filetags[goodmult,]
      }
      
      if(tagPresent & !is.null(polyLength)) {        
        polyBaseRemove <- unique(unlist(lapply(polyBase, grep, filetags[,tagcol])))
        goodpolyBase <- rep(TRUE, nrow(filetags))
        if(length(polyBaseRemove) > 0) goodpolyBase[polyBaseRemove] <- FALSE
        filterTags[4] <- filterInfo(goodpolyBase, "poly")
        filetags <- filetags[goodpolyBase,]
      }

      if(!is.null(filterReport)) {
        filDat <- matrix(c(file = files[ii], filterTags), nrow = 1)
        colnames(filDat) <- c("file", "chrs", "length", "multireads", "polyLength")
        write.table(filDat, file = paste(filterReport, ".txt", sep = ""), append = (ii > 1), row.names = FALSE, col.names = ii == 1, sep = "\t", quote = FALSE)
      }
      
      ir <- IRanges(start = as.integer(filetags[,startcol]), end = as.integer(filetags[, endcol]))      
      
      if(tagPresent & countPresent)
        {
          aln <- GRanges(seqnames = filetags[, chrcol], ir, tag = (filetags[, tagcol]), count = (as.integer(filetags[, countcol])))
          if(strandPresent) strand(aln) <- Rle(as.character(filetags[,strandcol]))
        } else if(tagPresent) {
          aln <- GRanges(seqnames = filetags[, chrcol], ir, tag = (filetags[, tagcol]), count = 1)
          if(strandPresent) strand(aln) <- Rle(as.character(filetags[,strandcol]))
          aln <- aln[order(as.factor(seqnames(aln)), as.integer(start(aln)), as.character(values(aln)$tag)),]          
          dupTags <- which(!(as.character(seqnames(aln)) == c(as.character(seqnames(aln))[-1], "!") &
                             start(aln) == c(start(aln)[-1], Inf) &
                             end(aln) == c(end(aln)[-1], Inf) &
                             as.character(strand(aln)) == c(as.character(strand(aln))[-1], "!") &
                             as.character(values(aln)$tag) == as.character(c(values(aln)$tag[-1], "!"))))
          aln <- aln[dupTags,]
          values(aln)$count <- diff(c(0, dupTags))
        } else aln <- GRanges(seqnames = filetags[, chrcol], ir)
      
      rm(filetags)
      gc()
      message(".", appendLF = FALSE)
      if(!is.null(chrs)) seqinfo(aln, new2old = match(seqlevels(seqinf), seqlevels(aln))) <- seqinf
      
      aln
    }, cols = cols, header = header)

    message(".done!")
    #if(!missing(tempFile)) save(Tags, paste("tags_", tempFile, sep = ""))

    if(is.null(chrs)) seqinf <- seqinfo(Tags[[1]])
            
    aD <- .processTags(Tags, verbose = verbose, estimationType = estimationType, seqinf = seqinf, libnames = libnames, replicates = replicates, discardTags = discardTags)
    aD
  }


.processTags <- function(GTags, estimationType, verbose = TRUE, seqinf, libnames, replicates, discardTags)
  {
    if(verbose)
      message("Analysing tags...", appendLF = FALSE)

    unqTags <- GTags[[1]]
    if(length(GTags) > 1) {
        for(ii in 2:length(GTags)) {
            message(".", appendLF = FALSE)
            unqTags <- GRanges(seqnames = c(seqnames(unqTags), seqnames(GTags[[ii]])),
                               IRanges(start = c(start(unqTags), start(GTags[[ii]])), end = c(end(unqTags), end(GTags[[ii]]))),
                               strand = c(strand(unqTags), strand(GTags[[ii]])),
                               tag = c(values(unqTags)$tag, values(GTags[[ii]])$tag))
            unqTags <- sort(unqTags[!.tagDuplicated(unqTags)])
            unqTags
        }      
    }                    

    #if(is.null(unqTags$tag)) unqTags$tag <- NA
    if(length(GTags) > 1) {
        #values(unqTags)$count <- Rle(0, length(unqTags))
        counts <- do.call("cbind",
                          lapply(GTags, function(GTag) {
                              GTag <- sort(GTag)

                              if(!discardTags) {
                                  mt <- selfmatchIntegerQuads(c(seqnames(unqTags), seqnames(GTag)),
                                                              c(start(unqTags), start(GTag)),
                                                              c(strand(unqTags), strand(GTag)),
                                                              match(c(unqTags$tag,GTag$tag), unique(unqTags$tag)))
                                  mt <- mt[-(1:length(unqTags))]
                              } else mt <- match(GTag, unqTags)
                              
                              counts <- Rle(0, length(unqTags))
                              counts[mt] <- GTag$count
                              DataFrame(counts)
                          })
                          )
        counts <- do.call("cbind", lapply(counts, as.integer))
    } else counts <- matrix(as.integer(values(GTags[[1]])$count), ncol = 1)
    
    colnames(counts) <- libnames

    if(!is.null(unqTags$tag)) {
        values(unqTags) <- values(unqTags)$tag
        colnames(values(unqTags)) <- "tag"

        ordTags <- order(as.factor(values(unqTags)$tag))
        dups <- which(!duplicated(as.factor(values(unqTags)$tag)[ordTags]))
        diffDups <- diff(c(dups, length(unqTags) + 1))
        values(unqTags)$multireads[ordTags] <- rep(diffDups, diffDups)
    } else {
        unqTags$tag <- NA
        unqTags$multireads <- 1
    }
    
    if(verbose) message(".done!")
    
#    sapply(chrs, function(x) if(any(uniqueTags$end[uniqueTags$chr == x] > chrlens[chrs == x]))
#           warning(paste("Chromosome", x, "has tags which extend over the given chromsome length.")))

    gc()
    aD <- new("alignmentData")
    aD@libnames = as.character(libnames)
#    aD@libsizes = libsizes
    aD@replicates = as.factor(replicates)
    aD@alignments = unqTags
    aD@data <- counts
        
    libSet <- !duplicated(as.character(values(aD@alignments)$tag)) | is.na(aD@alignments$tag)
    repSizes <- do.call("rbind", lapply(levels(aD@replicates), function(rep) {
        whichRep <- which(aD@replicates == rep)
        libRep <- getLibsizes(data = sapply(whichRep, function(ii) as.integer(aD@data[libSet,ii])), replicates = aD@replicates, estimationType = estimationType)
      cbind(whichRep, libRep)
    }))
    
    aD@libsizes[repSizes[,1]] <- repSizes[,2]
    
    chrmatch <- which(seqlevels(aD@alignments) %in% seqlevels(seqinf))[match(seqlevels(aD@alignments), seqlevels(seqinf))]
    chrmatch <- chrmatch[!is.na(chrmatch)]
    
    maxlens <- sapply(seqlevels(seqinf), function(x) max(end(aD@alignments[seqnames(aD@alignments) == x])))

    if(!all(is.na(seqlengths(seqinf)))) 
        if(any(seqlengths(seqinf) < maxlens)) {
            seqlengths(seqinf)[seqlengths(seqinf) < maxlens] <- maxlens[seqlengths(seqinf) < maxlens]
            warning("Some chromosomes contain tags extending over the given lengths!")
        }
    
    seqinfo(aD@alignments, new2old = chrmatch) <- seqinf

    aD <- aD[order(as.numeric(seqnames(aD@alignments)), start(aD@alignments), end(aD@alignments)),]

    aD
  }
