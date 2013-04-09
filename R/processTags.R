readBAM <-
function(files, dir = ".", replicates, libnames, chrs, chrlens, countID = NULL,
         gap = 200, minlen = 15, maxlen = 29, polyLength, estimationType = "quantile", verbose = TRUE)
  {
    if(missing(polyLength)) polyLength <- NULL
    if(!is.null(polyLength)) polyBase <- do.call("rbind", lapply(c("A", "C", "G", "T"), function(polybase) {
      suppressWarnings(apply(matrix(c(".", rep(polybase, polyLength + 1)),
                                    nrow = polyLength + 1, ncol = polyLength + 1), 1, paste, sep = "", collapse = ""))
    }))
    if(length(files) != length(replicates)) stop("The 'replicates' vector needs to be the same length as the 'files' vector")
    
    chrs <- as.character(chrs)
    replicates <- as.factor(replicates)

    seqinf <- Seqinfo(seqnames = chrs, seqlengths = chrlens)

    if(!all(levels(replicates) %in% replicates))
      stop("There appear to be additional levels in your (factor) replicates specification which are not present in the replicates vector.")
    
    if(any(chrlens != as.integer(chrlens)))
      stop("The 'chrlens' vector must be castable as an integer")
    chrlens <- as.integer(chrlens)    

    if(missing(libnames))
      libnames <- sub(".*/", "", files)

    files <- paste(dir, files, sep = "/")

    if(class(chrs) != "character")
      stop("'chrs' must be of type 'character'.")

    if(verbose)
      message("Reading files...", appendLF = FALSE)

    sampleNumbers <- 1:length(files)

    Tags <- lapply(sampleNumbers, function(ii) {

      tags <- scanBam(files[ii])[[1]]
      if(!is.null(countID)) {
        counts <- Rle(scanBam(files[ii], param = ScanBamParam(tag=countID))[[1]][[1]][[1]])
      } else counts <- Rle(rep(1L, length(tags$seq)))

      revcomp <- sapply(tags$flag, function(flag) intToBits(flag)[5] == 1)
      if(any(revcomp)) tags$seq[revcomp] <- reverseComplement(tags$seq[revcomp])

      keepReads <- which(!is.na(tags$pos) & tags$qwidth > minlen & tags$qwidth < maxlen)
      
      ir <- IRanges(start = as.integer(tags$pos[keepReads]), width = tags$qwidth[keepReads])
      aln <- GRanges(seqnames = tags$rname[keepReads], ir, strand = tags$strand[keepReads], tag = Rle(as.character(tags$seq[keepReads])), count = counts[keepReads])
      aln <- aln[as.integer(seqnames(aln)) %in% which(seqlevels(aln) %in% chrs),]

      if(!is.null(polyLength)) {
        polyBaseRemove <- unique(unlist(lapply(polyBase, grep, as.character(values(aln)$tag))))
        if(length(polyBaseRemove) > 0) aln <- aln[-polyBaseRemove,]
      }
        if(is.null(countID)) {
          aln <- aln[order(as.factor(seqnames(aln)), as.integer(start(aln)), as.character(values(aln)$tag)),]
          dupTags <- which(!(as.character(seqnames(aln)) == c(as.character(seqnames(aln))[-1], "!") &
                             start(aln) == c(start(aln)[-1], Inf) &
                             end(aln) == c(end(aln)[-1], Inf) &
                             as.character(strand(aln)) == c(as.character(strand(aln))[-1], "!") &
                             as.character(values(aln)$tag) == as.character(c(values(aln)$tag[-1], "!"))))
          aln <- aln[dupTags,]
          values(aln)$count <- diff(c(0, dupTags))
      }

      message(".", appendLF = FALSE)
      seqinfo(aln, new2old = match(seqlevels(seqinf), seqlevels(aln))) <- seqinf
      aln
      })

    message(".done!")
    
    .processTags(Tags, verbose = verbose, gap = gap, estimationType = estimationType, seqinf = seqinf, libnames = libnames, replicates = replicates)
    
  }


readGeneric <-
function(files, dir = ".", replicates, libnames, chrs, chrlens,
         cols, header = TRUE, gap = 200, minlen = 15, maxlen = 29, polyLength, estimationType = "quantile", verbose = TRUE, ...)
  {
    if(missing(polyLength)) polyLength <- NULL
    if(!is.null(polyLength)) polyBase <- do.call("rbind", lapply(c("A", "C", "G", "T"), function(polybase) {
      suppressWarnings(apply(matrix(c(".", rep(polybase, polyLength + 1)),
                                    nrow = polyLength + 1, ncol = polyLength + 1), 1, paste, sep = "", collapse = ""))
    }))
    
    chrs <- as.character(chrs)
                                        #    if(any(replicates != as.integer(replicates)))
                                        #      stop("The 'replicates' vector must be castable as an integer")
    replicates <- as.factor(replicates)

    seqinf <- Seqinfo(seqnames = chrs, seqlengths = chrlens)
    
    if(!all(levels(replicates) %in% replicates))
      stop("There appear to be additional levels in your (factor) replicates specification which are not present in the replicates vector.")

    
    if(any(chrlens != as.integer(chrlens)))
      stop("The 'chrlens' vector must be castable as an integer")
    chrlens <- as.integer(chrlens)    

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

    if(class(chrs) != "character")
      stop("'chrs' must be of type 'character'.")

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
                 
      if(tagPresent & !is.null(polyLength)) {        
        polyBaseRemove <- unique(unlist(lapply(polyBase, grep, filetags[,tagcol])))
        if(length(polyBaseRemove) > 0) filetags <- filetags[-polyBaseRemove,]
      }

      chrtags <- which(filetags[,chrcol] %in% chrs)
      
      ir <- IRanges(start = as.integer(filetags[chrtags,startcol]), end = as.integer(filetags[chrtags, endcol]))      
      
      if(tagPresent & countPresent)
        {
          aln <- GRanges(seqnames = filetags[chrtags, chrcol], ir, tag = Rle(filetags[chrtags, tagcol]), count = Rle(as.integer(filetags[chrtags, countcol])))
          if(strandPresent) strand(aln) <- Rle(as.character(filetags[chrtags,strandcol]))
        } else if(tagPresent) {
          aln <- GRanges(seqnames = filetags[chrtags, chrcol], ir, tag = Rle(filetags[chrtags, tagcol]), count = 1)
          if(strandPresent) strand(aln) <- Rle(as.character(filetags[chrtags,strandcol]))
          aln <- aln[order(as.factor(seqnames(aln)), as.integer(start(aln)), as.character(values(aln)$tag)),]          
          dupTags <- which(!(as.character(seqnames(aln)) == c(as.character(seqnames(aln))[-1], "!") &
                             start(aln) == c(start(aln)[-1], Inf) &
                             end(aln) == c(end(aln)[-1], Inf) &
                             as.character(strand(aln)) == c(as.character(strand(aln))[-1], "!") &
                             as.character(values(aln)$tag) == as.character(c(values(aln)$tag[-1], "!"))))
          aln <- aln[dupTags,]
          values(aln)$count <- diff(c(0, dupTags))
        } else aln <- GRanges(seqnames = filetags[chrtags, chrcol], ir)

      aln <- aln[width(aln) >= minlen & width(aln) <= maxlen]
      
      rm(filetags, chrtags)
      gc()
      message(".", appendLF = FALSE)
      seqinfo(aln, new2old = match(seqlevels(seqinf), seqlevels(aln))) <- seqinf
      aln
    }, cols = cols, header = header)

    message(".done!")
    
    .processTags(Tags, verbose = verbose, estimationType = estimationType, gap = gap, seqinf = seqinf, libnames = libnames, replicates = replicates)
  }


.processTags <- function(GTags, gap, estimationType, verbose = TRUE, seqinf, libnames, replicates)
  {
    if(verbose)
      message("Analysing tags...", appendLF = FALSE)

    unqTags <- GTags[[1]]
    if(length(GTags) > 1) {
      for(ii in 2:length(GTags)) {
        unqTags <- GRanges(seqnames = c(seqnames(unqTags), seqnames(GTags[[ii]])),
                           IRanges(start = c(start(unqTags), start(GTags[[ii]])), end = c(end(unqTags), end(GTags[[ii]]))),
                           strand = c(strand(unqTags), strand(GTags[[ii]])),
                           tag = c(values(unqTags)$tag, values(GTags[[ii]])$tag))
        ordTag <- order(as.integer(seqnames(unqTags)), as.integer(start(unqTags)), as.integer(end(unqTags)), as.character(values(unqTags)$tag), as.character(strand(unqTags)))

        purgeTags <- which(as.integer(seqnames(unqTags)[ordTag]) == c(as.integer(seqnames(unqTags)[ordTag[-1]]), -1) &
                           start(unqTags)[ordTag] == c(start(unqTags)[ordTag[-1]], Inf) &
                           end(unqTags)[ordTag] == c(end(unqTags)[ordTag[-1]], Inf) &
                           as.character(strand(unqTags)[ordTag]) == c(as.character(strand(unqTags)[ordTag[-1]]), "!") &
                           (values(unqTags)$tag[ordTag] == c(values(unqTags)$tag[ordTag[-1]], "!") | (is.na(values(unqTags)$tag[ordTag]) & is.na(c(values(unqTags)$tag[ordTag[-1]], NA))))) + 1
        if(length(purgeTags) > 0)
          unqTags <- unqTags[ordTag[-purgeTags],]
      }      
    }                    

    if(length(GTags) > 1) {
      values(unqTags)$count <- 0
      counts <- do.call("DataFrame",
                        lapply(GTags, function(GTag) {
                          
                          countTag <- c(GTag, unqTags)
                          countTag <- countTag[order(as.integer(seqnames(countTag)), as.integer(start(countTag)), as.integer(end(countTag)), as.character(values(countTag)$tag), as.character(strand(countTag))),]
                          
                          purgeTags <- which(as.integer(seqnames(countTag)) == c(as.integer(seqnames(countTag)[-1]), -1) &
                                             start(countTag) == c(start(countTag)[-1], Inf) &
                                             end(countTag) == c(end(countTag)[-1], Inf) &
                                             as.character(strand(countTag)) == c(as.character(strand(countTag)[-1]), "!") &
                                             (values(countTag)$tag == c(values(countTag)$tag[-1], "!") | (is.na(values(countTag)$tag) & is.na(c(values(countTag)$tag[-1], NA))))) + 1
                                  message(".", appendLF = FALSE)
                                  if(length(purgeTags) > 0) return(values(countTag)$count[-purgeTags]) else return(values(countTag)$count)
                        })
                        )
    } else counts <- DataFrame(values(GTags[[1]])$count)
    
    colnames(counts) <- libnames

    values(unqTags) <- values(unqTags)$tag
    colnames(values(unqTags)) <- "tag"

    ordTags <- order(as.factor(values(unqTags)$tag))
    dups <- which(!duplicated(as.factor(values(unqTags)$tag)[ordTags]))
    diffDups <- diff(c(dups, length(unqTags) + 1))
    values(unqTags)$matches[ordTags] <- rep(diffDups, diffDups)
    
    if(verbose) message(".done!")
    
#    sapply(chrs, function(x) if(any(uniqueTags$end[uniqueTags$chr == x] > chrlens[chrs == x]))
#           warning(paste("Chromosome", x, "has tags which extend over the given chromsome length.")))

    aD <- new("alignmentData")
    aD@libnames = libnames
#    aD@libsizes = libsizes
    aD@replicates = as.factor(replicates)
    aD@alignments = unqTags
    aD@data = counts    

    libSet <- !duplicated(as.character(values(aD@alignments)$tag))
    repSizes <- do.call("rbind", lapply(levels(aD@replicates), function(rep) {
      whichRep <- which(aD@replicates == rep)
      libRep <- getLibsizes(data = sapply(whichRep, function(ii) as.integer(aD@data[libSet,ii])), replicates = aD@replicates, estimationType = estimationType)
      cbind(whichRep, libRep)
    }))
    
    aD@libsizes[repSizes[,1]] <- repSizes[,2]
    
    if(!missing(gap))
      aD@alignments <- findChunks(aD@alignments, gap)
    
    chrmatch <- which(seqlevels(aD@alignments) %in% seqlevels(seqinf))[match(seqlevels(aD@alignments), seqlevels(seqinf))]
    chrmatch <- chrmatch[!is.na(chrmatch)]
    
    maxlens <- sapply(seqlevels(seqinf), function(x) max(end(aD@alignments[seqnames(aD@alignments) == x])))

    if(any(seqlengths(seqinf) < maxlens)) {
      seqlengths(seqinf)[seqlengths(seqinf) < maxlens] <- maxlens[seqlengths(seqinf) < maxlens]
      warning("Some chromosomes contain tags extending over the given lengths!")
    }
    
    seqinfo(aD@alignments, new2old = chrmatch) <- seqinf
    
    aD
  }
