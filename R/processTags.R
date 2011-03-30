processTags <-
function(files, dir = ".", replicates, libnames, chrs, chrlens,
         cols, header = TRUE, gap = 200, verbose = TRUE, ...)
  {
    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)

    chrs <- as.character(chrs)
    if(any(replicates != as.integer(replicates)))
      stop("The 'replicates' vector must be castable as an integer")
    replicates <- as.integer(replicates)
    if(any(chrlens != as.integer(chrlens)))
      stop("The 'chrlens' vector must be castable as an integer")
    chrlens <- as.integer(chrlens)    

    countPresent <- TRUE
    tagPresent <- TRUE

    if(missing(cols)) cols <- NULL
    
    if(!is.null(cols))
      {
        if(!(all(c("chr", "start", "end") %in% names(cols))))
          stop("'cols' argument must contain named values for 'chr', 'start', 'end' or be NULL")
        if(any(c(!("count" %in% names(cols)), is.na(cols[names(cols) == "count"]))))
          {
            countPresent <- FALSE
            warning("No 'count' column specified in 'cols' argument; 'processTags' will assume that the file contains non-redundant reads")
          }
        if(any(c(!("tag" %in% names(cols)), is.na(cols[names(cols) == "tag"]))))
          {
            tagPresent <- FALSE
            warning("No 'tag' column specified in 'cols' argument; the 'alignData' object will omit sequence information.")
          }
      } else if(is.null(cols) & header == FALSE) warning("No 'cols' argument supplied and 'header = FALSE'. Using default values for columns")

    if(missing(libnames))
      libnames <- sub(".*/", "", files)

    #files <- paste(dir, files, sep = "/")
    files <- file.path(dir, files)
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
              warning("No 'count' column found in file; 'processTags' will assume that the file contains non-redundant reads")
            }
        } else if(!header & is.null(cols))
          {
            chrcol <- 1L
            tagcol <- 2L
            countcol <- 3L
            startcol <- 4L
            endcol <- 5L
            tagPresent <- TRUE
            countPresent <- TRUE
          } else if(!is.null(cols)) {
            chrcol <- cols[names(cols) == "chr"]
            if(tagPresent) tagcol <- cols[names(cols) == "tag"] else tagcol <- NA          
            if(countPresent) countcol <- cols[names(cols) == "count"] else countcol <- NA
            startcol <- cols[names(cols) == "start"]
            endcol <- cols[names(cols) == "end"]
          }
      
      chrtags <- which(filetags[,chrcol] %in% chrs)

      if(tagPresent & countPresent)
        {
          aln <- data.frame(chr = I(filetags[chrtags,chrcol]), start = as.integer(filetags[chrtags,startcol]), end = as.integer(filetags[chrtags,endcol]), tag = I(filetags[chrtags, tagcol]), count = as.integer(filetags[chrtags, countcol]))
        } else if(tagPresent) {
          aln <- data.frame(chr = I(filetags[chrtags,chrcol]), start = as.integer(filetags[chrtags,startcol]), end = as.integer(filetags[chrtags,endcol]), tag = I(filetags[chrtags, tagcol]))
        } else aln <- data.frame(chr = I(filetags[chrtags,chrcol]), start = as.integer(filetags[chrtags,startcol]), end = as.integer(filetags[chrtags,endcol]))
      if(!countPresent) {
        aln <- aln[order(I(aln$chr), aln$start, aln$end),]
        alndups <- fastUniques(aln)
        aln <- data.frame(aln[alndups,], count = as.integer(diff(c(which(alndups), length(alndups) + 1))))
      }

      rm(filetags, chrtags)
      gc()
      message(".", appendLF = FALSE)
      
      aln
    }, cols = cols, header = header)

    gc()

    message(".done!")

    message("Processing files...", appendLF = FALSE)
    
    uniqueTags <- do.call("rbind", Tags)
    uniqueTags <- subset(uniqueTags, select = c(chr, start, end, tag))


    gc()
    
    if(tagPresent) {
      uniqueTags <- uniqueTags[order(as.factor(uniqueTags$tag), as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end),]
    } else uniqueTags <- uniqueTags[order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end, decreasing = TRUE),]
    uniqueTags <- uniqueTags[fastUniques(uniqueTags),]

    gc()

    if(tagPresent) {
      tagDups <- c(which(!duplicated(uniqueTags$tag)), nrow(uniqueTags) + 1)
      tagMatches <- diff(tagDups)
      uniqueTags <- data.frame(uniqueTags, matches = rep(tagMatches, tagMatches))      
    }

    message("...", appendLF = FALSE)
    
    tagCounts <- do.call("cbind", lapply(Tags, function(tagSet)
                                         {
                                           tagSet <- rbind(tagSet, data.frame(subset(uniqueTags, select = c(chr, start, end, tag)), count = 0L))
                                           if(tagPresent) {
                                             tagSet <- tagSet[order(as.factor(tagSet$tag), as.factor(tagSet$chr), tagSet$start, tagSet$end, -tagSet$count, decreasing = FALSE),]
                                           } else tagSet <- tagSet[order(as.factor(tagSet$chr), tagSet$start, tagSet$end, -tagSet$count, decreasing = FALSE),]
                                           tagSet <- tagSet[fastUniques(tagSet[,colnames(tagSet) != "count"]),]
                                           message(".", appendLF = FALSE)
                                           return(tagSet$count)
                                         }))

    if(tagPresent) libsizes <- colSums(tagCounts[!duplicated(uniqueTags$tag),]) else libsizes <- colSums(tagCounts)
    
    data <- tagCounts

    rordtags <- order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end)
    uniqueTags <- uniqueTags[rordtags,,drop = FALSE]
    data <- data[rordtags,, drop = FALSE]

    if(verbose) message(".done!")
    
    sapply(chrs, function(x) if(any(uniqueTags$end[uniqueTags$chr == x] > chrlens[chrs == x]))
           warning(paste("Chromosome", x, "has tags which extend over the given chromsome length.")))

    aD <- new("alignmentData")
    aD@libnames = libnames
    aD@libsizes = libsizes
    aD@chrs = data.frame(chr = chrs, len = chrlens)
    aD@replicates = replicates
    aD@alignments = uniqueTags
    aD@data = data 

    if(!missing(gap))
      aD <- findChunks(aD, gap)

    aD
  }

