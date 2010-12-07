processTags <-
function(files, dir = ".", replicates, libnames, chrs, chrlens,
         cols, header = TRUE, verbose = TRUE, ...)
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

    files <- paste(dir, files, sep = "/")

    if(class(chrs) != "character")
      stop("'chrs' must be of type 'character'.")

    if(verbose)
      message("processing...", appendLF = FALSE)

    sampleNumbers <- 1:length(files)
    
    Tags <- lapply(sampleNumbers, function(ii, cols, header, ...) {
      filetags <- read.table(files[ii], header = header, as.is = TRUE, ...)
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
      
      aln <- data.frame(chr = I(filetags[chrtags,chrcol]), start = as.integer(filetags[chrtags,startcol]), end = as.integer(filetags[chrtags,endcol]))
      if(tagPresent) aln <- data.frame(aln, tag = I(filetags[chrtags, tagcol]))
      if(!countPresent) {
        aln <- aln[order(as.factor(aln$chr), aln$start, aln$end),]
        alndups <- fastUniques(aln)
        aln <- data.frame(aln[alndups,], count = as.integer(diff(c(which(alndups), length(alndups) + 1))))
      } else aln <- data.frame(aln, count = as.integer(filetags[chrtags, countcol]))
      
      aln
    }, cols = cols, header = header)
    
    uniqueTags <- NULL
    for(ii in 1:length(Tags))
      {
        uniqueTags <- rbind(uniqueTags, Tags[[ii]][,which(colnames(Tags[[ii]]) != "count"), drop = FALSE])
        if(tagPresent) {
          uniqueTags <- uniqueTags[order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end, as.factor(uniqueTags$tag), decreasing = TRUE),]
        } else uniqueTags <- uniqueTags[order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end, decreasing = TRUE),]
        
        uniqueTags <- uniqueTags[fastUniques(uniqueTags),]
        message(".", appendLF = FALSE)
      }
    
    Tags <- lapply(Tags, function(tagSet)
                   {
                     tagSet <- rbind(tagSet, data.frame(uniqueTags, count = 0))
                     if(tagPresent) {
                       tagSet <- tagSet[order(as.factor(tagSet$chr), tagSet$start, tagSet$end, as.factor(tagSet$tag), tagSet$count, decreasing = TRUE),]
                     } else tagSet <- tagSet[order(as.factor(tagSet$chr), tagSet$start, tagSet$end, tagSet$count, decreasing = TRUE),]
                     tagSet <- tagSet[fastUniques(tagSet[,colnames(tagSet) != "count"]),]
                     tagSet
                   })
    
    libsizes <- sapply(Tags, function(tagSet) if(tagPresent) return(sum(tagSet$count[!duplicated(tagSet$tag)])) else return(sum(tagSet$count)))
    
    data <- matrix(unlist(lapply(Tags, function(x) x$count)), ncol = length(sampleNumbers), nrow = nrow(uniqueTags))

    rordtags <- order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end)
    uniqueTags <- uniqueTags[rordtags,,drop = FALSE]
    data <- data[rordtags,, drop = FALSE]

    if(tagPresent) alignments <- data.frame(uniqueTags, duplicated = uniqueTags$tag %in% uniqueTags$tag[duplicated(uniqueTags$tag)]) else alignments <- data.frame(uniqueTags, duplicated = FALSE)
    
    sapply(chrs, function(x) if(any(alignments$end[alignments$chr == x] > chrlens[chrs == x]))
           warning(paste("Chromosome", x, "has tags which extend over the given chromsome length.")))
    
    if(verbose) message(".done!")
    
    aD <- new("alignmentData", libnames = libnames, libsizes = libsizes, alignments = alignments, data = data, chrs = chrs, chrlens = chrlens, replicates = replicates)
    aD
  }

