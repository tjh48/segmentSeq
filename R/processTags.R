processTags <-
function(files, replicates, libnames, chrs, chrlens,
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

    
    
    if(!missing(cols))
      if(!(all(c("chr", "tag", "count", "start", "end") %in% names(cols))))
        stop("'cols' argument must contain named values for 'chr', 'tag', 'count', 'start', 'end' or be NULL")
        
    if(missing(cols))
      {
        if(!header)
          warning("No 'cols' argument supplied and 'header = FALSE'. Using default values for columns")
        chrcol <- 1L
        tagcol <- 2L
        countcol <- 3L
        startcol <- 4L
        endcol <- 5L
      } else {
        chrcol <- cols[names(cols) == "chr"]
        tagcol <- cols[names(cols) == "tag"]
        countcol <- cols[names(cols) == "count"]
        startcol <- cols[names(cols) == "start"]
        endcol <- cols[names(cols) == "end"]
      }

    if(is.null(libnames) | missing(libnames))
      libnames <- sub(".*/", "", files)

    if(class(chrs) != "character")
      stop("'chrs' must be of type 'character'.")

    if(verbose)
      message("processing...", appendLF = FALSE)

    sampleNumbers <- 1:length(files)
    
    Tags <- lapply(sampleNumbers, function(ii, cols, header, ...) {
      filetags <- read.table(files[ii], header = header, as.is = TRUE, ...)
      if(header & missing(cols))
        if(all(c("chr", "tag", "count", "start", "end") %in% names(filetags)))
          {
            chrcol <- which(names(filetags) == "chr")
            tagcol <- which(names(filetags) == "tag")
            countcol <- which(names(filetags) == "count")
            startcol <- which(names(filetags) == "start")
            endcol <- which(names(filetags) == "end")
          } else stop(paste("Couldn't find appropriate column names (and columns were not specified) in file:", files[ii]))
      
      chrtags <- which(filetags[,chrcol] %in% chrs)
      
      data.frame(chr = I(filetags[chrtags,chrcol]), start = as.integer(filetags[chrtags,startcol]), end = as.integer(filetags[chrtags,endcol]),
                         tag = I(filetags[chrtags,tagcol]), count = as.integer(filetags[chrtags,countcol]))
    }, cols = cols, header = header, ...)

    uniqueTags <- NULL
    for(ii in 1:length(Tags))
      {
        uniqueTags <- rbind(uniqueTags, Tags[[ii]][,which(colnames(Tags[[ii]]) != "count"), drop = FALSE])
        uniqueTags <- uniqueTags[order(as.factor(uniqueTags$tag), as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end, decreasing = TRUE),]
        uniqueTags <- uniqueTags[fastUniques(uniqueTags),]
        message(".", appendLF = FALSE)
      }

    Tags <- lapply(Tags, function(tagSet)
                   {
                     tagSet <- rbind(tagSet, data.frame(uniqueTags, count = 0))
                     tagSet <- tagSet[order(as.factor(tagSet$tag), as.factor(tagSet$chr), tagSet$start, tagSet$end, tagSet$count, decreasing = TRUE),]
                     tagSet <- tagSet[fastUniques(tagSet[colnames(tagSet) != "count"]),]
                     tagSet
                   })

    libsizes <- sapply(Tags, function(tagSet) sum(tagSet$count[!duplicated(tagSet$tag)]))
    
    data <- matrix(unlist(lapply(Tags, function(x) x$count)), ncol = length(sampleNumbers), nrow = nrow(uniqueTags))
    rordtags <- order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end)

    data <- data[rordtags,]
    alignments <- data.frame(uniqueTags, duplicated = uniqueTags$tag %in% uniqueTags$tag[duplicated(uniqueTags$tag)])[rordtags,]
    

    sapply(chrs, function(x) if(any(alignments$end[alignments$chr == x] > chrlens[chrs == x]))
           warning(paste("Chromosome", x, "has tags which extend over the given chromsome length.")))
    
    if(verbose) message(".done!")
    
    aD <- new("alignmentData", libnames = libnames, libsizes = libsizes, alignments = alignments, data = data, chrs = chrs, chrlens = chrlens, replicates = replicates)
    aD
  }

