setGeneric("cbind", function(..., deparse.level=1) standardGeneric("cbind"), signature = "...")

setMethod("cbind", "alignmentData", function(..., deparse.level = 1) {
    binds <- list(...)

    fastUniques <- function(x)
      if(nrow(x) > 1) return(c(TRUE, rowSums(x[-1L,, drop = FALSE] == x[-nrow(x),,drop = FALSE]) != ncol(x))) else return(TRUE)

    uniqueTags <- do.call("rbind", lapply(binds, function(x) x@alignments))
    uniqueTags <- uniqueTags[order(as.factor(uniqueTags$chr), uniqueTags$start, uniqueTags$end, as.factor(uniqueTags$tag), decreasing = TRUE),]
    uniqueTags <- subset(uniqueTags, select = c("chr", "start", "end", "tag"))
    uniqueTags <- uniqueTags[fastUniques(uniqueTags),]

    mergeUniques <- function(z, uniqueTags)
      {
        if("tag" %in% colnames(z@alignments)) uniqueSelect <- !(uniqueTags$tag %in% z@alignments$tag) else uniqueSelect <- rep(TRUE, nrow(uniqueTags))
        zalign <- rbind(subset(z@alignments, select = c("chr", "start", "end", "tag")), uniqueTags[uniqueSelect,])
        zdata <- rbind(z@data, matrix(0, nrow = sum(uniqueSelect), ncol = ncol(z)))
        ordAlign <- with(zalign, order(as.factor(chr), start, end, as.factor(tag), decreasing = TRUE))  
        zalign <- zalign[ordAlign,, drop = FALSE]
        zdata <- zdata[ordAlign,,drop = FALSE]
        zun <- fastUniques(zalign)
        zdata <- zdata[zun,,drop = FALSE]
        zalign <- zalign[zun,,drop = FALSE]
      
        list(data = zdata, alignments = zalign)
      }

    mergedData <- lapply(binds, mergeUniques, uniqueTags = uniqueTags)

    chrs <- do.call("rbind", lapply(binds, function(x) x@chrs))
    chrs <- chrs[!duplicated(chrs),]
    
    naD <- new("alignmentData", data = do.call("cbind", lapply(mergedData, function(x) x$data)),
               alignments = mergedData[[1]]$alignments,
               libsizes = unlist(lapply(binds, function(x) x@libsizes)),
               libnames = unlist(lapply(binds, function(x) x@libnames)),
               replicates = as.integer(rep(1, sum(sapply(binds, ncol)))),
               chrs = chrs)

    naD <- naD[with(naD@alignments, order(as.factor(chr), start, end)),]
    naD@alignments <- data.frame(naD@alignments, duplicated = with(naD@alignments, tag %in% tag[duplicated(tag)]))
  
  naD

  })


setMethod("[", "alignmentData", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      x@alignments <- x@alignments[i,,drop = FALSE]
      x@data <- x@data[i,, drop = FALSE]
    }
  
  if(!missing(j))
    {
      x@libnames <- x@libnames[j]
      x@libsizes <- x@libsizes[j]
      x@replicates <- as.integer(x@replicates[j])
      x@data <- x@data[,j,drop = FALSE]
    }
  x
})


setValidity("alignmentData", function(object) {
  valid <- TRUE
  validmess <- c()
  if(nrow(object@data) != nrow(object@alignments))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of the rows in the '@data' slot must equal the number of rows in the '@alignments' slot.")
    }
  if(ncol(object@data) != length(object@libnames))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of the rows in the '@data' slot must equal the length of the '@libnames' slot.")
    }
  if(length(object@libnames) != length(object@libsizes))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of library names defined in the '@libnames' slot must equal the length of the '@libsizes' slot.")
    }
  if(length(object@replicates) != length(object@libnames))
    {
      valid <- FALSE
      validmess <- c(validmess, "The length of the '@replicates' slot must equal the length of the '@libnames' slot.")
    }
  if(!all(unique(object@alignments$chr) %in% object@chrs$chr))
    {
      valid <- FALSE
      validmess <- c(validmess, "All the chromosomes defined in the '$chr' column of the '@alignments' slot must be described in the '@chrs' slot.")
    }
  if(!all(apply(object@data, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(valid) return(valid) else validmess
})
  

setMethod("dim", "alignmentData", function(x) {
  dim(x@data)
})



setMethod("show", "alignmentData", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  cat('\nSlot "alignments":\n')
  if(nrow(object) > 5)
    {
      print(object@alignments[1:5,])
      cat(paste(nrow(object@alignments) - 5), "more rows...\n")
    } else print(object@alignments)
  cat('\nSlot "data":\n')
  if(nrow(object) > 5)
    {
      print(object@data[1:5,])
      cat(paste(nrow(data) - 5), "more rows...\n")
    } else print(object@data)  
  cat('\nSlot "libnames":\n')
  print(object@libnames)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
  cat('\nSlot "replicates":\n')
  print(object@replicates)
  cat('\nSlot "chrs":\n')
  print(object@chrs)
  })


