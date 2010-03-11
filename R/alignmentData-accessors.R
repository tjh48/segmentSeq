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
      x@replicates <- x@replicates[j]
      x@data <- x@data[,j,drop = FALSE]
    }
  x
})


setValidity("alignmentData", function(object) {
  valid <- TRUE
  validmess <- c()
  if(!is.character(object@chrs))
    {
      valid = FALSE
      validmess <- c(validmess, "The '@chrs' slot must have class 'character'")
    }
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
  if(length(object@chrs) != length(object@chrlens))
    {
      valid <- FALSE
      validmess <- c(validmess, "The length of the '@chrs' slot must equal the length of the '@chrlens' slot.")
    }
  if(!all(unique(object@alignments$chr) %in% object@chrs))
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
  cat('\nSlot "chrlens":\n')
  print(object@chrlens)
  })


