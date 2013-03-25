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

setMethod("libsizes<-", signature = "alignmentData", function(x, value) {
  if(!is.numeric(value)) stop("All members of libsizes for an alignmentData object must be numeric.")
  if(length(value) != ncol(x)) stop("Length of libsizes must be identical to the number of columns of the alignmentData object.")
  if(any(value <= 0)) stop("Library sizes less than or equal to zero make no sense to me!")
  x@libsizes <- value
  x
})

setMethod("libsizes", signature = "alignmentData", function(x) {
  x@libsizes
})

setMethod("replicates<-", signature = "alignmentData", function(x, value) {
  x@replicates <- as.factor(value)
  x
})

setMethod("replicates", signature = "alignmentData", function(x) {
  x@replicates
})


setValidity("alignmentData", function(object) {
  valid <- TRUE
  validmess <- c()
  if(nrow(object@data) != length(object@alignments))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of rows in the '@data' slot must equal the number of rows in the '@alignments' slot.")
    }
  if(ncol(object@data) != length(object@libnames))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of columns in the '@data' slot must equal the length of the '@libnames' slot.")
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
  if(!all(sapply(1:ncol(object), function(ii) all(as.integer(object@data[,ii]) == object@data[,ii]))))
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
  print(object@alignments)
  cat('\nSlot "data":\n')
  if(nrow(object) > 5)
    {
      print(object@data[1:5,])
      cat(paste(nrow(object@data) - 5), "more rows...\n")
    } else print(object@data)  
  cat('\nSlot "libnames":\n')
  print(object@libnames)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
  cat('\nSlot "replicates":\n')
  print(object@replicates)
  })


