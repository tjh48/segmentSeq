setMethod("[", "alignmentData", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(!missing(i))
    {
      i <- as.vector(i)
      x@data <- x@data[i,, drop = FALSE]
    }
  
  if(!missing(j))
    {
      j <- as.vector(j)
      x@libsizes <- x@libsizes[j]
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

setValidity("alignmentData", function(object) {
  acValid <- callNextMethod(object)
  if(class(acValid) == "character") valid <- FALSE else valid = TRUE
  if(class(acValid) == "character") valid <- acValid else valid = ""
  if(nrow(object@data) != length(object@alignments))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of rows in the '@data' slot must equal the number of rows in the '@alignments' slot.")
    }
  if(length(object@libsizes) != length(object@libnames))
    {
      valid <- FALSE
      validmess <- c(validmess, "The number of library names defined in the '@libnames' slot must equal the length of the '@libsizes' slot.")
    }
  if(any(as.integer(object@data) != object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(valid) return(valid) else validmess
})
  

setMethod("show", "alignmentData", function(object) {
  callNextMethod(object)
  cat('\nSlot "data":\n')
  .printIRangesMatrix(object@data)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
  })
