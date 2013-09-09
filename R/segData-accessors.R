setMethod("[", "segData", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(!missing(i))
    {
      i <- as.vector(i)
      if(nrow(x@data) > 0) x@data <- x@data[i,, drop = FALSE]
    }

  if(!missing(j))
    {
      j <- as.vector(j)
      x@data <- x@data[,j,drop = FALSE]
      x@libsizes <- x@libsizes[j]
    }  
  x
})

#setMethod("dim", "segData", function(x) {
#  c(max(length(x@coordinates), nrow(x@data)), length(x@replicates))
#})


setValidity("segData", function(object) {
  validmess <- c()
  valid <- TRUE
  if(length(object@libsizes) != ncol(object))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@libsizes' slot must equal length of '@replicates' slot.")
    }
  if(nrow(object@data) != length(object@coordinates) & nrow(object@data) != 0 & length(object@coordinates) != 0)
    {
      valid <- FALSE
      validmess <- c(validmess, "Number of rows of '@data' slot (if not zero) must be the same as '@coordinates' slot (if not zero).")
    }
  if(ncol(object@data) > 0 && !all(as.integer(object@data) == object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(ncol(object@data) != length(object@replicates) & ncol(object@data) != 0 & length(object@replicates) != 0)
    {
      valid <- FALSE
      validmess <- c(validmess, "Number of columns of '@data' slot (if not zero) must be the same as the length of the '@replicates' slot (if not zero).")
    }
  if(valid) return(valid) else validmess
})


setMethod("show", "segData", function(object) {
  callNextMethod(object)
  cat('\nSlot "data":\n')
  cat("Matrix with ", nrow(object@data), " rows.")
  .printIRangesMatrix(object@data)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
})


