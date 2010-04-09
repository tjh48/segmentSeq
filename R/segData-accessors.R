setMethod("[", "segData", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      x@data <- x@data[i,, drop = FALSE]
      if(nrow(x@rightData) > 0) x@rightData <- x@rightData[i,, drop = FALSE]
      if(nrow(x@leftData) > 0) x@leftData <- x@leftData[i,, drop = FALSE]
      x@segInfo <- x@segInfo[i,,drop = FALSE]
    }

  if(!missing(j))
    {
      x@replicates <- as.integer(x@replicates[j])
      x@data <- x@data[,j,drop = FALSE]
      if(ncol(x@rightData) > 0) x@rightData <- x@rightData[,j,drop = FALSE]
      if(ncol(x@leftData) > 0) x@leftData <- x@leftData[,j,drop = FALSE]
      x@libsizes <- x@libsizes[j]
    }

  x
})

setMethod("dim", "segData", function(x) {
  dim(x@data)
})


setValidity("segData", function(object) {
  validmess <- c()
  valid <- TRUE
  if(length(object@libsizes) != ncol(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@libsizes' slot must equal number of columns of '@data' slot.")
    }
  if(nrow(object@rightData) > 0 & any(dim(object@rightData) != dim(object@data)))
    {
      valid <- FALSE
      validmess <- c(validmess, "If '@rightData' slot is non-empty, the dimensions of '@rightData' slot must equal those of the '@data' slot.")
    }
  if(nrow(object@leftData) > 0 & any(dim(object@leftData) != dim(object@data)))
    {
      valid <- FALSE
      validmess <- c(validmess, "If '@leftData' slot is non-empty, the dimensions of '@leftData' slot must equal those of the '@data' slot.")
    }
  if(nrow(object@segInfo) != nrow(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Number of rows of '@segInfo' slot not same as '@data' slot.")
    }
  if(length(object@replicates) != ncol(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@replicates' slot must equal number of columns of '@data' slot.")
    }
  if(!all(apply(object@data, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(nrow(object@rightData) > 0 & !all(apply(object@rightData, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@rightData' matrix must be castable as integers.")
    }
  if(nrow(object@leftData) > 0 & !all(apply(object@leftData, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@leftData' matrix must be castable as integers.")
    }
  if(!all(apply(object@data, c(1,2), function(x) x == as.integer(x))))
    {
      valid <- FALSE
      validmess <- c(validmess, "All members of the '@data' matrix must be castable as integers.")
    }
  if(valid) return(valid) else validmess
})


setMethod("show", "segData", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  cat('\nSlot "data":\n')
  if(nrow(object) > 5)
    {
      print(object@data[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@data)
  cat('\nSlot "libsizes":\n')
  print(object@libsizes)
  cat('\nSlot "replicates":\n')
  print(object@replicates)
  cat('\nSlot "segInfo":\n')
  if(nrow(object@segInfo) > 5)
    {
      print(object@segInfo[1:5,])
      cat(paste(nrow(object) - 5), "more rows...\n")
    } else print(object@segInfo)
  if(length(object@priorType) > 1)
    {
      cat('Slot "priors":\n')
      cat(paste('Priors are of type:', object@priorType), '\n')
    }
})


