setMethod("[", "segData", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      x@data <- x@data[i,, drop = FALSE]
      x@coordinates <- x@coordinates[i,,drop = FALSE]
#      if(nrow(x@seglens) > 0)
#        x@seglens <- x@seglens[i,,drop = FALSE]
      if(nrow(x@locLikelihoods) > 0)
        x@locLikelihoods <- x@locLikelihoods[i,,drop = FALSE]
    }

  if(!missing(j))
    {
      x@replicates <- x@replicates[j]
      x@data <- x@data[,j,drop = FALSE]
      x@libsizes <- x@libsizes[j]
#      if(ncol(x@seglens) > 1) x@seglens <- x@seglens[,j,drop = FALSE]
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
  if(length(object@coordinates) != nrow(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Number of rows of '@coordinates' slot not same as '@data' slot.")
    }
  if(length(object@replicates) != ncol(object@data))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@replicates' slot must equal number of columns of '@data' slot.")
    }
  if(!all(sapply(1:ncol(object), function(ii) all(as.integer(object@data[,ii]) == object@data[,ii]))))
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
  cat('\nSlot "coordinates":\n')
  if(length(object@coordinates) > 5)
    {
      print(object@coordinates)
    } else print(object@coordinates)
})


