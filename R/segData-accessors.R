setMethod("[", "segData", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      if(nrow(x@data) > 0) x@data <- x@data[i,, drop = FALSE]
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
  c(max(length(x@coordinates), nrow(x@data)), length(x@replicates))
})


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
  if(length(object@replicates) != ncol(object))
    {
      valid <- FALSE
      validmess <- c(validmess, "Length of '@replicates' slot must equal length of '@replicates' slot.")
    }
  if(ncol(object@data) > 0 && !all(sapply(1:ncol(object@data), function(ii) all(as.integer(object@data[,ii]) == object@data[,ii]))))
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
  if(nrow(object@data) > 5)
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


