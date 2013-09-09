setMethod("[", "segClass", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      i <- as.vector(i)
      x@coordinates <- x@coordinates[i,,drop = FALSE]
      if(nrow(x@locLikelihoods) > 0)
        x@locLikelihoods <- x@locLikelihoods[i,,drop = FALSE]
    }

  if(!missing(j))
    {
      j <- as.vector(j)
      x@replicates <- x@replicates[j]
    }  
  x
})

setMethod("dim", "segClass", function(x) {
  c(max(length(x@coordinates), nrow(x@locLikelihoods)), length(x@replicates))
})


setMethod("show", "segClass", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  cat('\nSlot "replicates":\n')
  print(object@replicates)
  cat('\nSlot "coordinates":\n')
  print(object@coordinates)
  .printLocLikes(object@locLikelihoods)
})


