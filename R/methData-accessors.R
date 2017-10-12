% modification on git from copied files

setMethod("[", "methData", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(missing(j))
    j <- 1:ncol(x)
  if(missing(i))
    i <- 1:nrow(x)

  i <- as.vector(i)
  j <- as.vector(j)
  
  if(length(x@coordinates) > 0)
    x@coordinates <- x@coordinates[i,, drop = FALSE]
  if(nrow(x@locLikelihoods) > 0)
    x@locLikelihoods <- x@locLikelihoods[i,, drop = FALSE]
  x
})

setMethod("show", "methData", function(object) {  
  cat('\nSlot "coordinates"\n')
  show(object@coordinates)
  callNextMethod()
  .printLocLikes(object@locLikelihoods)
})
