setMethod("show", "lociData", function(object) {  
  show(object@coordinates)
  callNextMethod()  
  if(nrow(object@locLikelihoods) > 0)
    {
      cat('\nSlot "locLikelihoods":\n')
      if(nrow(object@locLikelihoods) > 5)
        {
          print(exp(object@locLikelihoods[1:5,]))
          cat(paste(nrow(object) - 5), "more rows...\n")
        } else print(exp(object@locLikelihoods))
    }
})

setMethod("[", "lociData", function(x, i, j, ..., drop = FALSE) {
  if(missing(j))
    j <- 1:ncol(x@data)
  if(missing(i))
    i <- 1:nrow(x@data)
  x <- callNextMethod()
  if(nrow(x@locLikelihoods) > 0)
    x@locLikelihoods <- x@locLikelihoods[i,, drop = FALSE]
  if(length(x@coordinates) > 0)
    x@coordinates <- x@coordinates[i,, drop = FALSE]

  x
})
