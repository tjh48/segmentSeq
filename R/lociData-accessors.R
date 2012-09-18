setMethod("show", "lociData", function(object) {  
  show(object@coordinates)
  callNextMethod()  
  if(nrow(object@locLikelihoods) > 0)
    {
      if(any(exp(object@locLikelihoods) > 1, na.rm = TRUE))
        {
          cat('\nSlot "locLikelihoods":\n')
          modFunction <- identity
        } else {          
          cat('\nSlot "locLikelihoods" (stored on log scale):\n')        
          modFunction <- exp
        }
      if(nrow(object@locLikelihoods) > 5)
        {          
          print(modFunction(object@locLikelihoods[1:5,]))
          cat(paste(nrow(object) - 5), "more rows...\n")
        } else print(modFunction(object@locLikelihoods))
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
