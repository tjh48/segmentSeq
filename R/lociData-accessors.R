setMethod("show", "lociData", function(object) {  
  show(object@coordinates)
  callNextMethod()  
  .printLocLikes(object@locLikelihoods)  
  if(nrow(object@locLikelihoods) > 0) {
    cat("\nExpected number of loci in each replicate group\n")
    print(colSums(exp(object@locLikelihoods)))
  }
})

setMethod("[", "lociData", function(x, i, j, ..., drop = FALSE) {
  if(missing(j))
    j <- 1:ncol(x@data)
  if(missing(i))
    i <- 1:nrow(x@data)
  

  if(length(i) == 0) return(x)
  
  i <- as.vector(i)
  j <- as.vector(j)
  
  x <- callNextMethod()
  if(nrow(x@locLikelihoods) > 0)
    x@locLikelihoods <- x@locLikelihoods[i,, drop = FALSE]
  if(length(x@coordinates) > 0)
    x@coordinates <- x@coordinates[i,, drop = FALSE]

  x
})

setMethod("dim", "lociData", function(x) {
    dim <- dim(x@data)
    if(dim[1] == 0) dim[1] <- length(x@coordinates)
    dim
})
