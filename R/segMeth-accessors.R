% modification on git from copied files
#setMethod("dim", "segMeth", function(x) {
#  nrows <- c(nrow(x@Cs), nrow(x@Ts), length(x@coordinates), nrow(x@locLikelihoods))
#  ncols <- c(length(x@replicates), ncol(x@Cs), ncol(x@Ts))
#  if(any(nrows != 0)) nrow <- nrows[nrows != 0][1] else nrow <- 0
#  if(any(ncols != 0)) ncol <- ncols[ncols != 0][1] else ncol <- 0  
#  c(nrow, ncol)
#})

#setMethod("show", "segMeth", function(object) {
#  callNextMethod(object)

#  cat('\nSlot "nonconversion":\n')
#  print(object@nonconversion)
#  cat('\nSlot "Cs":\n')
#  .printIRangesMatrix(round(object@Cs))
#  cat('\nSlot "Ts":\n')
#  .printIRangesMatrix(round(object@Ts))
#  })

#setMethod("[", "segMeth", function(x, i, j, ..., drop = FALSE) {
#  x <- callNextMethod(x, i, j, ..., drop = FALSE)
#  if(!missing(j)) {
#    j <- as.vector(j)  
#    if(nrow(x@Cs) > 0) x@Cs <- x@Cs[,j, drop = FALSE]
#    if(nrow(x@Ts) > 0) x@Ts <- x@Ts[,j, drop = FALSE]
#    if(length(x@nonconversion)) x@nonconversion <- x@nonconversion[j]
#  }
#  if(!missing(i)) {
#    i <- as.vector(i)
#    if(nrow(x@Cs) > 0) x@Cs <- x@Cs[i,, drop = FALSE]
#    if(nrow(x@Ts) > 0) x@Ts <- x@Ts[i,, drop = FALSE]
#    } 
#  x
#})
