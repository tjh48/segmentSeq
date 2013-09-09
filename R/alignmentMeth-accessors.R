setMethod("show", "alignmentMeth", function(object) {
  callNextMethod(object)

  cat('\nSlot "Cs":\n')
  .printIRangesMatrix(object@Cs)
  cat('\nSlot "Ts":\n')
  .printIRangesMatrix(object@Ts)
  
  cat('\nSlot "nonconversion":\n')
  print(object@nonconversion)
  })

setMethod("[", "alignmentMeth", function(x, i, j, ..., drop = FALSE) {
  x <- callNextMethod(x, i, j, ..., drop = FALSE)
  if(!missing(i)) {
    i <- as.vector(i)
    x@Cs <- x@Cs[i,, drop = FALSE]
    x@Ts <- x@Ts[i,, drop = FALSE]
  }
  
  if(!missing(j))
    {
      j <- as.vector(j)
      x@nonconversion <- x@libsizes[j]
      x@Cs <- x@Cs[,j,drop = FALSE]
      x@Ts <- x@Ts[,j,drop = FALSE]
    }
  x        
})

setValidity("alignmentData", function(object) {
  acValid <- callNextMethod(object)
  if(class(acValid) == "character") valid <- FALSE else valid = TRUE
  if(class(acValid) == "character") valid <- acValid else valid = ""
  if(length(object@nonconversion) > 0 && any(object@nonconversion < 0) && any(object@nonconversion > 1))
    {
      valid <- FALSE
      validmess <- c(validmess, "Non-conversion rates must be between zero and one.")
    }
  if(valid) return(valid) else validmess
})
