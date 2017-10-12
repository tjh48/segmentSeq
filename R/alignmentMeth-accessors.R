% modification on git from copied files
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
      x@nonconversion <- x@nonconversion[j]
      x@Cs <- x@Cs[,j,drop = FALSE]
      x@Ts <- x@Ts[,j,drop = FALSE]
    }
  x        
})

setValidity("alignmentMeth", function(object) {
    valid <- TRUE
    validmess <- c()
    if(length(object@replicates) != length(object@libnames))
        {
            valid <- FALSE
            validmess <- c(validmess, "The length of the '@replicates' slot must equal the length of the '@libnames' slot.")
        }  
    if(nrow(object@Cs) != length(object@alignments))
        {
            valid <- FALSE
            validmess <- c(validmess, "The number of rows in the '@Cs' slot must equal the number of rows in the '@alignments' slot.")
        }
    if(nrow(object@Ts) != length(object@alignments))
        {
            valid <- FALSE
            validmess <- c(validmess, "The number of rows in the '@Ts' slot must equal the number of rows in the '@alignments' slot.")
        }
    if(any(as.integer(object@Cs) != object@Cs))
        {
            valid <- FALSE
            validmess <- c(validmess, "All members of the '@Cs' matrix must be castable as integers.")
        }
    if(any(as.integer(object@Ts) != object@Ts))
        {
            valid <- FALSE
            validmess <- c(validmess, "All members of the '@Cs' matrix must be castable as integers.")
        }

    if(length(object@nonconversion) > 0 && any(object@nonconversion < 0) && any(object@nonconversion > 1))
        {
            valid <- FALSE
            validmess <- c(validmess, "Non-conversion rates must be between zero and one.")
        }
    if(valid) return(valid) else validmess
})
