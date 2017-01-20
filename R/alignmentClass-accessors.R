setMethod("[", "alignmentClass", function(x, i, j, ..., drop = FALSE) {
  if(!missing(i))
    {
      i <- as.vector(i)
      x@alignments <- x@alignments[i,,drop = FALSE]
    }
  
  if(!missing(j))
    {
      j <- as.vector(j)
      x@libnames <- x@libnames[j]
      x@replicates <- droplevels(x@replicates[j])      
    }
  x
})

setMethod("replicates<-", signature = "alignmentClass", function(x, value) {
  x@replicates <- as.factor(value)
  x
})

setMethod("replicates", signature = "alignmentClass", function(x) {
  x@replicates
})


.valid.alignmentClass <- function(object) {
  valid <- TRUE
  validmess <- c()
  if(length(object@replicates) != length(object@libnames))
    {
      valid <- FALSE
      validmess <- c(validmess, "The length of the '@replicates' slot must equal the length of the '@libnames' slot.")
    }
  if(valid) return(valid) else return(validmess)
}

setValidity("alignmentClass", .valid.alignmentClass)
  

setMethod("dim", "alignmentClass", function(x) {
  c(length(x@alignments), length(x@replicates))
})


setMethod("show", "alignmentClass", function(object) {
  cat(paste('An object of class "', class(object), '"\n', sep = ""))
  cat(paste(nrow(object), 'rows and', ncol(object), 'columns\n'))
  cat('\nSlot "libnames":\n')
  print(object@libnames)
  cat('\nSlot "replicates":\n')
  print(object@replicates)
  cat('\nSlot "alignments":\n')
  print(object@alignments)
})


