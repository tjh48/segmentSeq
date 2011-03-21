setClass("segData", representation(data = "matrix", rightData = "matrix", leftData = "matrix", libsizes = "numeric", chrs = "data.frame", replicates = "integer", priorType = "character", priors = "list", segInfo = "data.frame"))

setClass("alignmentData", representation(alignments = "data.frame", data = "matrix", libnames = "character", libsizes = "numeric", chrs = "data.frame", replicates = "integer"))

setClass("postSeg", contains = "countData")
