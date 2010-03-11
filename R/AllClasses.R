setClass("segData", representation(data = "matrix", libsizes = "numeric", replicates = "integer", priorType = "character", priors = "list", segInfo = "data.frame"))

setClass("alignmentData", representation(alignments = "data.frame", data = "matrix", libnames = "character", libsizes = "numeric", chrs = "character", chrlens = "integer", replicates = "integer"))
