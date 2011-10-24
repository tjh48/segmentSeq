setClass("segData", representation(data = "DataFrame", libsizes = "numeric", replicates = "factor", coordinates = "GRanges"))#, seglens = "matrix"))

setClass("alignmentData", representation(alignments = "GRanges", data = "DataFrame", libnames = "character", libsizes = "numeric", replicates = "factor"))

setClass("lociData", representation(locLikelihoods = "matrix", coordinates = "GRanges"), contains = "countData")
