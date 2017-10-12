% modification on git from copied files
setClass("alignmentClass", representation(alignments = "GRanges", libnames = "character", replicates = "factor"))
setClass("alignmentData", representation(data = "matrix", libsizes = "numeric"), contains = "alignmentClass")
setClass("alignmentMeth", representation(Cs = "matrix", Ts = "matrix", nonconversion = "numeric"), contains = "alignmentClass")


#setClass("segClass", representation(coordinates = "GRanges", locLikelihoods = "DataFrame", replicates = "factor"))
#setClass("segData", representation(data = "matrix", libsizes = "numeric"), contains = "segClass")
#setClass("segMeth", representation(Cs = "matrix", Ts = "matrix", nonconversion = "numeric"), contains = "segClass")

#setClass("lociData", contains = "countData")

setClass("lociData", representation(locLikelihoods = "matrix", coordinates = "GRanges"), contains = "countData")
setClass("methData", representation(locLikelihoods = "matrix", coordinates = "GRanges"), contains = "countData")
#setClass("methData", representation(locLikelihoods = "matrix", coordinates = "GRanges"), contains = "pairedData")
