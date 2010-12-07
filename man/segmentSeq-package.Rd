\name{segmentSeq-package}
\alias{segmentSeq-package}
\alias{segmentSeq}
\docType{package}
\title{
Segmentation of the genome based on multiple samples of high-throughput
sequencing data.
}
\description{
The segmentSeq package is intended to take multiple samples of
high-throughput data (together with replicate information) and identify
regions of the genome which have a (reproducibly) high density of tags
aligning to them.
}
\details{
\tabular{ll}{
Package: \tab segmentSeq\cr
Type: \tab Package\cr
Version: \tab 0.0.2\cr
Date: \tab 2010-01-20\cr
License: \tab GPL-3 \cr
LazyLoad: \tab yes\cr
Depends: \tab baySeq, ShortRead\cr
}
To use the package, we construct an \code{\link{alignmentData}} object
(either explicitly or using the \code{\link{processTags}} function).
containing the alignment information for each sample. We then use the
\code{\link{processAD}} function to identify all potential subsegments
of the data and the number of tags that align to these subsegments. We
then empirically determine the prior parameters of the data using the
\code{\link{getPriors}} function, and finally identify all segments
to which a high density of tags align in at least one replicate group
using the \code{\link{segmentSeq}} function. The output from this
segmentation is designed to be usable by the
\code{\link[baySeq:baySeq-package]{baySeq}} package.

The package (optionally) makes use of the 'snow' package for
parallelisation of computationally intensive functions. This is highly
recommended for large data sets.

See the vignette for more details.
}
\author{
Thomas J. Hardcastle

Maintainer: Thomas J. Hardcastle <tjh48@cam.ac.uk>
}
\references{
Hardcastle T.J., and Kelly, K.A. (2010). Genome Segmentation from
High-Throughput Sequencing Data. In submission.
}
\keyword{ package }
\seealso{
  \code{\link[baySeq:baySeq-package]{baySeq}}
}
\examples{

# Define the chromosome lengths for the genome of interest.

chrlens <- c(2e6, 1e6)

# Define the files containing sample information.

datadir <- system.file("data", package = "segmentSeq")
libfiles <- c("SL9.txt", "SL10.txt", "SL26.txt", "SL32.txt")

# Establish the library names and replicate structure.

libnames <- c("SL9", "SL10", "SL26", "SL32")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(libfiles, dir = datadir, replicates, libnames, chrlens, chrs = c(">Chr1", ">Chr2"), header = TRUE)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, maxgaplen = 500, cl = NULL)

}


