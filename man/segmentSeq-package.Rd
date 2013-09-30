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
aligning to them. The package was developed for use in identifying small
RNA precursors from small RNA sequencing data, but may also be useful in
some mRNA-Seq and chIP-Seq applications.
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
from sets of alignment files using either the \code{\link{readGeneric}}
function to read text files or the \code{\link{readBAM}} function to
read from BAM format files.

We then use the \code{\link{processAD}} function to identify all
potential subsegments of the data and the number of tags that align to
these subsegments. We then use either a heuristic or empirical Bayesian
approach to segment the genome into `loci' and `null' regions. We can then
acquire posterior likelihoods for each set of replicates which tell us
whether a region is likely to be a locus or a null in that replicate group.

The segmentation is designed to be usable by the
\code{\link[baySeq:baySeq-package]{baySeq}} package to allow
differential expression analyses to be carried out on the discovered loci.

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
Hardcastle T.J., Kelly, K.A. and Balcombe D.C. (2011). Identifying small
RNA loci from high-throughput sequencing data. In press.
}
\keyword{ package }
\seealso{
  \code{\link[baySeq:baySeq-package]{baySeq}}
}
\examples{

# Define the chromosome lengths for the genome of interest.

chrlens <- c(2e6, 1e6)

# Define the files containing sample information.

datadir <- system.file("extdata", package = "segmentSeq")
libfiles <- c("SL9.txt", "SL10.txt", "SL26.txt", "SL32.txt")

# Establish the library names and replicate structure.

libnames <- c("SL9", "SL10", "SL26", "SL32")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- readGeneric(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, gap = 100, cl = NULL)

}


