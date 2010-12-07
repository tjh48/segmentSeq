\name{clustSeg}
\alias{clustSeg}
\title{
A (fast) heuristic method for creation of a genome segment map.
}
\description{
This method uses kmeans clustering (with user defined initial centres)
on the densities of sequences mapping to the genome in order to identify
segments of the genome which have a higher than background intensity,
and thus construct a genome segment map.
}
\usage{
clustSeg(sD, aD, subRegion = NULL, getLikes = TRUE, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{sD}{
    A \code{\linkS4class{segData}} object derived from the 'aD' object.
}
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{subRegion}{A \code{'data.frame'} object defining the subregions of
    the genome to be segmented. If NULL (default), the whole genome is segmented.}
  \item{getLikes}{Should posterior likelihoods for the new segmented
    genome (loci and nulls) be assessed?}
  \item{cl}{A SNOW cluster object, or NULL. See Details.}
}
\details{
  A \code{'cluster'} object (package: snow) may be used for
  parallelisation of parts of this function when examining large data sets.
  Passing NULL to this variable will cause the function to run in non-parallel mode.
}
\value{
  A \code{\link[baySeq:countData]{countData}} object, containing count
  information on all the segments discovered.
}
\references{
    Hardcastle T.J., and Kelly, K.A. (2010). Genome Segmentation From
    High-Throughput Sequencing Data. In preparation.
}
\author{
Thomas J. Hardcastle
}
\seealso{
  \code{\link{classifySeg}}, an alternative approach to this problem using an empirical
  Bayes approach to classify segments.
  \code{\link{plotGenome}}, a function for plotting the alignment of
  tags to the genome (together with the segments defined by this
  function).
  \code{\link[baySeq:baySeq-package]{baySeq}}, a package for discovering
  differential expression in \code{\link[baySeq:countData]{countData}} objects.
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

# Use the segData object to produce a segmentation of the genome.

segD <- clustSeg(sD = sD, aD = alignData, subRegion = data.frame(chr = ">Chr1", start =
1, end = 1e5), cl = NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{classif}
