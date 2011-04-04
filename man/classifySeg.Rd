\name{classifySeg}
\alias{classifySeg}
\title{
A method for defining a genome segment map by an empirical Bayesian
classification method
}
\description{
This function acquires empirical distributions of sequence tag density
from an already existing (or heuristically defined) segment map. It uses
these to classify potential segments as either segments or nulls in
order to define a new (and improved) segment map.
}
\usage{
classifySeg(sD, cD, aD, lociCutoff = 0.9, nullCutoff = 0.9, subRegion =
NULL, getLikes = TRUE, lR = FALSE, samplesize = 1e5, cl, ...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sD}{
    A \code{\linkS4class{segData}} object derived from the 'aD' object.
}
  \item{cD}{
    A \code{\link[baySeq:countData]{countData}} object containing an
    already existing segmentation map, or NULL.
}
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{lociCutoff}{The minimum posterior likelihood of being a locus
    for a region to be treated as a locus.}
  \item{nullCutoff}{The minimum posterior likelihood of being a null for
    a region to be treated as a null.}
  \item{subRegion}{A \code{'data.frame'} object defining the subregions of
    the genome to be segmented. If NULL (default), the whole genome is
    segmented.}
  \item{getLikes}{Should posterior likelihoods for the new segmented
    genome (loci and nulls) be assessed?}
  \item{lR}{If TRUE, locus and null calls are made on the basis of
    likelihood ratios rather than posterior likelihoods. Not
    recommended.}
  \item{samplesize}{The sample size to be used when estimating the prior
    distribution of the data with the
    \code{\link[baySeq:getPriors.NB]{getPriors.NB}} function.}
  \item{cl}{A SNOW cluster object, or NULL. See Details.}
  \item{...}{Any additional parameters to be passed to \code{\link{heuristicSeg}}.}
}
\details{
This function acquires empirical distributions of sequence tag density
from the segmentation map defined by the 'cD' argument (if 'cD = NULL' then the
\code{\link{heuristicSeg}} function is used to define a segmentation map. It
uses these empirical distributions to acquire posterior likelihoods on
each potential segment being either a true segment or a null
region. These posterior likelihoods are then used to define the segment map.
}
\value{
  A \code{\link{postSeg}} object, containing the
  segmentation map discovered. 
}
\references{
  Hardcastle T.J., and Kelly, K.A. (2010). Genome Segmentation From
  High-Throughput Sequencing Data. In preparation.
}
\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{heuristicSeg}} fast heuristic alternative to this function.
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

datadir <- system.file("extdata", package = "segmentSeq")
libfiles <- c("SL9.txt", "SL10.txt", "SL26.txt", "SL32.txt")

# Establish the library names and replicate structure.

libnames <- c("SL9", "SL10", "SL26", "SL32")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 200)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, cl = NULL)

# Use the classifySeg function on the segData object to produce a postSeg object.

pS <- classifySeg(aD = alignData, sD = sD, subRegion = data.frame(chr = ">Chr1", start = 1, end = 1e5), getLikes = TRUE, cl = NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{classif}% __ONLY ONE__ keyword per line
