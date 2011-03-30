\name{lociLikelihoods}
\alias{lociLikelihoods}
\title{Evaluates the posterior likelihoods of each region defined by a
  segmentation map as a locus.}
\description{
  An empirical Bayesian approach that takes a segmentation map and uses
  this to bootstrap posterior likelihoods on each region being a locus
  for each replicate group.
}
\usage{
lociLikelihoods(cD, aD, newCounts = FALSE, bootStraps = 1,
                inferNulls = TRUE, nasZero = FALSE, cl)
}
\arguments{
  \item{cD}{
    A \code{\link[baySeq:countData]{countData}} or
    \code{\link{postSeg}} object that defines a segmentation map.
  }
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{newCounts}{Should new counts be evaluated for the segmentation
    map in 'cD' before calculating loci likelihoods? Defaults to FALSE}
  \item{bootStraps}{What level of bootstrapping should be carried out on
    the inference of posterior likelihoods? See
    \code{\link[baySeq:getLikelihoods.NB]{getLikelihoods.NB}}.}
  \item{inferNulls}{Should null regions be inferred from the gaps
    between segments defined by the 'cD' object?}
  \item{nasZero}{If FALSE, any locus with a posterior likelihood 'NA' in
    the existing segmentation map is treated as a null region for the
    first bootstrap; If TRUE, it is ignored for the first bootstrap.}
  \item{cl}{
    A SNOW cluster object, or NULL. See Details.
}
}
\details{
  A \code{'cluster'} object (package: snow) may be used for
  parallelisation of this function when examining large data sets.
  Passing NULL to this variable will cause the function to run in non-parallel mode.
}
\value{
A \code{\link{postSeg}} object.
}
\author{
Thomas J. Hardcastle
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

alignData <- processTags(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 200)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, cl = NULL)

# Use the segData object to produce a segmentation of the genome, but
# without evaluating posterior likelihoods.

segD <- heuristicSeg(sD = sD, aD = alignData,
    subRegion = data.frame(chr= ">Chr1", start = 1, end = 1e5),
    getLikes = FALSE, cl = NULL) 

# Use the postSeg function to evaluate the posterior likelihoods directly.

postSeg <- lociLikelihoods(segD, aD = alignData, bootStraps = 5,
inferNulls = TRUE, cl = NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
