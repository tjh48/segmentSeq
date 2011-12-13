\name{heuristicSeg}
\alias{heuristicSeg}
\title{
A (fast) heuristic method for creation of a genome segment map.
}
\description{
This method identifies by heuristic methods a set of loci from a
\code{segData} object. It does this by identifying within replicate groups
regions of the genome that satisfy the criteria for being a locus and
have no region within them that satisfies the criteria for being a
null. These criteria can be defined by the user or inferred from the
data. 
}
\usage{
heuristicSeg(sD, aD, bimodality = TRUE, RKPM = 300, gap = 100, subRegion
= NULL, getLikes = TRUE, verbose = TRUE, cl)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{sD}{
    A \code{\linkS4class{segData}} object derived from the `aD' object.
  }
  \item{bimodality}{Should the criteria for loci be inferred from the
    (likely) bimodal structure of the data?}
  \item{RKPM}{What RKPM (reads per kilobase per million reads) distinguishes
    between a locus and a null region? Ignored if \code{bimodality = TRUE}.}
  \item{gap}{What is the minimum length of a null region? Ignored if
    \code{bimodality = TRUE}.}
  \item{subRegion}{A \code{'data.frame'} object defining the subregions of
    the genome to be segmented. If NULL (default), the whole genome is segmented.}
  \item{getLikes}{Should posterior likelihoods for the new segmented
    genome (loci and nulls) be assessed?}
  \item{verbose}{Should the function be verbose? Defaults to TRUE.}
  \item{cl}{A SNOW cluster object, or NULL. See Details.}
}
\details{
  A \code{'cluster'} object (package: snow) may be used for
  parallelisation of parts of this function when examining large data sets.
  Passing NULL to this variable will cause the function to run in non-parallel mode.
}
\value{
  A \code{\link{lociData}} object, containing count
  information on all the segments discovered.
}
\references{
Hardcastle T.J., Kelly, K.A. and Balcombe D.C. (2011). Identifying small
RNA loci from high-throughput sequencing data. In press.

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
  differential expression in \code{\link{lociData}} objects.
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

# Process the files to produce an `alignmentData' object.

alignData <- readGeneric(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 100)

# Process the alignmentData object to produce a `segData' object.

sD <- processAD(alignData, cl = NULL)

# Use the segData object to produce a segmentation of the genome.

segD <- heuristicSeg(sD = sD, aD = alignData, subRegion = data.frame(chr = ">Chr1", start =
1, end = 1e5), cl = NULL)
}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{classif}
