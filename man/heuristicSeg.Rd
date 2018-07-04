% modification on git from copied files
\name{heuristicSeg}
\alias{heuristicSeg}
\title{
A (fast) heuristic method for creation of a genome segment map.
}
\description{
This method identifies by heuristic methods a set of loci from a
\code{segData} or \code{segMeth} object. It does this by identifying
within replicate groups regions of the genome that satisfy the criteria
for being a locus and have no region within them that satisfies the
criteria for being a null. These criteria can be defined by the user or
inferred from the data. 
}
\usage{
heuristicSeg(sD, aD, gap = 50, RKPM = 1000, prop, coverage = 1, locCutoff =
0.9, nullCutoff = 0.9, subRegion = NULL, largeness = 1e8, getLikes =
TRUE, verbose = TRUE, tempDir = NULL, cl = NULL, recoverFromTemp =
FALSE, trimMeth = FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} or
    \code{\linkS4class{methData}} object.
  }
  \item{sD}{
    A \code{\linkS4class{segData}} or \code{\linkS4class{segMeth}}
    object derived from the `aD' object.
  }
%  \item{bimodality}{Should the criteria for loci be inferred from the
%    (likely) bimodal structure of the data?
\item{gap}{What is the minimum length of a null region?}
% Ignored if \code{bimodality = TRUE}.}
  \item{RKPM}{For analysis of a \code{\linkS4class{segData}} object,
    what RKPM (reads per kilobase per million reads) distinguishes 
    between a locus and a null region?}
% Ignored if \code{bimodality = TRUE}.}
  \item{prop}{For analysis of a \code{\linkS4class{segMeth}} object,
    what proportion of methylated cytosines distinguishes between a
    locus and a null region? (see Details).}
  \item{coverage}{For analysis of a \code{\linkS4class{segMeth}} object,
    what is the minimum coverage required to make inferences on the
    presence/absense of a methylation locus?}
  \item{locCutoff}{For analysis of a \code{\linkS4class{segMeth}}
    object, with what likelihood must the proportion of methylated
    cytosines exceed the `prop' option to define a locus? Defaults to 0.9.}
    \item{nullCutoff}{For analysis of a \code{\linkS4class{segMeth}}
    object, with what likelihood must the proportion of methylated
    cytosines be less than the `prop' option to define a null region? Defaults to 0.9.}
  \item{subRegion}{A \code{'data.frame'} object defining the subregions of
    the genome to be segmented. If NULL (default), the whole genome is
    segmented.}
  \item{largeness}{The maximum size for a split analysis.}
  \item{getLikes}{Should posterior likelihoods for the new segmented
    genome (loci and nulls) be assessed?}
  \item{verbose}{Should the function be verbose? Defaults to TRUE.}
  \item{tempDir}{A directory for storing temporary files produced during
    the segmentation.}
  \item{cl}{A SNOW cluster object, or NULL. Defaults to NULL. See
    Details.}
  \item{recoverFromTemp}{If TRUE, will attempt to recover the position
    saved in 'tempDir'. Defaults to FALSE. See Details.}
  \item{trimMeth}{Should putative methylation regions be trimmed?
    Defaults to FALSE; see Details.}
}

\details{
  A \code{'cluster'} object (package: snow) may be used for
  parallelisation of parts of this function when examining large data sets.
  Passing NULL to this variable will cause the function to run in
  non-parallel mode.

  If \code{recoverFromTemp = TRUE}, the function will attempt to recover
  a crashed position from the temporary files in \code{tempDir}. At
  present, the function assumes you know what you are doing, and will
  perform no checking that these files are suitable for the specified
  recovery. Use with caution.

  The \code{prop} variable can be used to set the proportion of
  methylation required to identify a locus by giving a numerical value
  between 0-1. It can also be determined automatically (see
  thresholdFinder).

  Due to the way that methylation loci are identified, it is possible
  that the cytosines at the edges of methylation loci have limited
  evidence for methylation. The 'trimMeth' option trims cytosines at the
  edge of predicted methylation loci that have less than 50\% likelihood
  of being above the required threshold.
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
  \code{\link{thresholdFinder}}, a function for determining a suitable
  threshold on methylation by examining the data.
  \code{\link{plotGenome}}, a function for plotting the alignment of
  tags to the genome (together with the segments defined by this
  function).
  \code{\link[baySeq:baySeq-package]{baySeq}}, a package for discovering
  differential expression in \code{\link{lociData}} objects.
}

\examples{
# Define the files containing sample information.

datadir <- system.file("extdata", package = "segmentSeq")
libfiles <- c("SL9.txt", "SL10.txt", "SL26.txt", "SL32.txt")

# Establish the library names and replicate structure.

libnames <- c("SL9", "SL10", "SL26", "SL32")
replicates <- c(1,1,2,2)

# Process the files to produce an `alignmentData' object.

alignData <- readGeneric(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, gap = 100)

# Process the alignmentData object to produce a `segData' object.

sD <- processAD(alignData, gap = 100, cl = NULL)

# Use the segData object to produce a segmentation of the genome.

segD <- heuristicSeg(sD = sD, aD = alignData, prop = 0.2,
subRegion = data.frame(chr = ">Chr1", start = 1, end = 1e5),
cl = NULL)
}
% Add one or more standard keywords, see file `KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{classif}
