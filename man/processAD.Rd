\name{processAD}
\alias{processAD}
\title{Processes an `alignmentData' or `alignmentMeth' object into a
  `segData' or `segMeth' object for segmentation.
}
\description{
In order to discover segments of the genome with a high density of
sequenced data, a `segData' object must be produced. This is an object
containing a set of potential segments, together with the counts for
each sample in each potential segment.
}
\usage{
processAD(aD, gap, squeeze = 0, filterProp = 0.1, strandSplit = FALSE, verbose = TRUE, cl)
}
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} or
    \code{\linkS4class{alignmentMeth}} object.
  }
  \item{gap}{
    The maximum gap between aligned tags that should be allowed in
    constructing potential segments. See Details.
  }
  \item{squeeze}{
    If greater than zero, the minimum gap between aligned tags that
    should be allowed in constructing potential segments. See Details.
  }
  \item{filterProp}{
    If 'aD' is a \code{\linkS4class{alignmentMeth}} object and this is
    given, the minimum proportion of methylation at a base below which
    the base will be filtered out before constructing potential
    segments (but not during counting).}
  \item{strandSplit}{
    If TRUE, the data will be split by strand and segments will be
    constructed separately for each strand. Defaults to FALSE.}
  \item{verbose}{
    Should processing information be displayed? Defaults to TRUE.
  }
  
  \item{cl}{A SNOW cluster object, or NULL. See Details.}
}
  
  \details{
    This function takes an \code{\linkS4class{alignmentData}} or
    \code{\linkS4class{alignmentMeth}} object and
    constructs a \code{\linkS4class{segData}} or
    \code{\linkS4class{segMeth}} object from it. The function
    creates a set of potential segments by looking for all locations on
    the genome where the start of a region of overlapping alignments
    (or, if `squeeze' is non-zero, those alignments separated by no more
    than `squeeze'.) exists in the \code{\linkS4class{alignmentData}} object. A potential
    segment then exists from this start point to the end of all regions
    of overlapping alignments such that there is no region in the
    segment of at least length `gap'  where no tag aligns. The number
    of potential segments can therefore be increased by increasing this 
    limit, or (usually more usefully) decreased by decreasing this limit
    in order to save computational effort.

    A \code{'cluster'} object (package: snow) is recommended for
    parallelisation of this function when using large data sets.
    Passing NULL to this variable will cause the function to run in non-parallel mode.
}
\value{
  A \code{\linkS4class{segData}} object.
}

\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{getCounts}}, which produces the count data for each
  potential segment.
  \code{\link{heuristicSeg}} and \code{\link{classifySeg}}, which segment the
  genome based on the \code{segData} object produced by this function
  \code{\linkS4class{segData}}
  \code{\linkS4class{alignmentData}}
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

sD <- processAD(alignData, gap = 100, cl = NULL)

}
\keyword{manip}
