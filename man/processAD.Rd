\name{processAD}
\alias{processAD}
\title{Processes an 'alignmentData' object into a 'segData' object for segmentation.
}
\description{
In order to discover segments of the genome with a high density of
sequenced data, a 'segData' object must be produced. This is an object
containing a set of potential segments, together with the counts for
each sample in each potential segment.
}
\usage{
processAD(aD, gap = NULL, verbose = TRUE, cl)
}
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{gap}{
    The maximum gap between aligned tags that should be allowed in
    constructing potential segments. See Details.
  }
  
  \item{verbose}{
    Should processing information be displayed? Defaults to TRUE.
  }
  
  \item{cl}{A SNOW cluster object, or NULL. See Details.}
}
  
  \details{
    This function takes an \code{\linkS4class{alignmentData}} object and
    constructs a \code{\linkS4class{segData}} object from it. The function
    creates a set of potential segments by looking for all locations on
    the genome where the start of a region of overlapping alignments
    exists in the \code{\linkS4class{alignmentData}} object. A potential
    segment then exists from this start point to the end of all regions
    of overlapping alignments such that there is no region in the
    segment of at least length 'gap'  where no tag aligns. The
    \code{'gap'} argument thus defines the maximum gap that can exist
    between tags in a segment of high density of alignments. The number
    of potential segments can therefore be increased by increasing this 
    limit, or (usually more usefully) decreased by decreasing this limit
    in order to save computational effort.

    The 'gap' argument is now by default specified in the \code{\link{processTags}}
    function used to create the 'aD' object, and so 'gap' can be left as
    NULL providing this has been done.

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

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 200)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, gap = 200, cl = NULL)

}
\keyword{manip}
