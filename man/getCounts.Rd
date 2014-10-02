\name{getCounts}
\alias{getCounts}
\title{Gets counts from alignment data from a set of genome segments.}

\description{A function for extracting count data from an
  \code{alignmentData} object given a set of segments defined on the genome.}
\usage{
getCounts(segments, aD, preFiltered = FALSE, adjustMultireads = TRUE,
  useChunk = FALSE, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segments}{
    A \code{GRanges} object which defines a set of segments for
    which counts are required.
  }
  \item{aD}{
    An \code{\link{alignmentData}} object.
  }
  \item{preFiltered}{The function internally cleans the data; however,
    this may not be needed and omitting these steps may save computational
    time. See Details.}
  \item{adjustMultireads}{If working with methylation data, this option
    toggles an adjustment for reads that align to multiple locations on
    the genome. Defaults to TRUE.}
  \item{useChunk}{If all segments are within defined `chunks' of the
  alignmentData object, speed increases if this is set to
  TRUE. Otherwise, counts may be inaccurate. Defaults to FALSE.}
  \item{cl}{A SNOW cluster object, or NULL. See Details.}    
}
\details{
  The function extracts count data from \code{alignmentData} object
  'aD' given a set of segments. The non-trivial aspect of this function
  is that at a segment which contains a tag that matches to multiple
  places in that segment (and thus appears multiple times in the
  \code{alignmentData} object) should count it only once.

  If \code{preFiltered = FALSE} then the function allows for missing (NA)
  data in the segments, unordered segments and duplicated segments. If
  the segment list has no missing data, is already ordered, and contains
  no duplications, then computational time can be saved by setting
  \code{preFiltered = TRUE}.
  
  A \code{cluster} object (package: snow) is recommended for
  parallelisation of this function when using large data sets.
  Passing NULL to this variable will cause the function to run in non-parallel mode.
  
  In general, this function will probably not be accessed by the user as
  the \code{\link{processAD}} function includes a call to \code{getCounts} as
  part of the standard processing of an \code{alignmentData} object into
  a \code{segData} object.
  }

\value{
If `as.matrix', a matrix, each column of which corresponds to a library in the
\code{alignmentData} object `aD' and each row to the
segment defined by the corresponding row in `segments'. Otherwise an
  equivalent \code{DataFrame} object. 
}

\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{processAD}}
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
chrlens, gap = 100)

# Get count data for three arbitrarily chosen segments on chromosome 1.

getCounts(segments = GRanges(seqnames = c(">Chr1"),
          IRanges(start = c(1,100,2000), end = c(40,3000,5000))), 
          aD = alignData, cl = NULL) 

}

\keyword{manip}
