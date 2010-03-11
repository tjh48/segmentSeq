\name{getCounts}
\alias{getCounts}
\title{Gets counts from alignment data from a set of genome segments.}

\description{A function for extracting count data from an
  \code{'alignmentData'} object given a set of segments defined on the genome.}
\usage{
getCounts(segments, aD, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segments}{
    A \code{'data.frame'} object which defines a set of segments for
    which counts are required.
  }
  \item{aD}{
    An \code{\link{alignmentData}} object.
  }
  \item{cl}{A SNOW cluster object, or NULL. See Details.}
    
}
\details{
  The function extracts count data from \code{alignmentData} object
  'aD' given a set of segments. The non-trivial aspect of this function
  is that at a segment which contains a tag that matches to multiple
  places in that segment (and thus appears multiple times in the
  \code{alignmentData} object should count it only once.

  A \code{'cluster'} object (package: snow) is recommended for
  parallelisation of this function when using large data sets.
  Passing NULL to this variable will cause the function to run in non-parallel mode.

  In general, this function will probably not be accessed by the user as
  the \code{\link{processAD}} function includes a call to 'getCounts' as
  part of the standard processing of an \code{alignmentData} object into
  a \code{segData} object.
}
\value{
A matrix, each column of which corresponds to a library in the
  \code{alignmentData} object 'aD' and each row to the
  segment defined by the corresponding row in the data.frame 'segments'.
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

datadir <- system.file("data", package = "segmentSeq")
libfiles <- dir(datadir, pattern = ".txt", full.names = TRUE)

# Establish the library names and replicate structure.

libnames <- c("SL10", "SL26", "SL32", "SL9")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(libfiles, replicates, libnames, chrlens, chrs = c(">Chr1", ">Chr2"), header = TRUE)

# Get count data for three arbitrarily chosen segments on chromosome 1.

getCounts(segments = data.frame(chr = ">Chr1", start = c(1,100,2000), end =
c(40, 3000, 5000)), aD = alignData, cl = NULL)

}

\keyword{manip}
