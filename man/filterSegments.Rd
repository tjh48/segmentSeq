\name{filterSegments}
\alias{filterSegments}

\title{
Filters an set of segments (given an ordering on the segments) such that
no segments overlap.
}
\description{
This function takes a set of segments, plus an ordering on that set, and
filters the set such that no segments overlap, preferentially keeping
the segments first in ordering.
}
\usage{
filterSegments(segs, orderOn, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{
    A \code{'data.frame'} containing columns 'chr', 'start' and 'end'.
  }
  \item{orderOn}{
    An vector of some statistic that can be used to create an ordering
    on the 'segs' data.frame.
  }
  \item{...}{Additional parameters that can be passed to \code{order}
    when ordering the segments using the 'orderOn' parameter.}

}
\details{
  This function takes the set of segments defined by the
  \code{data.frame} 'segs', together with some statistic (e.g.,
  likelihood of similarity with background) defined by the
  'orderOn' vector. Additional options can be passed to the 'order'
  function (for example, relating to the direction of the ordering)
  through the '...' parameter.

  The function takes the segment first in the ordering and discards any
  segments that overlap with it. It then proceeds to the next remaining
  segment in the ordering and discards any segments that overlap with
  this. This process continues until we have a set of non-overlapping
  segments.

  This function can be used to create a random sample of non-overlapping
  segments by providing a randomly chosen set of values for the
  'orderOn' vector.

}
\value{
  A vector giving the rows of the \code{data.frame} object 'segs' which
  form a non-overlapping set.
}
\author{
Thomas J. Hardcastle
}


\examples{

# Define the chromosome lengths for the genome of interest.

chrlens <- c(2e6, 1e6)

# Define the files containing sample information.

datadir <- system.file("data", package = "segmentSeq")
libfiles <- dir(datadir, pattern = ".txt", full.names = TRUE)

# Establish the library names and replicate structure.

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

# Create random sampling of non-overlapping segments for chromosome 1 of sD object.

filterSegments(subset(sD@segInfo, select = c(chr, start, end)), runif(nrow(sD)))

}
\keyword{manip}
\keyword{misc}
