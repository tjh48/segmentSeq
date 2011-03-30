\name{getOverlaps}
\alias{getOverlaps}
\title{Identifies overlaps between two sets of genomic coordinates}
\description{
This function identifies which of a set of genomic segments overlaps
with another set of coordinates; either with partial overlap or with the
segments completely contained within the coordinates. The function is
used within the 'segmentSeq' package for various methods of constructing
a segmentation map, but may also be useful in downstream analysis
(e.g. annotation analyses).
}
\usage{
getOverlaps(coordinates, segments, overlapType = "overlapping", whichOverlaps = TRUE, cl)
}
\arguments{
  \item{coordinates}{
    A \code{data.frame} object, with columns 'chr', 'start' and 'end',
    defining the set of coordinates with which the segments may overlap.
}
  \item{segments}{
    A \code{data.frame} object, with columns 'chr', 'start' and 'end',
    defining the set of segments which may overlap within the coordinates.
  }
  \item{overlapType}{
    Which kind of overlaps are being sought? Can be one of
    'overlapping', 'contains' or 'within'. See Details.
  }
  \item{whichOverlaps}{If TRUE, returns the 'segments' overlapping with
    the 'coordinates'. If FALSE, returns a boolean vector specifying
    which of the 'coordinates' overlap with the 'segments'.
    }
  \item{cl}{
    A SNOW cluster object, or NULL. See Details.
}
}
\details{
  If \code{overlapType = "overlapping"} then any overlap between the
  'coordinates' and the 'segments' is sufficient. If \code{overlapType =
    "contains"} then a region defined in 'coordinates' must completely
  contain at least one of the 'segments' to count as an overlap. If
  \code{overlapType = "within"} then a region defined in 'coordinates'
  must be completely contained by at least one of the 'segments' to
  count as an overlap.
  
  A \code{'cluster'} object (package: snow) may be used for
  parallelisation of this function when examining large data sets.
  Passing NULL to this variable will cause the function to run in non-parallel mode.
}
\value{
If \code{whichOverlaps = TRUE}, then the function returns a list object
  with length equal to the number of rows of the 'coordinates'
  argument. The 'i'th member of the list will be a numeric vector giving
  the row numbers of the 'segments' data.frame which overlap with the
  'i'th row of the 'coordinates' data.frame, or NA if no segments
  overlap with this coordinate region.

  If \code{whichOverlaps = FALSE}, then the function returns a boolean
  vector with length equal to the number of rows of the 'coordinates'
  argument, indicating which of the regions defined in coordinates have
  the correct type of overlap with the 'segments'.
}
\author{
Thomas J. Hardcastle
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

# Find which tags overlap with an arbitrary set of coordinates.

getOverlaps(coordinates = data.frame(chr = ">Chr1", start =
c(1,100,2000), end = c(40, 3000, 5000)), segments =
alignData@alignments, overlapType = "overlapping", whichOverlaps = TRUE,
cl = NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
