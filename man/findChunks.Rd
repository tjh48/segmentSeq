\name{findChunks}
\alias{findChunks}
\title{Identifies `chunks' of data within a set of aligned reads.}
\description{
This function identifies chunks of data within a set of aligned reads by looking for gaps within the alignments; regions where no reads align. If we assume that a locus should not contain a gap of sufficient length, then we can separate the analysis of the data into chunks defined by these gaps, reducing the complexity of the problem of segmentation.
}
\usage{
findChunks(alignments, gap, checkDuplication = TRUE)
}
\arguments{
  \item{alignments}{A \code{\linkS4class{GRanges}} object
    defining a set of aligned reads.}
  \item{gap}{The minimum length of a gap across which it is assumed that
    no locus can exist.}
  \item{checkDuplication}{Should we check whether or not reads are
    duplicated within a chunk? Defaults to TRUE.}
}

\details{This function is called by the \code{\link{readGeneric}} and
  \code{\link{readBAM}} functions but may usefully be called again if
  filtering of an \code{linkS4class{alignmentData}} object has altered
  the data present, or to increase the computational effort required for
  subsequent analysis. The lower the `gap' parameter used to define the
  chunks, the faster (though potentially less accurate) any subsequent
  analyses will be.
}
\value{
A modified \code{\link{GRanges}} object, now containing columns `chunk'
  and `chunkDup' (if 'checkDuplication' is TRUE), identifying the chunk
  to which the alignment belongs and whether the alignment of the tag is
  duplicated within the chunk respectively.
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

# Read the files to produce an `alignmentData' object.

alignData <- readGeneric(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 100)

# Filter the data on number of matches of each tag to the genome

alignData <- alignData[values(alignData@alignments)$matches < 5,]

# Redefine the chunking structure of the data.

alignData <- findChunks(alignData@alignments, gap = 100)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
