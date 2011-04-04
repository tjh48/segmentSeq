\name{findChunks}
\alias{findChunks}
\title{Identifies 'chunks' of data within an 'alignmentData' object.}
\description{
This function identifies chunks of data within an 'alignmentData' object by looking for gaps within the alignments; regions where no tags align. If we assume that a locus should not contain a gap of sufficient length, then we can separate the analysis of the data into chunks defined by these gaps, reducing the complexity of the problem of segmentation.
}
\usage{
findChunks(aD, gap)
}
\arguments{
  \item{aD}{An \code{\linkS4class{alignmentData}} object.}
  \item{gap}{The minimum length of a gap across which it is assumed that no locus can exist.}
}

\details{This function is called by the \code{\link{processTags}} function but may usefully be called again if filtering of an \code{linkS4class{alignmentData}} object has altered the data present, or to increase the computational effort required for subsequent analysis. The lower the 'gap' parameter used to define the chunks, the faster any subsequent analyses will be.
}
\value{
A modified \code{\linkS4class{alignmentData}} object, in which the '@alignments' slot contains columns 'chunk' and 'chunkDup', identifying the chunk to which the alignment belongs and whether the alignment of the tag is duplicated within the chunk respectively.
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

# Filter the data on number of matches of each tag to the genome

alignData <- alignData[alignData@alignments$matches < 5,]

# Redefine the chunking structure of the data.

alignData <- findChunks(alignData, gap = 200)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
