\name{mergeSD}
\alias{mergeSD}
\title{Merges multiple 'segData' objects.}
\description{This file takes two or more segData objects, which may have
  different co-ordinates defined in the '@segInfo' slot and merges them
  into one for subsequent analysis. This is not quite equivalent to
  'cbind' as it may change the co-ordinate structure of the data. 
}
\usage{
mergeSD(..., aD, replicates, gap = 200, cl = NULL)
}
\arguments{
  \item{...}{Two or more 'segData' objects.}
  \item{aD}{An 'alignmentData' object that describes (in the same order)
    all samples that appear in the 'segData' objects described in '...'}
  \item{replicates}{The replicate structure of the new 'segData' object'.}
    \item{gap}{
    The maximum gap between aligned tags that should be allowed in
    constructing potential segments. See Details.
  }
  \item{cl}{A SNOW cluster object, or NULL.}
}

\value{A 'segData' object.}
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

# Process the first two files to produce an 'alignmentData' object.

alignData12 <- processTags(file = libfiles[1:2], dir = datadir, replicates =
replicates[1:2], libnames = libnames[1:2], chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 200)

# Process the alignmentData object to produce a 'segData' object.
sD12 <- processAD(alignData12, cl = NULL)

# Repeat with the third and fourth files

alignData34 <- processTags(file = libfiles[3:4], dir = datadir, replicates =
replicates[3:4], libnames = libnames[3:4], chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens, gap = 200)
sD34 <- processAD(alignData34, cl = NULL)

# Bind the two alignData objects together
alignData <- cbind(alignData12, alignData34)

# Merge the two segData objects
sD <- mergeSD(sD12, sD34, aD = alignData, replicates = c(1,1,2,2), gap =200, cl = NULL)

# Note missing data in new 'segData' object.
sD

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
