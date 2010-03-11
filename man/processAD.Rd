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
processAD(aD, maxgaplen = 500, maxloclen = NULL, verbose = TRUE, cl = cl)
}
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{maxgaplen}{
    The maximum gap between aligned tags that should be allowed in
    constructing potential segments. See Details.
  }
  
  \item{maxloclen}{
    The maximum length that any potential segment may be. If NULL
    (recommended) no such limit exits. See Details.
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
    segment of at least length 'maxgaplen' where no tag aligns. The
    \code{'maxgaplen'} argument thus defines the maximum gap that can exist
    between tags in a segment of high density of alignments. The number
    of potential segments can therefore be increased by increasing this
    limit, or (usually more usefully) decreased by decreasing this limit
    in order to save computational effort.

    The number of potential segments may also be decreased by setting a
    maximum length for any potential segments by setting the argument
    'maxloclen'. The use of this argument is not recommended as it
    appears to substatially degrade the results.

    The number of potential segments created can be further reduced by
    setting a limit on the maximum length that any segment may be with
    the \code{'maxloclen'} argument. The use of this limit tends to have
    severe negative effects on the final segmentation results, however,
    and is therefore not recommended.
    
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
  \code{\link{segmentSequences}}, which segments the genome based on the
  \code{segData} object produced.
  \code{\linkS4class{segData}}
  \code{\linkS4class{alignmentData}}
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

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, maxgaplen = 500, cl = NULL)

}
\keyword{manip}
