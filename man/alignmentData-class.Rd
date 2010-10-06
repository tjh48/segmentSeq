\name{alignmentData-class}
\Rdversion{1.1}
\docType{class}
\alias{alignmentData-class}
\alias{alignmentData}
\alias{[,alignmentData-method}
\alias{[,alignmentData,ANY,ANY-method}
\alias{dim,alignmentData-method}
\alias{initialize,alignmentData-method}
\alias{show,alignmentData-method}

\title{Class "alignmentData"}
\description{The \code{alignmentData} class records information about
  a set of alignments of high-throughput sequencing data to a
  genome. Details include the alignments themselves, details on the
  chromosomes of the genome to which the data are aligned, and
  information on the libraries from which the data come.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("alignmentData",
  ...)}, but more usually by using the \code{\link{processTags}} function.
}
\section{Slots}{
  \describe{
    \item{\code{alignments}:}{Object of class \code{"data.frame"}.
      Stores information about the alignments. See Details.}
    \item{\code{data}:}{Object of class \code{"matrix"}. For each
      alignment described in the \code{alignments} slot, contains the
      number of times the alignment is seen in each sample.}
    \item{\code{libnames}:}{Object of class \code{"character"}. The
      names of the libraries for which alignment data exists.}
    \item{\code{libsizes}:}{Object of class \code{"numeric"}. The
      library sizes (see Details) for each of the libraries.}
    \item{\code{chrs}:}{Object of class \code{"character"}. The
      chromosome names. Must be given as \code{"character"} to
      accomodate non-standard chromosome names (such as "X"}
    \item{\code{chrlens}:}{Object of class \code{"numeric"}. The
      lengths of each of the chromosomes defined by slot \code{chrs}.}
    \item{\code{replicates}:}{Object of class \code{"numeric"}.
      Replicate information for each of the libraries. See Details.}
  }
}

\section{Details}{
  The \code{alignments} slot is the key element of this class. This is a
  \code{"data.frame"} object that contains the columns 'chr', 'start',
  'end', 'duplicated', 'tag', 'count', 'sampleNumber' and
  'replicate'. Columns 'chr', 'start' and 'end' define the chromosome,
  start and end point of the tag. 'duplicated' indicates whether or not
  the tag uniquely matches this location (\code{FALSE}) or whether the
  tag matches some other location on the genome (\code{TRUE}). The 'tag'
  column gives the sequence of the tag as a factor. The 'count' column
  gives the number of times the tag appears in the library. Which
  library is involved is specified by the 'sampleNumber' column,
  and the 'replicate' column gives the replicate group that this library
  is associated with.

  The library sizes, defined in the \code{libsizes} slot, provide some
  scaling factor for the observed number of counts of a tag in different
  samples. One method of calculating this, for example, would be to take
  the number of sequences read from the high-throughput sequencing
  machine that align to the reference genome.

  The \code{replicates} slot should take the form of a vector of
  integers such that if and only if the ith sample is a replicate of the
  jth sample then \code{@replicates[i] == @replicates[j]}. In addition,
  values in the replicates slot should take values from \code{1:n} where
  \code{n} is the number of replicate groups.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "alignmentData")}: ... }
    \item{dim}{\code{signature(x = "alignmentData")}: ... }
    \item{initialize}{\code{signature(.Object = "alignmentData")}: ... }
    \item{show}{\code{signature(object = "alignmentData")}: ... }
	 }
}

\author{Thomas J. Hardcastle}
\note{Methods 'new', 'dim', '[' and 'show' have been defined for these classes.}

\seealso{
  \code{\link{processTags}}, which will produce a \code{'alignmentData'}
  object from appropriately formatted tab-delimited files.
  \code{\link{processAD}}, which will convert an \code{'alignmentData'}
  object into a 'segData' object for segmentation.
}
\examples{

# Define the chromosome lengths for the genome of interest.

chrlens <- c(2e6, 1e6)

# Define the files containing sample information.

datadir <- system.file("data", package = "segmentSeq")
libfiles <- c("SL9.txt", "SL10.txt", "SL26.txt", "SL32.txt")

# Establish the library names and replicate structure.

libnames <- c("SL9", "SL10", "SL26", "SL32")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(libfiles, dir = datadir, replicates, libnames, chrlens, chrs = c(">Chr1", ">Chr2"), header = TRUE)

}
\keyword{classes}

