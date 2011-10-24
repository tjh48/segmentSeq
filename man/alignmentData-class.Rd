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
\alias{cbind}
\alias{cbind,alignmentData-method}

\title{Class "alignmentData"}
\description{The \code{alignmentData} class records information about
  a set of alignments of high-throughput sequencing data to a
  genome. Details include the alignments themselves, details on the
  chromosomes of the genome to which the data are aligned, and
  information on the libraries from which the data come.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("alignmentData",
  ...)}, but more usually by using one of \code{\link{readBAM}} or
\code{\link{readGeneric}} functions to generate the object from a set of
alignment files.
}
\section{Slots}{
  \describe{
    \item{\code{alignments}:}{Object of class \code{"GRanges"}.
      Stores information about the alignments. See Details.}
    \item{\code{data}:}{Object of class \code{"DataFrame"}. For each
      alignment described in the \code{alignments} slot, contains the
      number of times the alignment is seen in each sample.}
    \item{\code{libnames}:}{Object of class \code{"character"}. The
      names of the libraries for which alignment data exists.}
    \item{\code{libsizes}:}{Object of class \code{"numeric"}. The
      library sizes (see Details) for each of the libraries.}
    \item{\code{replicates}:}{Object of class \code{"factor"}.
      Replicate information for each of the libraries. See Details.}
  }
}

\section{Details}{
  The \code{alignments} slot is the key element of this class. This is a
  \code{GRanges} object that, in addition to the usual elements defining
  the location of aligned objects to a reference genome, also describes
  the values `tag', giving the sequence of the tag aligning to the
  location, `matches', indicating in how many places that tag matches to
  the genome, `chunk', an identifier for the sets of tags that align
  close enough together to form a potential locus, and `chunkDup',
  indicating whether that tag matches to multiple places within the chunk.

  The library sizes, defined in the \code{libsizes} slot, provide some
  scaling factor for the observed number of counts of a tag in different
  samples.

  The \code{replicates} slot is a vector of factors such that the ith
  sample is a replicate of the jth sample if and only if \code{@replicates[i] ==
    @replicates[j]}.
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

\seealso{
  \code{\link{readGeneric}}, which will produce a \code{'alignmentData'}
  object from appropriately formatted tab-delimited files.
  \code{\link{readBAM}}, which will produce a \code{'alignmentData'}
  object from BAM files.
  \code{\link{processAD}}, which will convert an \code{'alignmentData'}
  object into a 'segData' object for segmentation.
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
chrlens)

}
\keyword{classes}

