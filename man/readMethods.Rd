\name{readMethods}
\alias{readGeneric}
\alias{readBAM}


\title{Functions for processing files of various formats into an `alignmentData' object.
}
\description{
These functions take alignment files of various formats to produce an object
(see Details) describing the alignment of sequencing tags from
  different libraries. At present, BAM and text files are supported.

}
\usage{
readGeneric(files, dir = ".", replicates, libnames, chrs, chrlens, cols,
            header = TRUE, minlen = 15, maxlen = 1000,
            multireads = 1000, polyLength, estimationType = "quantile",
            discardTags = FALSE, verbose = TRUE, filterReport = NULL, ...)

readBAM(files, dir = ".", replicates, libnames, chrs, chrlens, countID = NULL,
        minlen = 15, maxlen = 1000, multireads = 1000,
        polyLength, estimationType = "quantile", discardTags = FALSE, verbose = TRUE,
        filterReport = NULL) 

}
\arguments{
  \item{files}{
    Filenames of the files to be read in.
  }
  \item{dir}{
    Directory (or directories) in which the files can be found.
  }
  \item{replicates}{
      A vector defining the replicate structure if the group. If and
      only if the ith library is a replicate of the jth library 
      then \code{@replicates[i] == @replicates[j]}. This argument
      may be given in any form but will be stored as a factor.
}
  \item{libnames}{
    Names of the libraries defined by the file names.
}
  \item{chrs}{
    A chracter vector defining (a selection of) the chromosome names
    used in the alignment files.
}
  \item{chrlens}{
    Lengths of the chromosomes to which the alignments were made.
  }
  \item{cols}{
    A named character vector which describes which column of the input
    files contains which data. See Details.}
  \item{countID}{
    A (two-character) string used by the BAM file to identify the
    `counts' of individual sequenced reads; that is, how many times a
    given read appears in the sequenced library. If NULL, it is assumed
    that the data are redundant (see Details).}
  \item{header}{Do the input files have a header line? Defaults
    to TRUE. See Details.}
  \item{minlen}{Minimum length for reads.}
  \item{maxlen}{Maximum length for reads.}
  \item{multireads}{The functions will discard any read that aligns to
    the genome in more locations than given by this value. Set to Inf to
    keep everything. Defaults to 1000.}
  \item{polyLength}{
    If given, an integer value N defining the length of (approximate)
    homopolymers which will be removed from the data. If a tag contains
    a sequence of N+1 reads consisting of at least N identical bases, it
    will be removed. If not given, all data is used.
  }
  \item{estimationType}{
    The estimationType that will be used by the `baySeq' function
    \code{\link[baySeq]{getLibsizes}} to infer the library sizes of the
    samples.
  }
  \item{discardTags}{
    If TRUE, information about the sequence of the aligned reads will be
    discarded. Useful for very large data sets. Defaults to FALSE.
    }
  \item{verbose}{
    Should processing information be displayed? Defaults to TRUE.
  }
  \item{filterReport}{
    If not NULL, this should be a string defining a file to which will
    be written those data filtered on the basis of chromsome choices,
    widths of sequences, multireads or polyBase.}
  \item{...}{Additional parameters to be passed to
    \code{\link{read.table}}. In particular, the `sep' and `skip'
    arguments may be useful.}
}
\details{

  readBAM:
  This function takes a set of BAM files and generates the
  \code{'alignmentData'} object from these. If a character string for
  `countID' is given, the function assumes the data are non-redundant
  and that `countID' identifies the count data (i.e., how many times
  each read appears in the sequenced library) in each BAM file. If
  `countID' is NULL, then it is assumed that the data are redundant, and
  the count data are inferred from the file.
  
  readGeneric:  
  The purpose of this function is to take a set of plain text files
  and produce an \code{'alignmentData'} object. The function uses
  \code{\link{read.table}} to read in the columns of data in the files
  and so by default columns are separated by any white
  space. Alternative separators can be used by passing the appropriate
  value for \code{'sep'} to \code{\link{read.table}}.

  The files may contain columns with column names
  \code{'chr'}, \code{'tag'}, \code{'count'}, \code{'start'},
  \code{'end'}, \code{'strand'} in which case the `cols' argument can be
  ommitted and `header' set to TRUE. If this is the case, there is no
  requirement for all the files to have the same ordering of columns
  (although all must have these column names).

  Alternatively, the columns of data in the input files can be specified by
  the `cols' argument in the form of a named character vector (e.g;
  \code{'cols = c(chr = 1, tag = 2, count = 3, start = 4, end = 5,
  strand = 6)'} would cause the function to assume that the first column
  contains the chromosome information, the second column contained the
  tag information, etc. If `cols' is specified then information in the
  header is ignored. If  `cols' is missing and `header' is FALSE, then it
  is assumed that the data takes the form described in the example above.

  The \code{'tag'}, \code{'count'} and \code{'strand'} columns may optionally be
  omitted from either the file column headers or the `cols' argument. If
  the \code{'tag'} column is omitted, then the data will not account for
  duplicated sequences when estimating the number of counts in loci. If
  the \code{'count'} column is omitted, the \code{'readGeneric'} function
  will assume that the file contains the alignments of each copy of each
  sequence tag, rather than an aggregated alignment of each unique
  sequence. The unique alignments will be identified and the number of
  sequence tags aligning to each position will be calculated. If
  \code{'strand'} is omitted, the strand will simply be ignored.

}
\value{
An \code{alignmentData} object.
}
\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{alignmentData}}
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

# Process the files to produce an `alignmentData' object.

alignData <- readGeneric(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames, chrs = c(">Chr1", ">Chr2"), chrlens =
chrlens)

}
\keyword{files}
