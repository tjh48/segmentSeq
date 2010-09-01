\name{processTags}
\alias{processTags}

\title{Convenience function for processing tab-delimited files in a
  certain format into an 'alignmentData' object.
}
\description{
This function takes files in a text format with defined columns
(see Details) that describe the alignment of sequencing tags from
  different libraries.

}
\usage{
processTags(files, replicates, libnames, chrs, chrlens, cols, header =
TRUE, verbose = TRUE, ...)
}
\arguments{
  \item{files}{
    Filenames of the files to be read in.
}
  \item{replicates}{
    Replicate information on the libraries. See Details. 
}
  \item{libnames}{
    Names of the libraries defined by the file names.
}
  \item{chrs}{
    Chromosome names (as \code{'character'}) used in the alignment files.
}
  \item{chrlens}{
    Lengths of the chromosomes to which the alignments were made.
  }
  \item{cols}{
    A named character vector which describes which column of the input
    files contains which data. See Details.}
  
  \item{header}{Do the input files have a header line? Defaults
    to TRUE. See Details.}
  
  \item{verbose}{
    Should processing information be displayed? Defaults to TRUE.
  }

  
  \item{...}{Additional parameters to be passed to \code{\link{read.table}}.}
}
\details{
  The purpose of this function is to take a set of plain text files
  and produce an \code{'alignmentData'} object. The function uses
  \code{\link{read.table}} to read in the columns of data in the files
  and so by default columns are separated by any white
  space. Alternative separators can be used by passing the appropriate
  value for \code{'sep'} to \code{\link{read.table}}.

  The files may contain columns with column names
  \code{'chr'}, \code{'tag'}, \code{'count'}, \code{'start'},
  \code{'end'}, in which case the 'cols' argument can be ommitted and
  'header' set to TRUE. If this is the case, there is no requirement for
  all the files to have the same ordering of columns (although all
  must have these column names).

Alternatively, the columns of data in the input files can be specified by
the 'cols' argument in the form of a named character vector (e.g;
\code{'cols = c(chr = 1, tag = 2, count = 3, start = 4, end = 5)'} would
cause the function to assume that the first column contains the
chromosome information, the second column contained the tag information,
&c. If 'cols' is specified then information in the header is ignored. If
'cols' is missing and 'header = FALSE' then it is assumed that the
data takes the form described in the example above.

The \code{'tag'} and \code{'count'} columns may optionally be
omitted from either the file column headers or the 'cols' argument. If
the \code{'tag'} column is omitted, then the data will not account for
duplicated sequences when estimating the number of counts in loci. If
the \code{'count'} column is omitted, the \code{'processTags'} function
will assume that the file contains the alignments of each copy of each
sequence tag, rather than an aggregated alignment of each unique
sequence. The unique alignments will be identified and the number of
sequence tags aligning to each position will be calculated.

The \code{replicates} argument should take the form of a vector of
integers such that if and only if the ith library is a replicate of the
jth library then \code{@replicates[i] == @replicates[j]}. In addition,
values in the replicates slot should take values from \code{1:n} where
\code{n} is the number of replicate groups.

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

datadir <- system.file("data", package = "segmentSeq")
libfiles <- dir(datadir, pattern = ".txt", full.names = TRUE)

# Establish the library names and replicate structure.

libnames <- c("SL10", "SL26", "SL32", "SL9")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(libfiles, replicates, libnames, chrlens, chrs = c(">Chr1", ">Chr2"), header = TRUE)

}
\keyword{files}
