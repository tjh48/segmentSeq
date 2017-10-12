% modification on git from copied files
\name{readMeths}
\alias{readMeths}
\title{
A function for reading data from the YAMA methylation aligner (or
similarly parsed data) from which to identify methylation loci and/or
differentially methylated regions.
}
\description{
This function takes as input a set of files that describe the number of
times a set of cytosines are observed to be methylated or unmethylated
in some high-throughput sequencing data. It merges the data from these
files into  an object of `\code{\link{alignmentMeth}}' class which can
then be further processed to identify methylation loci.
}
\usage{
readMeths(files, dir = ".", libnames, replicates, nonconversion, chrs)
}
\arguments{
  \item{files}{
    A character vector defining the file names of the alignment files to
    be read in.
}
  \item{dir}{
    The directory in which the files are located.
}
  \item{libnames}{
    A character vector giving the names of the samples to be read in.
}
  \item{replicates}{
    A vector defining the replicate structure of the data. The `i'th and
    `j'th libraries are treated as replicates if and only if
    replicates[i] == replicates[j].
}
  \item{nonconversion}{
    A numeric vector (all members should lie between 0 and 1) defining the
    non-conversion rate of each library. See
    \code{\link{alignmentMeth-class}} for details.
  }
  \item{chrs}{
    An (optional) character vector giving the names of the chromosomes
    to be read from the files. If ommitted, all chromosomes will be read
    in.
  }
}

\value{
  An object of class \code{\link{alignmentMeth}}.
}
\author{
Thomas J. Hardcastle
}

\seealso{
\code{\link{alignmentMeth-class}}.
}
\examples{
datadir <- system.file("extdata", package = "segmentSeq")
files <- c("short_18B_C24_C24_trim.fastq_CG_methCalls.gz",
"short_Sample_17A_trimmed.fastq_CG_methCalls.gz",
"short_13_C24_col_trim.fastq_CG_methCalls.gz",
"short_Sample_28_trimmed.fastq_CG_methCalls.gz")

mD <- readMeths(files = files, dir = datadir,
libnames = c("A1", "A2", "B1", "B2"), replicates = c("A","A","B","B"),
nonconversion = c(0.004777, 0.005903, 0.016514, 0.006134))
}
\keyword{files}
