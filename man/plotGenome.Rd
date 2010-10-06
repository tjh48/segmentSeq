\name{plotGenome}
\alias{plotGenome}
\title{Plots the alignment of sequence tags on the genome given an
  'aligmentData' object and (optionally) a set of segments found.
}
\description{
  Plots the data from an \code{alignmentData} object for a given set of
  samples. Can optionally include in the plot the annotation data from a
  \code{countData} object containing segment information.
}
\usage{
plotGenome(aD, sD, chr = 1, limits = c(0, 10^4), samples = NULL,
plotType = "pileup", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{sD}{
    A \code{\linkS4class{countData}} object (produced by the
    \code{\link{similaritySeg}} function and therefore) containing
    appropriate annotation information. Can be omitted if this
    annotation is not known/required.
  }
  \item{chr}{
    The name of the chromosome (translated into \code{'character'} type
    if given in any other form) to be plotted. Should correspond to a
    chromosome name in the \code{alignmentData} object.
  }
  \item{limits}{
    The start and end point of the region to be plotted.
  }
  \item{samples}{
    The sample numbers of the samples to be plotted. If NULL, plots all samples.
  }
  \item{plotType}{
    The manner in which the plot is created. Currently only
    \code{'plotType = pileup'} is supported.
  }
  \item{...}{
    Any additional graphical parameters for passing to \code{plot}.
  }
}

\value{
  Plotting function.
}
\author{
Thomas J. Hardcastle
}

\seealso{
\code{\linkS4class{alignmentData}}, \code{\link{similaritySeg}}
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

# Plot the alignments to the genome on chromosome 1 between bases 1 and 10000

plotGenome(alignData, chr = ">Chr1", limits = c(1, 1e5))

}

\keyword{hplot}
