\name{plotGenome}
\alias{plotGenome}
\title{Plots the alignment of sequence tags on the genome given an
  `aligmentData' object and (optionally) a set of segments found.
}
\description{
  Plots the data from an \code{alignmentData} object for a given set of
  samples. Can optionally include in the plot the annotation data from a
  \code{lociData} object containing segment information.
}
\usage{
plotGenome(aD, locData, chr = 1, limits = c(0, 1e4), samples = NULL,
plotType = "pileup", plotDuplicated = FALSE, density = 0, showNumber =
TRUE, logScale = FALSE, cap = Inf, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{aD}{
    An \code{\linkS4class{alignmentData}} object.
  }
  \item{locData}{
    A \code{\linkS4class{lociData}} object (produced by the
    \code{\link{heuristicSeg}} or \code{\link{classifySeg}} function and therefore) containing
    appropriate annotation information. Can be omitted if this
    annotation is not known/required.
  }
  \item{chr}{
    The name of the chromosome to be plotted. Should correspond to a
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
    \code{plotType = pileup} is recommended.
  }
  \item{plotDuplicated}{
    If TRUE, then any duplicated sequence tags (i.e., sequence tags that
    match to multiple places in the genome) in the `aD' object will be plotted on
    a negative scale for each sample. Defaults to FALSE (recommended).}
  \item{density}{The density of the shading lines to be used in plotting
    each segment.}
  \item{showNumber}{Should the row number of each segment be shown?}
  \item{logScale}{Should a log scale be used for the number of sequence
    tags found at each base?}
  \item{cap}{A numeric value defining a cap on the maximum number of
    reads to be plotted at any one point. Useful if a large number of
    reads at one location prevent a clear signal being seen elsewhere.}
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
\code{\linkS4class{alignmentData}}, \code{\link{heuristicSeg}}, \code{\link{classifySeg}}
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
chrlens, gap = 100)

# Plot the alignments to the genome on chromosome 1 between bases 1 and 10000

plotGenome(alignData, chr = ">Chr1", limits = c(1, 1e5))

}

\keyword{hplot}
