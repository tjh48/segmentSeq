% modification on git from copied files
\name{segData-class}
\Rdversion{1.1}
\docType{class}
\alias{segData}
\alias{segData-class}
\alias{[,segData-method}
\alias{[,segData,ANY,ANY-method}
\alias{[,segData,ANY-method}
\alias{dim,segData-method}
\alias{initialize,segData-method}
\alias{show,segData-method}

\title{Class "segData"}
\description{  The \code{segData} class inherits from the
  \code{segClass} class and contains data about potential
  segments on the genome, together with counts for each of those segments.}
\section{Objects from the class}{
Objects can be created by calls of the form \code{new("segData", ...,
  seglens)}. However, more usually they will be created by calling the
\code{\link{processAD}} function.
}
\section{Slots}{
  \describe{
    \item{\code{coordinates}:}{A \code{\link[GenomicRanges]{GRanges}}
      object defining the coordinates of the segments.}
    
    \item{\code{replicates}:}{Object of class \code{"factor"}. The
      replicate structure for the samples.}
    \item{\code{locLikelihoods}:}{Object of class \code{"DataFrame"}
      describing estimated likelihoods that each region defined in
      `coordinates' is a locus in each replicate group.}
    \item{\code{data}:}{Object of class \code{matrix}. Contains the
      number of counts observed for each sample in each potential
      segment.}
    \item{\code{libsizes}:}{Object of class \code{"numeric"}. The
      library sizes for each sample.}
  }
}

\section{Details}{
  The \code{@coordinates} slot contains information on each of the potential
  segments; specifically, chromosome, start and end of the segment,
  together. Each row of the \code{@coordinates} slot should correspond to
  the same row of the \code{@data} slot.

  In almost all cases objects of this class should be produced by the
  \code{\link{processAD}} function.
}

\section{Methods}{
  Methods 'new', 'dim', '[' and 'show' have been defined for this class.
}

\author{Thomas J. Hardcastle}

\seealso{
\code{\link{processAD}}, the function that will most often be used to
create objects of this class.
\code{\link{classifySeg}}, an empirical Bayesian method for defining a
segmentation based on a segData object.

}
\examples{
# Define the files containing sample information.

datadir <- system.file("extdata", package = "segmentSeq")
libfiles <- c("SL9.txt", "SL10.txt", "SL26.txt", "SL32.txt")

# Establish the library names and replicate structure.

libnames <- c("SL9", "SL10", "SL26", "SL32")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- readGeneric(file = libfiles, dir = datadir, replicates =
replicates, libnames = libnames)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, gap = 100, cl = NULL)

}
\keyword{classes}
