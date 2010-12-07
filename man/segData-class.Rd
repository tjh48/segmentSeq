\name{segData-class}
\Rdversion{1.1}
\docType{class}
\alias{segData-class}
\alias{[,segData-method}
\alias{[,segData,ANY,ANY-method}
\alias{dim,segData-method}
\alias{initialize,segData-method}
\alias{show,segData-method}

\title{Class "segData"}
\description{  The \code{segData} class contains data about potential
  segments on the genome containing data about each potential subsegment.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("segData", ...,
  seglens)}. However, more usually they will be created by calling the
\code{\link{processAD}} function.
}
\section{Slots}{
  \describe{
    \item{\code{data}:}{Object of class \code{"matrix"}. Contains the
      number of counts observed for each sample in each potential
      segment.}
    \item{\code{leftData}:}{Object of class \code{"matrix"}. Contains
      the number of counts observed for the region to the left of the
      potential segment.}
    \item{\code{rightData}:}{Object of class \code{"matrix"}. Contains
      the number of counts observed for the region to the right of the
      potential segment.}
    \item{\code{libsizes}:}{Object of class \code{"numeric"}. The
      library sizes for each sample.}
    \item{\code{replicates}:}{Object of class \code{"numeric"}. The
      replicate structure for the samples. This should be a vector of
      consecutive integers starting with 1.}
    \item{\code{priorType}:}{Character string describing the type of prior
  information available in slot \code{'priors'}.}
    \item{\code{priors}:}{Prior parameter information, estimated from the
    data (or otherwise acquired). See Details.}
    \item{\code{segInfo}:}{Object of class \code{"data.frame"}. A data.frame containing the following columns;
    'chr', 'start', 'end', 'leftSpace', 'rightSpace'. See Details.}
  }
}

\section{Details}{
  The \code{@segInfo} slot contains information on each of the potential
  segments; specifically, chromosome, start and end of the segment,
  together with the distance from each segment to the next segment on the
  left and right hand sides. These data are contained in the columns
  'chr', 'start', 'end', 'leftSpace', 'rightSpace' respectively. Each
  row of the \code{@segInfo} slot should correspond to the same row of
  the \code{@data} slot.

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

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, maxgaplen = 500, cl = NULL)

# Estimate prior parameters for the segData object.

}
\keyword{classes}
