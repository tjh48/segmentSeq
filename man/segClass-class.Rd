\name{segClass-class}
\Rdversion{1.1}
\docType{class}
\alias{segClass}
\alias{segClass-class}
\alias{[,segClass-method}
\alias{[,segClass,ANY,ANY-method}
\alias{dim,segClass-method}
\alias{initialize,segClass-method}
\alias{show,segClass-method}

\title{Class "segClass"}
\description{  The \code{segClass} class contains data about potential
  segments on the genome.}
\section{Objects from the class}{
Objects can be created by calls of the form \code{new("segClass", ...,
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
\code{\link{segData}}, which inherits from this class.
\code{\link{segMeth}}, which inherits from this class.
}

\keyword{classes}