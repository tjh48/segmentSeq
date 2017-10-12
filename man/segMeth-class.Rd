% modification on git from copied files
\name{segMeth-class}
\Rdversion{1.1}
\docType{class}
\alias{segMeth}
\alias{segMeth-class}
\alias{[,segMeth-method}
\alias{[,segMeth,ANY,ANY-method}
\alias{[,segMeth,ANY-method}
\alias{dim,segMeth-method}
\alias{initialize,segMeth-method}
\alias{show,segMeth-method}

\title{Class "segMeth"}
\description{  The \code{segMeth} class inherits from the
  \code{segClass} class and contains data about potential
  segments on the genome, together with counts for each of those segments.}
\section{Objects from the class}{
Objects can be created by calls of the form \code{new("segMeth", ...,
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
    \item{\code{Cs}:}{Object of class \code{matrix}. Contains the
      number of methylated cytosines (which are sequenced as a `C')
      observed for each sample in each potential segment.}
    \item{\code{Ts}:}{Object of class \code{matrix}. Contains the
      number of unmethylated cytosines (which are sequenced as a `T')
      observed for each sample in each potential segment.}
    \item{\code{nonconversion}:}{Object of class \code{"numeric"}. The
      (estimated) nonconversion rate (see Details) for each of the libraries.}
  }
}

\section{Details}{
  The \code{@coordinates} slot contains information on each of the potential
  segments; specifically, chromosome, start and end of the segment,
  together. Each row of the \code{@coordinates} slot should correspond to
  the same row of the \code{@C} and \code{@T} slots.

  The nonconversion slot is an estimate of the rate (for each library)
  at which an unmethylated cytosine has failed to be converted by sodium
  bisulphite treatment into thymine, and is thus recorded (incorrectly) as
  methylated. In some cases, this can be estimated from considering
  observed methylation rates on regions known to be unmethylated (e.g.,
  chloroplasts) or by introducing unmethylated control sequences.

  
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
\code{\link{segClass}}, from which this class inherits.
}

\keyword{classes}
