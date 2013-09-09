\name{methData-class}
\docType{class}
\alias{methData}
\alias{methData-class}
\alias{dim,methData-method}
\alias{[,methData-method}
\alias{[,methData,ANY,ANY-method}
\alias{show,methData-method}

\title{Class "methData"}

\description{
  The \code{methData} class is based on the \code{\link{pairedData}}
  class defined in the `baySeq' package, but includes a `coordinates'
  slot giving the coordinates of genomic loci and a `locLikelihoods'
  slot which contains the estimated likelihoods that each annotated
  region is a locus in each replicate group and a coordinates structure
  giving the locations of the loci.
}

\section{Slots}{
  \describe{
    \item{\code{locLikelihoods}:}{Object of class \code{"matrix"}
      describing estimated likelihoods that each region defined in
      `coordinates' is a locus in each replicate group.}
    \item{\code{coordinates}:}{Object of class \code{"GRanges"} defining
      the coordinates of the genomic loci.}
    \item{\code{data}:}{Object of class \code{"matrix"} defining the
      number of methylated cytosines observed for each locus defined in
      `coordinates'}
    \item{\code{data}:}{Object of class \code{"matrix"} defining the
      number of un-methylated cytosines observed for each locus defined in
      `coordinates'}
    \item{\code{replicates}:}{Object of class \code{"factor"} defining
      the replicate structure of the data.}
    \item{\code{libsizes}:}{Object of class \code{"numeric"} describing
      the library size (scaling factor) for each sample.}
    \item{\code{groups}:}{Object of class \code{"list"} defing the group
      (model) structure of the data (see the
      \code{\link[baySeq:baySeq-package]{baySeq} package})}.
    \item{\code{annotation}:}{Object of class \code{"data.frame"} giving
      any additional annotation information for each locus.}
    \item{\code{priorType}:}{Object of class \code{"character"}
      describing the type of prior information available in slot \code{'priors'}.}
    \item{\code{priors}:}{Object of class \code{"list"} defing the prior
      parameter information. Calculated by
      the \link[baySeq:baySeq-package]{baySeq} package.}
    \item{\code{posteriors}:}{Object of class \code{"matrix"} giving the
      estimated posterior likelihoods for each replicate group. Calculated by
      the \link[baySeq:baySeq-package]{baySeq} package.}
    \item{\code{nullPosts}:}{Object of class \code{"numeric"} which, if
      calculated, defines the posterior likelihoods for the data having no
      true expression of any kind. Calculated by
      the \link[baySeq:baySeq-package]{baySeq} package.}
    \item{\code{estProps}:}{Object of class \code{"numeric"} giving the
      estimated proportion of tags belonnging to each group. Calculated by
      the \link[baySeq:baySeq-package]{baySeq} package.}
    \item{\code{seglens}:}{Object of class \code{"matrix"} defining the
      lengths of each segment containing the counts described in the
      \code{'data'} slot. May be initialised with a vector, or ignored altogether.}
  }
}

\section{Extends}{
  Class \code{"\linkS4class{pairedData}"}, directly.
}

\section{Methods}{
  Methods `new', `dim', `[' and `show' have been defined for this class.
}

\author{Thomas J. Hardcastle}

\keyword{classes}
