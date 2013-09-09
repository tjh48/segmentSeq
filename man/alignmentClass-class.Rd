\name{alignmentClass-class}
\Rdversion{1.1}
\docType{class}
\alias{alignmentClass-class}
\alias{alignmentClass}
\alias{[,alignmentClass-method}
\alias{[,alignmentClass,ANY,ANY-method}
\alias{dim,alignmentClass-method}
\alias{initialize,alignmentClass-method}
\alias{show,alignmentClass-method}
\alias{cbind,alignmentClass-method}

\title{Class "alignmentClass"}
\description{The \code{alignmentClass} class records information about
  a set of alignments of high-throughput sequencing data to a
  genome.}

\section{Slots}{
  \describe{
    \item{\code{alignments}:}{Object of class \code{"GRanges"}.
      Stores information about the alignments. See Details.}
    \item{\code{libnames}:}{Object of class \code{"character"}. The
      names of the libraries for which alignment data exists.}
    \item{\code{replicates}:}{Object of class \code{"factor"}.
      Replicate information for each of the libraries. See Details.}
  }
}

\section{Details}{
  The \code{alignments} slot is a \code{GRanges} object defining the location of aligned
  objects to a reference genome.

  The \code{replicates} slot is a vector of factors such that the ith
  sample is a replicate of the jth sample if and only if \code{@replicates[i] ==
    @replicates[j]}.

  The \code{libnames} slot is a vector defining the names of the
  libraries described by the object.
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "alignmentClass")}: ... }
    \item{dim}{\code{signature(x = "alignmentClass")}: ... }
    \item{initialize}{\code{signature(.Object = "alignmentClass")}: ... }
    \item{show}{\code{signature(object = "alignmentClass")}: ... }
	 }
}

\author{Thomas J. Hardcastle}

\seealso{
  \code{\link{alignmentData}}
  \code{\link{alignmentMeth}}
}

\keyword{classes}

