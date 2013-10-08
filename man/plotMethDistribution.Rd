\name{plotMethDistribution}
\alias{plotMethDistribution}

\title{
Plots the distribution of methylation on the genome.
}
\description{
Plots the distribution of methylation (as defined in an
\code{\link{alignmentMeth}} object upon the genome.
}
\usage{
plotMethDistribution(meth, samples, bw = 1e-3, subtract, centromeres,
add = FALSE, ...)
}
\arguments{
  \item{meth}{
    An object of class \code{\link{alignmentMeth}} containing the
    methylation data.
}
  \item{samples}{
    A numeric vector defining the columns of data in the `meth' object
    from which to estimate proportions of methylation.
}
  \item{subtract}{
    A numeric vector giving values to be subtracted from the density of
    methylation. See Details.
  }
  \item{bw}{A scaling factor for the estimation of density of methylation.}
  \item{centromeres}{
    If given, the location of the centromeres will be plotted.
  }
  \item{add}{
    Should the distribution curve be added to an existing plot? Defaults
    to FALSE.
}
  \item{\dots}{
    Any additional parameters to be passed to \code{\link{plot}}.
}
}
\details{
  The function returns the density of methylation calculated.
  This can be used in further plots as the 'subtract' parameter, which
  allows one methylation profile to be subtracted from another.
}
\value{
An object of class \code{\link{density}} describing the plotted distribution.
}
\author{
Thomas J. Hardcastle
}

\seealso{
\code{\link{alignmentMeth}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
