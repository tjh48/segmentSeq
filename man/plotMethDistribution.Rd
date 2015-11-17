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
plotMethDistribution(meth, samples, bw = 1e-3, subtract, chrs, centromeres,
add = FALSE, col, ylim = NULL, legend = TRUE, ...)
}
\arguments{
  \item{meth}{
    An object of class \code{\link{alignmentMeth}} containing the
    methylation data.
}
  \item{samples}{
    A numeric vector defining the columns of data in the `meth' object
    from which to estimate proportions of methylation, or a list object
    containing numeric vectors if multiple distributions are to be
    derived from the `meth' object, or a factor in which each level of
    the factor defines a set of columns for the `meth' object. If
    missing, derived from the `@replicates' slot of the `meth' object.
  }
  \item{bw}{
    Gives the bandwidth of the density plots; analagous to the `bw'
    parameter in \code{\link{density}}.
    }
  \item{subtract}{
    A numeric vector giving values to be subtracted from the density of
    methylation. See Details.
  }
  \item{chrs}{
    The names of the chromosomes for which a distribution should be
    plotted. If mising, derived from the `@alignments' slot of the
    `meth' object.
  }
  \item{centromeres}{
    If given, a numeric vector defining the position of the centromeres
    on the chromosomes. These will be then be plotted.
  }
  \item{add}{
    Should the distribution curve be added to an existing plot? Defaults
    to FALSE.
  }
  \item{col}{
    A vector of colours to be used to plot the distributions. If
    missing, generated from \code{\link{rainbow}}.
  }
  \item{ylim}{
    Limits on the y-axis. Defaults to NULL, in which case limits are
    automatically set.
    }
  \item{legend}{
    Should a legend be added to the plot? Defaults to TRUE.
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
