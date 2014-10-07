\name{averageProfiles}
\alias{averageProfiles}
\alias{plotAverageProfile}

\title{
Computes and plots the average distribution of aligned reads (taken from an
alignmentData object) or methylation (taken from an
alignmentMeth object) over a set of coordinates (and optionally the
surrounding regions).
}
\description{
Given an alignmentData or alignmentMeth object and a set of coordinates,
plots the average distribution of coverage/methylation over those
coordinates. The plotted distributions can be split up into different
sample groups by the user. 
}
\usage{
averageProfiles(mD, samples, coordinates, cuts, maxcuts = 200, bw =
5000, surrounding = 0, add = FALSE, col, ylim, meanOnly = TRUE, plot =
TRUE, ...)
plotAverageProfile(position, profiles, col, surrounding, ylim, add =
FALSE, meanOnly = TRUE, legend = TRUE, titles, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mD}{
    The \code{\linkS4class{alignmentData}} or
    \code{\linkS4class{alignmentMeth}} object defining the 
    coverage or methylation on the genome.
}
  \item{samples}{
    A factor or list defining the different groups of samples for which
    to plot different distributions. If a list, each member should
    consist of integer values defining the column numbers of the 'mD'
    object. If missing, will use the mD@replicates value.
}
  \item{coordinates}{
    A \code{\linkS4class{GRanges}} object defining the coordinates of
    interest (e.g. genic regions).
}
  \item{cuts}{
    Optionally, the number of subdivisions to cut the coordinates in
    when calculating the average coverage/methylation density.
  }
  \item{maxcuts}{
    The maximum number of subdivisions permitted when calculating the
    average coverage/methylation density.
  }
  \item{bw}{
    If 'cuts' is missing, this factor divides the product of the length of the
    'coordinates' object and the median width of the coordinates to
    defines the number of cuts (minimum twenty cuts).
}
  \item{surrounding}{
    If non-zero, the size of the region up- and down-stream of the given
    coordinates for which average coverage/methylation should also be calculated.
}
  \item{add}{
    If TRUE, the plotted distribution will be added to the current plot.
}
  \item{col}{
    If given, a vector of colours for each of the groups defined in
    'samples'. Defaults to 'rainbow(length(samples))'.
  }
  \item{ylim}{
    See `ylim' option for plot. If missing, will be calculated from
    data.
  }
  \item{meanOnly}{
    If TRUE, the mean methylation profile for each member of the
    'samples' parameter is plotted on a single graph. If FALSE, every
    5th percentile is plotted for each member of the sample parameters,
    each on a separate graph.
  }
  \item{plot}{Should the profile be plotted? Defaults to TRUE.}
  \item{position}{
    A vector describing the position of each point to be plotted. Take
    from the `$position' element in the list object returned by
    `averageProfiles'.
  }
  \item{profiles}{
    A matrix describing the profiles to be plotted. Take
    from the `$profiles element in the list object returned by
    `averageProfiles'.
  }
  \item{legend}{
    If TRUE, a legend describing the samples is included on the plot.
  }
  \item{titles}{
    If given, and `meanOnly = FALSE', a vector of titles for the
    quantile plots.
    }  
  \item{\dots}{
    Additional arguments to be passed to the 'plotAverageProfile'
    function, and hence to the 'plot' or 'lines' methods.    
  }
}
\value{
  Invisibly, a list containing the coordinates of the lines plotted.
}

\author{
  Thomas J. Hardcastle
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
