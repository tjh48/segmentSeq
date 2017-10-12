% modification on git from copied files
\name{plotMeth}
\alias{plotMeth}

\title{
Plots a map of cytosine methylation (and optionally, methylation loci).
}
\description{
This function takes an \code{\link{alignmentMeth}} object and plots the
observed levels of methylation within some defined region on the
genome. If a \code{\link{methData}} object is also supplied, loci of
methylation will also be shown.
}
\usage{
plotMeth(aM, loci, chr, limits, samples, showNumber = TRUE, rgb = c(1, 0, 0), angle = 45, cap, add = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{aM}{
    An \code{\linkS4class{alignmentMeth}}.
  }
  \item{loci}{
    %%     ~~Describe \code{loci} here~~
  }
  \item{chr}{
    The name of the chromosome to be plotted. Should correspond to a
    chromosome name in the \code{alignmentMeth} object. 
  }
  \item{limits}{
    The start and end point of the region to be plotted.
  }
  \item{samples}{
    The sample numbers of the samples to be plotted. If NULL, plots all samples.
}
  \item{showNumber}{
    Should the row number of each segment be shown?
  }
  \item{rgb}{
    The rgb code (\code{\link{rgb}}) with which to colour the loci.
}
  \item{angle}{
    The angle at which loci are shaded (see \code{\link{rect}}).
}
  \item{cap}{
    Caps the maximum level of coverage shown on the plot; thus, if a
    base has been sequenced at a level greater than the cap, the data
    for that base will be shown as if it has a coverage of cap.
}
  \item{add}{
    If TRUE, adds the plot of methylation level to the current
    plot. Defaults to FALSE.
}
}
\value{
  Plotting function.
}
\author{
  Thomas J. Hardcastle
}
\seealso{
\code{\linkS4class{alignmentMeth}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
