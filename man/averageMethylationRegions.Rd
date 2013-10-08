\name{averageMethylationRegions}
\alias{averageMethylationRegions}

\title{
Plots the average distribution of methylation (taken from an
alignmentMeth object) over a set of coordinates (and optionally a
surrounding region).
}
\description{
Given an alignmentMeth object and a set of coordinates, plots the
average distribution of methylation over those coordinates. The
distributions can be split up into different sample groups by the user.
}
\usage{
averageMethylationRegions(mD, samples, coordinates, cuts, bw = 0.01, surrounding = 0, add = FALSE, col, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mD}{
    The \code{\linkS4class{alignmentMeth}} object object defining the
    methylation at each sequenced base
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
    when calculating the average methylation density.
}
  \item{bw}{
    If 'cuts' is missing, this factor, multiplied by the length of the
    'coordinates' vector, defines the number of cuts (minimum twenty cuts).
}
  \item{surrounding}{
    If non-zero, the size of the region up- and down-stream of the given
    coordinates for which average methylation should also be calculated.
}
  \item{add}{
    If TRUE, the plotted distribution will be added to the current plot.
}
  \item{col}{
    If given, a vector of colours for each of the groups defined in
    'samples'. Defaults to 'rainbow(length(samples))'.
}
  \item{\dots}{
    Additional arguments to be passed to the 'plot' or 'lines' functions.
}
}
\value{
  Invisibly, a list containing the coordinates of the lines plotted.
}

\author{
  Thomas J. Hardcastle
}

\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}
