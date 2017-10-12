% modification on git from copied files
\name{mergeMethSegs}
\alias{mergeMethSegs}

\title{
Merges neighbouring methylation loci with the same pattern of expression.
}
\description{
Within a region of cytosine methylation, there may be some cytosines
which show no evidence of methylation. The presence of these cytosines
may lead to the region being split into multiple methylation loci. This
function merges neighbouring loci if the pattern of expression is the
same in each locus, and if they are not separated by too great a genomic
distance.
}
\usage{
mergeMethSegs(segs, aD, gap, cl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{
    A \code{\link{methData}} object defining the loci to be merged.
  }
  \item{aD}{
    An \code{\link{alignmentMeth}} object from which the loci have been derived.
  }
  \item{gap}{
    The maximum gap below which neighbouring loci may be merged.
  }
  \item{cl}{
    A cluster object, or NULL.
  }
}
\value{
An object of type \code{\link{methData}}.
}
\author{
  Thomas J. Hardcastle
}

\seealso{
  \code{\link{methData}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
