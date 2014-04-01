\name{selectLoci}
\alias{selectLoci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filters a `lociData' object based on given selection criteria.
}
\description{
  Selects loci from a `lociData' object based on likelihood, false
  discovery rate or familywise error rate for downstream processing.
}
\usage{
selectLoci(cD, likelihood, FDR, FWER, perReplicate = TRUE)
}
\arguments{
  \item{cD}{
    The \code{\link{lociData}} object to be filtered.
  }
\item{likelihood}{
  If provided, all loci with a likelihood greater than this criterion
  will be selected.
}
\item{FDR}{
  If provided (and likelihood is not provided), the maximal set of
  loci which controls the FDR at this level is selected.
}
\item{FWER}{
  If provided (and likelihood and FDR are not provided), the maximal set of
  loci which controls the FWER at this level is selected.
}
\item{perReplicate}{
  If TRUE, selection of loci is done on a replicate by replicate
  basis. If FALSE, selection will be done on the likelihood that the
  locus represents a true locus in at least one replicate group.
}
}
\value{
  A \code{\link{lociData}} object.
}
\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{lociLikelihoods}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
