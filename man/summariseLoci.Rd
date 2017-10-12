% modification on git from copied files
\name{summariseLoci}
\alias{summariseLoci}
\title{
Summarise the expected number of loci in a `lociData' object.
}
\description{
  Summarises the expected number of loci, either in toto or on a per
  replicate group basis.
}
\usage{
summariseLoci(cD, perReplicate = TRUE)
}
\arguments{
  \item{cD}{
    A `lociData' object with calculated values in the `@lociLikelihoods' slot.
}
  \item{perReplicate}{
    Should the expectation be calculated on a per replicate group basis,
    or the total number of loci identified in the dataset?
}
}
\value{
  A numeric vector summarising the expectated number of loci in the cD object.
}
\author{
Thomas J. Hardcastle
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{mainip}

