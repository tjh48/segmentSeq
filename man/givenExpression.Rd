\name{givenExpression}
\alias{givenExpression}

\title{
  Adjusts posterior likelihoods of differential expression by the
  likelihood that a locus is expressed.
}
\description{For certain types of differential behaviour, the likelihood
  of difference can be high even for regions for which the likelihood
  that the region represents a true locus is low in all replicate
  groups. This function corrects the likelihood of differential
  behaviour by the likelihood that the locus is expressed in at least
  one replicate group.
}
\usage{
givenExpression(cD)
}

\arguments{
  \item{cD}{A lociData object with computed posteriors and
  locLikelihoods.}
}
\value{
A lociData object with updated posteriors.
}
\author{
  Thomas J. Hardcastle
}

\keyword{mainip }% use one of  RShowDoc("KEYWORDS")
