\name{thresholdFinder}
\alias{thresholdFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determines threshold for the proportion of methylation at which a
methylation locus may be identified.
}
\description{
    This function offers a variety of methods for the analysis of
    methylation data to determine a suitable threshold for the
    proportion of methylation at which to distinguish a methylation
    locus from a non-methylated locus.
}
\usage{
thresholdFinder(method, aM, subset, minprop = 0.05, bootstrap = 100,
verbose = FALSE, cl = NULL, processAD.args = list(), heuristicSeg.args = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    Character string defining method to use for threshold
    estimation. Available options are 'varsum', 'minden', 'beta' and
    'abc'. See Details.
}
  \item{aM}{
    An \code{\link{alignmentMeth}} object containing observed
    methylation counts.
  }
  \item{subset}{Numeric vector defining a subset on \code{aM} object for
    use in threshold estimation.}
  \item{minprop}{
    For 'minden' method, a minimum proportion permitted for choice of threshold.
  }
  \item{bootstrap}{
    The maximum number of bootstraps to be permitted in estimating a
    threshold. Defaults to 100. See Details.
    }
  \item{verbose}{
    Verbose reporting. Defaults to FALSE.
  }
  \item{cl}{A cluster object, or NULL. Defaults to NULL.}

  \item{processAD.args}{Arguments to be passed to processAD function if
    bootstrapping.}
  \item{heuristicSeg.args}{Arguments to be passed to heuristicSeg function if
    bootstrapping.}

}
\details{
  This function operates on the data observed within each replicate
    group, and then takes the mean of the thresholds calculated for each
    group.
  
  Methods currently available for threshold estimation are 'varsum',
    'minden', 'beta' and 'abc'. The 'varsum' method attempts to split
    the vector of proportions of methylation observed at each cytosine
    into two sets of minimal total variance. The 'minden' method
    finds the minimum point on a smoothed kernel density of the
    proportions of methylation. The 'beta' method estimates for each
    cytosine a posterior distribution on proportions of methylation
    based on the beta-binomial conjugacy, takes the average of these
    distributions and finds the minumum. The 'abc' method performs like
    the beta method, but estimates the posterior distribution through
    approximate Bayesian computation.

    Bootstrapping uses the estimated threshold to define
    loci. Based on the defined loci, cytosines are then only included in
    a re-estimation of the thresholds if they are identified as
    belonging to an expressed locus within the current replicate group,
    or if they are not expressed in any replicate group.
}
\value{
A numeric value defining a threshold on methylation.
}
\author{
Thomas J. Hardcastle
}

\seealso{
\code{\link{heuristicSeg}}
}
\examples{
datadir <- system.file("extdata", package = "segmentSeq")
files <- c("short_18B_C24_C24_trim.fastq_CG_methCalls.gz",
"short_Sample_17A_trimmed.fastq_CG_methCalls.gz",
"short_13_C24_col_trim.fastq_CG_methCalls.gz",
"short_Sample_28_trimmed.fastq_CG_methCalls.gz")

mD <- readMeths(files = files, dir = datadir,
libnames = c("A1", "A2", "B1", "B2"), replicates = c("A","A","B","B"),
nonconversion = c(0.004777, 0.005903, 0.016514, 0.006134))

thresholdFinder("beta", mD, cl = NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}% use one of  RShowDoc("KEYWORDS")
