\name{getPriors}
\alias{getPriors}

\title{
Finds a set of priors for a 'segData' object.
}
\description{
This function creates a random selection of non-overlapping segments
that can be used to estimate prior parameters for the 'segData'
object.
}
\usage{
getPriors(sD, type = "Pois", verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sD}{
    A \code{\linkS4class{segData}} object.
}
  \item{type}{
    A string describing the type of priors to be estimated. Currently
    only "Pois" (Poisson-Gamma priors) and "NB" (Negative Binomial) are supported.
  }

  \item{verbose}{
    Should processing information be displayed? Defaults to TRUE.
  }
  
  \item{\dots}{
    Additional arguments to be passed to one of the
    \code{\link[baySeq:getPriors]{getPriors}} functions. See Details.
}
}
\details{
This function takes a random sample of non-overlapping potential
subsegments from the 'segData' object and uses these to construct a
\code{\link[baySeq:countData]{countData}} object which is then passed to
one of the \code{\link[baySeq:getPriors]{getPriors}} functions belonging
to the 'baySeq' package. Which function is specified depends on the
string given in the \code{priorType} argument; currently only
\code{priorType = "Pois"} and \code{priorType = "NB"} are supported.
Additional arguments can be passed to whichever function is being used
via the \code{'...'} argument.
}
\value{
  A \code{\linkS4class{segData}} object with a \code{\linkS4class{priorData}}
  structure in the \code{'@priors'} slot.
}
\references{
Hardcastle T.J., and Kelly, K.A. (2010). Genome Segmentation from
High-Throughput Sequencing Data. In submission.
}
\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\linkS4class{segData}}, \code{\link[baySeq:getPriors]{getPriors}}
}
\examples{

# Define the chromosome lengths for the genome of interest.

chrlens <- c(2e6, 1e6)

# Define the files containing sample information.

datadir <- system.file("data", package = "segmentSeq")
libfiles <- dir(datadir, pattern = ".txt", full.names = TRUE)

# Establish the library names and replicate structure.

libnames <- c("SL10", "SL26", "SL32", "SL9")
replicates <- c(1,1,2,2)

# Process the files to produce an 'alignmentData' object.

alignData <- processTags(libfiles, replicates, libnames, chrlens, chrs = c(">Chr1", ">Chr2"), header = TRUE)

# Process the alignmentData object to produce a 'segData' object.

sD <- processAD(alignData, maxgaplen = 500, cl = NULL)

# Estimate prior parameters for the segData object.

sDP <- getPriors(sD, type = "Pois", samplesize = 100, perSE = 0.1, maxit
= 1000, cl = NULL)


}
\keyword{distribution}
\keyword{models}
