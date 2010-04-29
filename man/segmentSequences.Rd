\name{segmentSequences}
\alias{segmentSequences}

\title{Takes the 'segData' object and uses it to define segments on the genome.
}
\description{
This function takes the 'segData' object and produces a set of segments
on the genome based on the data contained within it.
}
\usage{
segmentSequences(sDP, pcut = 0.5, priorDE =
1e-2, verbose =
TRUE, ..., cl)
}

\arguments{
  \item{sDP}{
    A \code{\linkS4class{segData}} object.
}
  \item{pcut}{
    The maximum acceptable likelihood that a potential segment is
    similar to background or to the regions to either the left and right
    of the potential segment. See Details.
  }
  \item{priorDE}{Prior likelihood of similarity.}
  \item{verbose}{
    Should processing information be displayed? Defaults to TRUE.
  }

\item{...}{Additional parameters to be passed to the
  \code{\link[baySeq:getLikelihoods]{getLikelihoods}} function of the
  'baySeq' package.}

\item{cl}{A SNOW cluster object, or NULL. See Details.}

}
\details{
This function takes each potential segment defined by the
  \code{\linkS4class{segData}} object and evaluates the likelihood that it is
  similar to either background or the empty region to the left and right
  of the potential segment in all replicate groups. See Hardcastle &
  Kelly (2010) for more details on how this likelihood is evaluated.

  The potential segments are then ranked by increasing likelihood of
  similarity. Any potential segment with a likelihood of similarity
  greater than \code{'pcut'} is discarded.

  If \code{'estimatePriors = TRUE'} then an attempt will be made to
  empirically re-estimate the prior likelihoods of similarity from the
  data. This may improve the segmentation slightly at some computational
  cost.

  There is then a filtration step (\code{\link{filterSegments}}).
  The segment with the lowest likelihood of similarity is kept, and any
  segments that have overlap with this segment are discarded. The
  segment with the next lowest likelihood of similarity is then kept,
  and any segments that have overlap with this segment are
  discarded. This process continues until we have a set of
  non-overlapping segments.
}
\value{
  A \code{\link[baySeq:countData]{countData}} object, containing count
  information on all the segments discovered.
}
\references{
  Hardcastle T.J., and Kelly, K.A. (2010). Genome Segmentation From
  High-Throughput Sequencing Data. In submission.

  Hardcastle T.J., and Kelly, K.A. (2010). Empirical Bayesian
  Methods For Identifying Patterns of Differential Expression in Count
  Data. In submission.
}

\author{
Thomas J. Hardcastle
}

\seealso{
  \code{\link{getPriors}}, a function for establishing prior parameters
  used to estimate posterior likelihoods.
  \code{\link{plotGenome}}, a function for plotting the alignment of
  tags to the genome (together with the segments defined by this
  function).
  \code{\link[baySeq:baySeq-package]{baySeq}}, a package for discovering
  differential expression in \code{\link[baySeq:countData]{countData}} objects.
      
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

# Use the segData object to produce a segmentation of the genome.

segD <- segmentSequences(sDP, pcut = 0.1, cl = NULL)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}
\keyword{classif}
