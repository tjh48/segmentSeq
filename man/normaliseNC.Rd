% modification on git from copied files
\name{normaliseNC}
\alias{normaliseNC}
\title{
  A function providing adjustment of cytosine methylated/unmethylated
  counts based on a nonconversion rate.
}
\description{
  This function adjusts the observed cytosine methylated/unmethylated
  counts at each cytosine site based on the reported nonconversion rates
  for each samples.
}
\usage{
normaliseNC(mD, nonconversion)
}
\arguments{
  \item{mD}{
    Either an \code{\link{alignmentMeth}} or \code{\link{segMeth}}
    object, or a \code{\link{lociData}} object (for which
    nonconversion must be explicitly supplied).
}
  \item{nonconversion}{
    A vector defining nonconversion rates for each sample, required if a
    \code{\link{lociData}} object is supplied in `mD' and ignored otherwise.
}
}
\details{
This function operates by estimating the expected number of
unconverted cytosines at each site and subtracting this from the
reported methylated cytosines and adding to the reported unmethylated
cytosines. It should not be used on data that will be analysed in a way
that accounts for nonconversion; e.g., using the `bbNCDist'
densityFunction object. 
}
\value{
A modified version of the object supplied as `mD'.
}
\references{
Hardcastle T.J. Discovery of methylation loci and analyses of
differential methylation from replicated high-throughput sequencing
data. bioRxiv (http://dx.doi.org/10.1101/021436)
}
\author{
Thomas J. Hardacastle
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

mD <- normaliseNC(mD)
}
\keyword{manip}
