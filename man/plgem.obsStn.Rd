\name{plgem.obsStn}
\alias{plgem.obsStn}
\alias{plgem.resampledStn}
\title{Computation of Observed and Resampled PLGEM-STN statistics}
\description{
These functions compute the observed and resampled signal to noise ratio (STN) statistic
using PLGEM fitting parameters (obtained using the function \sQuote{plgem.fit}) to detect differential
expression in a microarray exprSet \sQuote{data}.
}
\usage{plgem.obsStn(data, plgemFit, baseline.condition = 1, verbose = FALSE)
plgem.resampledStn(data, plgemFit, baseline.condition = 1,
iterations = "automatic", verbose=FALSE)
}
\arguments{
  \item{data}{an object of class \sQuote{exprSet} with a \sQuote{conditionName} covariate, see details.}
  \item{plgemFit}{list; the output of \sQuote{plgem.fit}.}
  \item{baseline.condition}{the condition to be treated as the baseline.}
  \item{verbose}{logical; if TRUE, comments are printed out while running.}
  \item{iterations}{number of iterations for the resampling step; if "automatic" it is automatically determined.}
}
\details{
The exprSet \sQuote{data} must have at least two conditions; by default the first one is considered the baseline.
\sQuote{data} must have a phenoData slot with a covariate called \sQuote{conditionName}. The values of this covariate
must be sample labels, that have to be identical for samples to be treated as replicates.

PLGEM-STN statistic determines the degree of differential expression between a condition and the baseline:

PLGEM-STN = [mean(condition)-mean(baseline)] / [modeledSpread(condition)+modeledSpread(baseline)],

where: ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept

\sQuote{plgem.obsStn} determines the observed PLGEM-STN values for each probeset of \sQuote{data}.
\sQuote{plgem.resampledStn} determines the resampled PLGEM-STN values for each probeset of \sQuote{data} using
a resampling approach; see References for details. The number of iterations should be correlated with
the number of replicates of the condition used for model fitting.
}
\value{
\sQuote{plgem.obsStn} returns a matrix of observed PLGEM-STN. The rownames of the matrix are the \sQuote{data} probeset ids.

\sQuote{plgem.resampledStn} returns a list with two items:
  \item{RESAMPLED.STN}{matrix of resampled PLGEM-STN whose rownames of the matrix are the \sQuote{data} probeset ids,
and the columns represent the different replication schema of the provided comparisons; see References for details.}
  \item{REPL.NUMBER}{the number of replicates found for each experimental condition; see References for details.}
}
\references{N. Pavelka et al., BMC Bioinformatics, 2004 Dec 17;5(1):203; \url{http://www.genopolis.it} }
\author{Mattia Pelizzola \email{mattia.pelizzola@unimib.it} and Norman Pavelka \email{norman.pavelka@unimib.it}}
\seealso{\code{\link{plgem.fit}},\code{\link{plgem.deg}},\code{\link{run.plgem}} }
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
LPSobsStn<-plgem.obsStn(data = LPSeset, plgemFit = LPSfit)
LPSresampledStn<-plgem.resampledStn(data = LPSeset, plgemFit = LPSfit)
}
\keyword{models}
