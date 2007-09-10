\name{plgem.obsStn}
\alias{plgem.obsStn}
\alias{plgem.resampledStn}
\title{Computation of Observed and Resampled PLGEM-STN statistics}
\description{
These functions compute the observed and resampled signal to noise ratio (STN) statistic
using PLGEM fitting parameters (obtained using the function \sQuote{plgem.fit}) to detect differential
expression in a microarray ExpressionSet \sQuote{data}.
}
\usage{plgem.obsStn(data, plgemFit, covariateNumb = 1, baseline.condition = 1, verbose = FALSE)
plgem.resampledStn(data, plgemFit,  covariateNumb = 1, baseline.condition = 1,
iterations = "automatic", verbose=FALSE)
}
\arguments{
  \item{data}{an object of class \sQuote{ExpressionSet}}
  \item{plgemFit}{list; the output of \sQuote{plgem.fit}.}
  \item{covariateNumb}{number; the covariate used to determine on which samples to fit plgem.}
  \item{baseline.condition}{the condition to be treated as the baseline.}
  \item{verbose}{logical; if TRUE, comments are printed out while running.}
  \item{iterations}{number of iterations for the resampling step; if "automatic" it is automatically determined.}
}
\details{

The covariateNumb covariate (the 1st one by default) of the phenoData of the ExpressionSet \sQuote{data} is
expected to contain the necessary information about te experimental design. The values of this covariate
must be sample labels, that have to be identical for samples to be treated as replicates.
In particular, the ExpressionSet \sQuote{data} must have at least two conditions in the \sQuote{covariateNumb} covariate;
by default the first one is considered the baseline.

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
\author{Mattia Pelizzola \email{mattia.pelizzola@gmail.com} and Norman Pavelka \email{NXP@stowers-institute.org}}
\seealso{\code{\link{plgem.fit}},\code{\link{plgem.deg}},\code{\link{run.plgem}} }
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
LPSobsStn<-plgem.obsStn(data = LPSeset, plgemFit = LPSfit)
LPSresampledStn<-plgem.resampledStn(data = LPSeset, plgemFit = LPSfit)
}
\keyword{models}
