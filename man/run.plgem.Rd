\name{run.plgem}
\alias{run.plgem}
\title{Wrapper for Power Law Global Error Model PLGEM}
\description{
This function automatically performs model fitting and evaluation, observed and resampled PLGEM-STN determination and
selection of differentially expressed genes using the PLGEM method.
}
\usage{
run.plgem(esdata, signLev = 0.001, rank = 100, baselineCondition = 1, Iterations = "automatic",
fitting.eval = TRUE, plotFile = FALSE, Verbose = FALSE)
} 
\arguments{
  \item{esdata}{an object of class \sQuote{exprSet} with a \sQuote{conditionName} covariate; see details.}
  \item{signLev}{number or array; significance level(s) for the DEG selection, value(s) must be in (0,1).}
  \item{rank}{number; the number of genes to be selected according to their PLGEM-STN rank; see details.}
  \item{baselineCondition}{number; the condition to be treated as the baseline.}
  \item{Iterations}{number of iterations for the resampling step; if "automatic" it is automatically determined.}
  \item{fitting.eval}{logical; if TRUE, the fitting is evaluated generating a diagnostic plot.}
  \item{plotFile}{logical; if TRUE, the generated plot is written on a file.}
  \item{Verbose}{logical; if TRUE, comments are printed out while running.}
}
\details{
The exprSet \sQuote{data} must have at least two conditions; by default the first one is considered the baseline.
\sQuote{data} must have a phenoData slot with a covariate called \sQuote{conditionName}. The values of this covariate
must be sample labels, that have to be identical for samples to be treated as replicates.

The model is fitted on the most replicated condition. When more conditions exist with the max number of
replicates, the condition providing the best fit is chosen.

If less than 3 replicates are provided for the condition used for fitting, then the selection is based on ranking according to the
observed PLGEM-STN statistics. In this case the first \sQuote{rank} genes are selected for each comparison.

Otherwise DEG are selected comparing the observed and resampled PLGEM-STN at the \sQuote{sign.lev} significance level(s), that
can be treated as an estimate of the false positive rate. See References for details. 
}
\value{
This function returns a list with a number of items that is equal to the number of different significance levels
(\sQuote{delta}) used as input. Each item is again a list, whose number of items correspond to the number of
performed comparisons, i.e. the number of conditions in the starting exprSet minus the baseline. In each list-item the values are
the observed PLGEM-STN and the names are the DEG probeset ids.
}
\references{N. Pavelka et al., BMC Bioinformatics, 2004 Dec 17;5(1):203; \url{http://www.genopolis.it} }
\author{Mattia Pelizzola \email{mattia.pelizzola@unimib.it} and Norman Pavelka \email{norman.pavelka@unimib.it}}
\seealso{\code{\link{plgem.fit}},\code{\link{plgem.obsStn}},\code{\link{plgem.resampledStn}} }
\examples{
data(LPSeset)
LPSdegList <- run.plgem(esdata = LPSeset)
}
\keyword{models}
