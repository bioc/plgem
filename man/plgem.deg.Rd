\name{plgem.deg}
\alias{plgem.deg}
\title{Selection of Differentially Expressed Genes using PLGEM}
\description{
This function compares observed PLGEM-derived STN statistics with resampled test statistics, and selects
the differentially expressed genes (DEG) with a given significance level \sQuote{delta}.  
}
\usage{plgem.deg(observedStn, plgemResampledStn, delta = 0.001, verbose = FALSE)}
\arguments{
  \item{observedStn}{matrix of observed signal to noise ratios (STN); output of function \sQuote{plgem.obsStn}.}
  \item{plgemResampledStn}{list; the output of function \sQuote{plgem.resampledStn}.}
  \item{delta}{number or array; the significance level(s) used for the DEG selection, value(s) must be in (0,1).}
  \item{verbose}{logical; if TRUE, comments are printed out while running.}
}
\details{
This function allows for the selection of DEG by comparing their observed and resampled STN.
\sQuote{delta} roughly represents the false positive rate of the DEG selection, e.g. if a \sQuote{delta} of 0.001
is chosen in a microarray with 10000 entries, 10 of the selected DEG are expected to be false positives.
}
\value{
This function returns a list with a number of items is equal to the number of different significance levels
(\sQuote{delta}) used as input. Each item is again a list, whose number of items correspond to the number of
performed comparisons, i.e. the number of conditions in the starting exprSet minus the baseline. In each list-item the values are
the observed STN and the names are the DEG probeset ids.
}
\references{N. Pavelka et al., BMC Bioinformatics, 2004 Dec 17;5(1):203; \url{http://www.genopolis.it} }
\author{Mattia Pelizzola \email{mattia.pelizzola@unimib.it} and Norman Pavelka \email{norman.pavelka@unimib.it}}
\seealso{\code{\link{plgem.fit}},\code{\link{plgem.obsStn}},\code{\link{plgem.resampledStn}},\code{\link{run.plgem}} }
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
LPSobsStn<-plgem.obsStn(data = LPSeset, plgemFit = LPSfit)
LPSresampledStn<-plgem.resampledStn(data = LPSeset, plgemFit = LPSfit)
LPSdegList<-plgem.deg(observedStn = LPSobsStn, plgemResampledStn = LPSresampledStn, delta = 0.001)
}
\keyword{models}
