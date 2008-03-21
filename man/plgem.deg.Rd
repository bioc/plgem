\name{plgem.deg}
\alias{plgem.deg}
\title{
  Selection of differentially expressed genes/proteins using PLGEM
}
\description{
  This function selects differentially expressed genes/proteins (DEG) at a given
  significance level \sQuote{delta}, based on observed \bold{PLGEM}
  signal-to-noise ratio (STN) values (typically obtained via a call to
  \code{\link{plgem.obsStn}}) and pre-computed p-values (typically obtained via
  a call to \code{\link{plgem.pValue}}).
}
\usage{
  plgem.deg(observedStn, plgemPval, delta=0.001, verbose=FALSE)
}
\arguments{
  \item{observedStn}{\code{matrix} of observed STN values; output of function
    \code{\link{plgem.obsStn}}.}
  \item{plgemPval}{\code{matrix} of p-values; output of function
    \code{\link{plgem.pValue}}.}
  \item{delta}{numeric vector; the significance level(s) to be used for the
    selection of DEG; value(s) must be between 0 and 1 (excluded).}
  \item{verbose}{\code{logical}; if TRUE, comments are printed out while running.}
}
\details{
  This function allows for the selection of DEG by setting a significance
  cut-off on pre-calculated p-values. The significance level \sQuote{delta}
  roughly represents the false positive rate of the DEG selection, e.g. if a
  \sQuote{delta} of 0.001 is chosen in a microarray dataset with 10000 genes,
  on average 10 of the selected DEG are expected to be false positives.
}
\value{
  This function returns a list with a number of items equal to the number of
  different significance levels (\sQuote{delta}) used as input. Each item of
  this list is again a list, whose number of items correspond to the number of
  performed comparisons (i.e. the number of conditions in the starting
  ExpressionSet minus the baseline). Each of these second level list-items is a
  vector of observed STN values of the genes or proteins that passed the
  corresponding significance threshold in the corresponding comparison.
}
\references{
  Pavelka N, Pelizzola M, Vizzardelli C, Capozzoli M, Splendiani A, Granucci F,
  Ricciardi-Castagnoli P. A power law global error model for the identification
  of differentially expressed genes in microarray data. BMC Bioinformatics. 2004
  Dec 17;5:203.; \url{http://www.biomedcentral.com/1471-2105/5/203}

  Pavelka N, Fournier ML, Swanson SK, Pelizzola M, Ricciardi-Castagnoli P,
  Florens L, Washburn MP. Statistical similarities between transcriptomics and
  quantitative shotgun proteomics data. Mol Cell Proteomics. 2007 Nov 19;
  \url{http://www.mcponline.org/cgi/content/abstract/M700240-MCP200v1}
}
\author{
  Mattia Pelizzola \email{mattia.pelizzola@gmail.com}
  Norman Pavelka \email{nxp@stowers-institute.org}
}
\seealso{
  \code{\link{plgem.fit}}, \code{\link{plgem.obsStn}},
  \code{\link{plgem.resampledStn}}, \code{\link{plgem.pValue}},
  \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset, fittingEval=TRUE)
  LPSobsStn <- plgem.obsStn(data=LPSeset, plgemFit=LPSfit)
  set.seed(123)
  LPSresampledStn <- plgem.resampledStn(data=LPSeset, plgemFit=LPSfit)
  LPSpValues <- plgem.pValue(LPSobsStn, LPSresampledStn)
  LPSdegList <- plgem.deg(observedStn=LPSobsStn, plgemPval=LPSpValues, delta=0.001)
}
\keyword{models}
