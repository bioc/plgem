\name{plgem.pValue}
\alias{plgem.pValue}
\title{Computation of PLGEM p-values}
\description{
  This function computes p-values for observed PLGEM signal-to-noise ratio (STN)
  values (typically obtained via a call to \code{\link{plgem.obsStn}}) from
  resampled STN values (typically obtained via a call to
  \code{\link{plgem.resampledStn}}).
}
\usage{
  plgem.pValue(observedStn, plgemResampledStn, verbose=FALSE)
}
\arguments{
  \item{observedStn}{\code{matrix} of observed PLGEM STN values; output of
    function \code{\link{plgem.obsStn}}.}
  \item{plgemResampledStn}{\code{list}; output of function
    \code{\link{plgem.resampledStn}}.}
  \item{verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  The p-value of each given observed STN value is computed based on the quantile
  that the given value occupies in the corresponding distribution of
  resampled PLGEM STN values, based on the following relationship:
  
  p-value = min(2*quantile, 2*(1-quantile))
}
\value{
  \code{plgem.pValue} returns a matrix with the same \code{\link{dim}}ensions
  and \code{\link{dimnames}} as the input \code{observedStn}, where each entry
  represents the p-value of the corresponding observed PLGEM STN value.
}
\references{
  Pavelka N, Pelizzola M, Vizzardelli C, Capozzoli M, Splendiani A, Granucci F,
  Ricciardi-Castagnoli P. A power law global error model for the identification
  of differentially expressed genes in microarray data. BMC Bioinformatics. 2004
  Dec 17; 5:203; \url{http://www.biomedcentral.com/1471-2105/5/203}.

  Pavelka N, Fournier ML, Swanson SK, Pelizzola M, Ricciardi-Castagnoli P,
  Florens L, Washburn MP. Statistical similarities between transcriptomics and
  quantitative shotgun proteomics data. Mol Cell Proteomics. 2008 Apr;
  7(4):631-44; \url{http://www.mcponline.org/cgi/content/full/7/4/631}.
}
\author{
  Mattia Pelizzola \email{mattia.pelizzola@gmail.com}
  Norman Pavelka \email{nxp@stowers.org}
}
\seealso{
  \code{\link{plgem.fit}}, \code{\link{plgem.obsStn}},
  \code{\link{plgem.resampledStn}}, \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset)
  LPSobsStn <- plgem.obsStn(data=LPSeset, plgemFit=LPSfit)
  head(LPSobsStn)
  set.seed(123)
  LPSresampledStn <- plgem.resampledStn(data=LPSeset, plgemFit=LPSfit)
  LPSpValues <- plgem.pValue(LPSobsStn, LPSresampledStn)
  head(LPSpValues)
}
\keyword{models}
