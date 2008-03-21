\name{plgem.obsStn}
\alias{plgem.obsStn}
\alias{plgem.resampledStn}
\title{Computation of Observed and Resampled PLGEM-STN statistics}
\description{
  These functions compute the observed and resampled signal to noise ratio (STN)
  statistic using \bold{PLGEM} fitting parameters (obtained via a call to
  function \sQuote{plgem.fit}) to detect differential expression in an
  ExpressionSet \sQuote{data}, containing either microarray or proteomics data.
}
\usage{
  plgem.obsStn(data, plgemFit, covariateNumb=1, baseline.condition=1,
    verbose=FALSE)
  plgem.resampledStn(data, plgemFit, covariateNumb=1, baseline.condition=1,
    iterations="automatic", verbose=FALSE)
}
\arguments{
  \item{data}{an object of class \sQuote{ExpressionSet}}
  \item{plgemFit}{\code{list}; the output of \sQuote{plgem.fit}.}
  \item{covariateNumb}{\code{integer} (or coercible to \code{integer}); the
    covariate used to determine on which samples to fit the \bold{PLGEM}.}
  \item{baseline.condition}{\code{integer} (or coercible to \code{integer}); the
    condition to be treated as the baseline.}
  \item{verbose}{\code{logical}; if TRUE, comments are printed out while running.}
  \item{iterations}{number of iterations for the resampling step; if
  "automatic" it is automatically determined.}
}
\details{
  The \sQuote{covariateNumb} covariate (the 1st one by default) in the
  \sQuote{pData} of the ExpressionSet \sQuote{data} is expected to contain the
  necessary information about the experimental design. The values of this
  covariate must be sample labels, that have to be identical for samples to be
  treated as replicates. In particular, the ExpressionSet \sQuote{data} must
  have at least two conditions in the \sQuote{covariateNumb} covariate;
  by default the first one is considered the baseline.

  PLGEM-STN statistic determines the degree of differential expression between a
  condition and the baseline:

  PLGEM-STN = [mean(condition)-mean(baseline)] / [modeledSpread(condition)+modeledSpread(baseline)],

  where: ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept

  \sQuote{plgem.obsStn} determines the observed PLGEM-STN values for each gene
  or protein in \sQuote{data}. \sQuote{plgem.resampledStn} determines the
  resampled PLGEM-STN values for each gene or protein in \sQuote{data} using
  a resampling approach; see References for details. The number of iterations
  should be chosen depending on the number of replicates of the condition used
  for fitting the model.
}
\value{
  \sQuote{plgem.obsStn} returns a \code{matrix} of observed PLGEM STN values.
  The \code{\link{rownames}} of this matrix are identical to the
  \code{\link{rownames}} of \sQuote{data}. The \code{\link{colnames}} represent
  the different experimental conditions that were compared to the baseline.

  \sQuote{plgem.resampledStn} returns a list with two items:
    \item{RESAMPLED.STN}{\code{matrix} of resampled PLGEM STN values, with
    \code{\link{rownames}} identical to those in \sQuote{data}, and
    \code{\link{colnames}} representing the different number of replicates
    found in the different comparisons; see References for details.}
    \item{REPL.NUMBER}{the number of replicates found for each experimental
    condition; see References for details.}
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
  \code{\link{plgem.fit}}, \code{\link{plgem.pValue}}, \code{\link{plgem.deg}},
  \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
  LPSobsStn<-plgem.obsStn(data = LPSeset, plgemFit = LPSfit)
  set.seed(123)
  LPSresampledStn<-plgem.resampledStn(data = LPSeset, plgemFit = LPSfit)
}
\keyword{models}
