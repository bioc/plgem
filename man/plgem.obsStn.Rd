\name{plgem.obsStn}
\alias{plgem.obsStn}
\title{Computation of Observed PLGEM-STN Statistics}
\description{
  This function computes observed signal-to-noise ratio (STN) values using
  \bold{PLGEM} fitting parameters (obtained via a call to function
  \code{\link{plgem.fit}}) to detect differential expression in an
  \code{ExpressionSet}, containing either microarray or proteomics data.
}
\usage{
  plgem.obsStn(data, plgemFit, covariate=1, baselineCondition=1,
    verbose=FALSE)
}
\arguments{
  \item{data}{an object of class \code{ExpressionSet}; see Details for important
    information on how the \code{phenoData} slot of this object will be
    interpreted by the function.}
  \item{plgemFit}{\code{list}; the output of function \code{\link{plgem.fit}}.}
  \item{covariate}{\code{integer}, \code{numeric} or \code{character}; specifies
    the covariate to be used to distinguish the various experimental conditions
    from one another. See Details for how to specify the \code{covariate}.}
  \item{baselineCondition}{\code{integer}, \code{numeric} or \code{character};
    specifies the condition to be treated as the baseline. See Details for how
    to specify the \code{baselineCondition}.}
  \item{verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  The \code{phenoData} slot of the \code{ExpressionSet} given as input is
  expected to contain the necessary information to distinguish the various
  experimental conditions from one another. The columns of the \code{pData} are
  referred to as \sQuote{covariates}. There has to be at least one covariate
  defined in the input \code{ExpressionSet}. The sample attributes according to
  this covariate must be distinct for samples that are to be treated as distinct
  experimental conditions and identical for samples that are to be treated as
  replicates.
  
  There is a couple different ways how to specify the \code{covariate}: If an
  \code{integer} or a \code{numeric} is given, it will be taken as the covariate
  number (in the same order in which the covariates appear in the
  \code{colnames} of the \code{pData}). If a \code{character} is given, it will
  be taken as the covariate name itself (in the same way the covariates are
  specified in the \code{colnames} of the \code{pData}). By default, the first
  covariate appearing in the \code{colnames} of the \code{pData} is used.
  
  Similarly, there is a couple different ways how to specify which experimental
  condition to treat as the baseline. The available \sQuote{condition names} are
  taken from \code{unique(as.character(pData(data)[, covariate]))}. If
  \code{baselineCondition} is given as a \code{character}, it will be taken as
  the condition name itself. If \code{baselineCondition} is given as an
  \code{integer} or a \code{numeric} value, it will be taken as the condition
  number (in the same order of appearance as in the \sQuote{condition names}).
  By default, the first condition name is used.
  
  PLGEM-STN values are a measure of the degree of differential expression
  between a condition and the baseline: \deqn{
  STN = \frac{mean_{condition}-mean_{baseline}}{modeledSpread_{condition}+modeledSpread_{baseline}},}{%
  STN = [mean(condition)-mean(baseline)] / [modeledSpread(condition)+modeledSpread(baseline)],}
  where: \deqn{\log{(modeledSpread)}=PLGEMslope*\log{(mean)}+PLGEMintercept}{%
  ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept}

  \code{plgem.obsStn} determines the observed PLGEM-STN values for each gene
  or protein in \code{data}; see References for details.
}
\value{
  A \code{list} of two elements:
  \item{fit}{the input \code{plgemFit}.}
  \item{PLGEM.STN}{a \code{matrix} of observed PLGEM-STN values. The
    \code{\link{rownames}} of this matrix are identical to the \code{rownames}
    of \code{data}. The \code{\link{colnames}} represent the different
    experimental conditions that were compared to the baseline.}
}
\references{
  Pavelka N, Pelizzola M, Vizzardelli C, Capozzoli M, Splendiani A, Granucci F,
  Ricciardi-Castagnoli P. A power law global error model for the identification
  of differentially expressed genes in microarray data. BMC Bioinformatics. 2004
  Dec 17; 5:203; \url{http://www.biomedcentral.com/1471-2105/5/203}.

  Pavelka N, Fournier ML, Swanson SK, Pelizzola M, Ricciardi-Castagnoli P,
  Florens L, Washburn MP. Statistical similarities between transcriptomics and
  quantitative shotgun proteomics data. Mol Cell Proteomics. 2008 Apr;
  7(4):631-44; \url{http://www.mcponline.org/cgi/content/abstract/7/4/631}.
}
\author{
  Mattia Pelizzola \email{mattia.pelizzola@gmail.com}

  Norman Pavelka \email{normanpavelka@gmail.com}
}
\seealso{
  \code{\link{plgem.fit}}, \code{\link{plgem.resampledStn}},
  \code{\link{plgem.pValue}}, \code{\link{plgem.deg}}, \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset)
  LPSobsStn <- plgem.obsStn(data=LPSeset, plgemFit=LPSfit)
  head(LPSobsStn[["PLGEM.STN"]])
}
\keyword{models}
