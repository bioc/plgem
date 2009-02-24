\name{plgem.obsStn}
\alias{plgem.obsStn}
\alias{plgem.resampledStn}
\title{Computation of Observed and Resampled PLGEM-STN statistics}
\description{
  These functions compute observed and resampled signal-to-noise ratio (STN)
  values using \bold{PLGEM} fitting parameters (obtained via a call to
  function \code{\link{plgem.fit}}) to detect differential expression in an
  \code{ExpressionSet} \sQuote{data}, containing either microarray or proteomics
  data.
}
\usage{
  plgem.obsStn(data, plgemFit, covariate=1, baselineCondition=1,
    verbose=FALSE)
  plgem.resampledStn(data, plgemFit, covariate=1, baselineCondition=1,
    iterations="automatic", verbose=FALSE)
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
  \item{iterations}{number of iterations for the resampling step; if
    \code{"automatic"} it is automatically determined.}
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
  be taken as the covariate name (in the same way the covariates are specified
  in the \code{colnames} of the \code{pData}). By default, the first covariate
  appearing in the \code{colnames} of the \code{pData} is used.
  
  Similarly, there is a couple different ways how to specify which experimental
  condition to treat as the baseline. The available \sQuote{condition names} are
  taken from \code{unique(as.character(pData(data)[, covariate]))}. If
  \code{baselineCondition} is given as a \code{character}, it will be taken as
  the condition name itself. If \code{baselineCondition} is given as an
  \code{integer} or a \code{numeric} value, it will be taken as the condition
  number (in the same order of appearance in the \sQuote{condition names}). By
  default, the first appearing condition name is used.
  
  PLGEM-STN values are a measure of the degree of differential expression
  between a condition and the baseline:

  PLGEM-STN = [mean(condition)-mean(baseline)] / [modeledSpread(condition)+modeledSpread(baseline)],

  where: ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept

  \code{plgem.obsStn} determines the observed PLGEM-STN values for each gene
  or protein in \code{data}. \code{plgem.resampledStn} determines the resampled
  PLGEM-STN values for each gene or protein in \code{data} using a resampling
  approach; see References for details. The number of iterations should be
  chosen depending on the number of available replicates of the condition used
  for fitting the model.
}
\value{
  \code{plgem.obsStn} returns a \code{matrix} of observed PLGEM-STN
  values. The \code{\link{rownames}} of this matrix are identical to the
  \code{rownames} of \code{data}. The \code{\link{colnames}} represent
  the different experimental conditions that were compared to the baseline.

  \code{plgem.resampledStn} returns a list with two items:
    \item{RESAMPLED.STN}{\code{matrix} of resampled PLGEM-STN values, with
    \code{\link{rownames}} identical to those in \code{data}, and
    \code{\link{colnames}} representing the different number of replicates
    found in the different comparisons; see References for details.}
    \item{REPL.NUMBER}{the number of replicates found for each experimental
    condition; see References for details.}
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
  \code{\link{plgem.fit}}, \code{\link{plgem.pValue}}, \code{\link{plgem.deg}},
  \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset)
  LPSobsStn <- plgem.obsStn(data=LPSeset, plgemFit=LPSfit)
  set.seed(123)
  LPSresampledStn <- plgem.resampledStn(data=LPSeset, plgemFit=LPSfit)
  plot(density(LPSresampledStn[["RESAMPLED.STN"]], bw=0.01), col="black", lwd=2,
    xlab="PLGEM STN values",
    main="Distribution of observed and resampled PLGEM STN values")
  lines(density(LPSobsStn, bw=0.01), col="red")
  legend("topright", legend=c("resampled", "observed"), col=c("black", "red"),
    lwd=2:1)
}
\keyword{models}
