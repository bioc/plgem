\name{run.plgem}
\alias{run.plgem}
\title{
  Wrapper for Power Law Global Error Model (PLGEM) analysis method
}
\description{
  This function automatically performs \bold{PLGEM} fitting and evaluation,
  determination of observed and resampled \bold{PLGEM} STN values, and
  selection of differentially expressed genes/proteins (DEG) using the
  \bold{PLGEM} method.
}
\usage{
  run.plgem(esdata, signLev=0.001, rank=100, covariateNumb=1,
    baselineCondition=1, Iterations="automatic", fitting.eval=TRUE,
    plotFile=FALSE, writeFiles=FALSE, Verbose=FALSE)
} 
\arguments{
  \item{esdata}{an object of class \code{ExpressionSet}; see Details for
    important information on how the \code{phenoData} slot of this object will
    be interpreted by the function.}
  \item{signLev}{numeric vector; significance level(s) for the DEG selection.
    Value(s) must be in (0,1).}
  \item{rank}{\code{integer} (or coercible to \code{integer}); the number of
    genes or proteins to be selected according to their PLGEM-STN rank. Only
    used if number of available replicates is too small to perform resampling
    (see Details).}
  \item{covariateNumb}{\code{integer} (or coercible to \code{integer}); the
    covariate used to determine on which samples to fit plgem.}
  \item{baselineCondition}{\code{integer} (or coercible to \code{integer}); the
    condition to be treated as the baseline.}
  \item{Iterations}{number of iterations for the resampling step; if "automatic"
    it is automatically determined.}
  \item{fitting.eval}{\code{logical}; if \code{TRUE}, the fitting is evaluated
    generating a diagnostic plot.}
  \item{plotFile}{\code{logical}; if \code{TRUE}, the generated plot is written on a
    file.}
  \item{writeFiles}{\code{logical}; if \code{TRUE}, the generated list of DEG is
    written on disk file(s).}
  \item{Verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  The \sQuote{covariateNumb} covariate (the first one by default) of the
  \code{phenoData} of the \code{ExpressionSet} \sQuote{data} is expected to
  contain the necessary information about the experimental design. The values of
  this covariate must be sample labels, that have to be identical for samples to
  be treated as replicates. In particular, the \code{ExpressionSet}
  \sQuote{esdata} must have at least two conditions in the \sQuote{covariateNumb}
  covariate; by default the first one is considered the baseline.

  The model is fitted on the most replicated condition. When more conditions
  exist with the max number of replicates, the condition providing the best fit
  is chosen.

  If less than 3 replicates are provided for the condition used for fitting,
  then the selection is based on ranking according to the observed \bold{PLGEM}
  STN values. In this case the first \sQuote{rank} genes or proteins are
  selected for each comparison.

  Otherwise DEG are selected comparing the observed and resampled \bold{PLGEM}
  STN values at the \sQuote{signLev} significance level(s), based on p-values
  obtained via a call to function \code{\link{plgem.pValue}}. See References for
  details.
}
\value{
  This function returns a list with a number of items that is equal to the
  number of different significance levels (\sQuote{signLev}) used as input. Each
  item is again a list, whose number of items correspond to the number of
  performed comparisons, i.e. the number of conditions defined in the
  \code{phenoData} of \sQuote{esdata} minus the baseline. In each list-item the
  values are the observed \bold{PLGEM} STN values of the significantly changing
  genes or proteins, named according to the \code{rownames} of the \code{exprs}
  of \sQuote{esdata}.
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
  \code{\link{plgem.resampledStn}}, \code{\link{plgem.pValue}},
  \code{\link{plgem.write.summary}}
}
\examples{
  data(LPSeset)
  set.seed(123)
  LPSdegList <- run.plgem(esdata=LPSeset)
}
\keyword{models}
