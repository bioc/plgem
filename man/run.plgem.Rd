\name{run.plgem}
\alias{run.plgem}
\title{
  Wrapper for Power Law Global Error Model (PLGEM) analysis method
}
\description{
  This function automatically performs \bold{PLGEM} fitting and evaluation,
  determination of observed and resampled PLGEM-STN values, and selection of
  differentially expressed genes/proteins (DEG) using the \bold{PLGEM} method.
}
\usage{
  run.plgem(esdata, signLev=0.001, rank=100, covariate=1,
    baselineCondition=1, Iterations="automatic", trimAllZeroRows=FALSE,
    zeroMeanOrSD=c("replace", "trim"), fitting.eval=TRUE,
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
  \item{covariate}{\code{integer}, \code{numeric} or \code{character}; specifies
    the covariate to be used to distinguish the various experimental conditions
    from one another. See Details for how to specify the \code{covariate}.}
  \item{baselineCondition}{\code{integer}, \code{numeric} or \code{character};
    specifies the condition to be treated as the baseline. See Details for how
    to specify the \code{baselineCondition}.}
  \item{Iterations}{number of iterations for the resampling step; if
    \code{"automatic"} it is automatically determined.}
  \item{trimAllZeroRows}{\code{logical}; if \code{TRUE}, rows in the data set
    containing only zero values are trimmed before fitting \bold{PLGEM}. See
    help page of function \code{\link{plgem.fit}} for details.}
  \item{zeroMeanOrSD}{either \code{NULL} or \code{character}; what should be
    done if a row with non-positive mean or zero standard deviation is
    encountered before fitting \bold{PLGEM}? Current options are one of
    \code{"replace"} or \code{"trim"}. Partial matching is used to switch
    between the options and setting the value to \code{NULL} will cause the
    default behaviour to be enforced, i.e. to \code{"replace"}. See help page of
    function \code{\link{plgem.fit}} for details.}
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
  
  The model is fitted on the most replicated condition. When more conditions
  exist with the maximum number of replicates, the condition providing the best
  fit is chosen (based on the adjusted \eqn{r^{2}}{r^2}). If there is again a
  tie, the first one is arbitrarily taken.

  If less than 3 replicates are provided for the condition used for fitting,
  then the selection is based on ranking according to the observed PLGEM-STN
  values. In this case the first \code{rank} genes or proteins are selected for
  each comparison.

  Otherwise DEG are selected comparing the observed and resampled PLGEM-STN
  values at the \code{signLev} significance level(s), based on p-values obtained
  via a call to function \code{\link{plgem.pValue}}. See References for details.
}
\value{
  This function returns a list with a number of items that is equal to the
  number of different significance levels (\code{signLev}) used as input. Each
  item is again a list, whose number of items correspond to the number of
  performed comparisons, i.e. the number of conditions defined in the
  \code{phenoData} of the input \code{ExpressionSet} minus the baseline. In each
  list-item the values are the observed PLGEM-STN values of the significantly
  changing genes or proteins, named according to the \code{featureNames} of the
  input \code{ExpressionSet}.
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
  LPSdegList <- run.plgem(esdata=LPSeset, fitting.eval=FALSE)
}
\keyword{models}
