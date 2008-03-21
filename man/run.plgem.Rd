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
  \item{esdata}{an object of class \code{ExpressionSet}.}
  \item{signLev}{number or array; significance level(s) for the DEG selection, value(s) must be in (0,1).}
  \item{rank}{number; the number of genes to be selected according to their PLGEM-STN rank; see details.}
  \item{covariateNumb}{\code{integer} (or coercible to \code{integer}); the covariate used to determine on which samples to fit plgem.}
  \item{baselineCondition}{\code{integer} (or coercible to \code{integer}); the condition to be treated as the baseline.}
  \item{Iterations}{number of iterations for the resampling step; if "automatic" it is automatically determined.}
  \item{fitting.eval}{\code{logical}; if TRUE, the fitting is evaluated generating a diagnostic plot.}
  \item{plotFile}{\code{logical}; if TRUE, the generated plot is written on a file.}
  \item{writeFiles}{\code{logical}; if TRUE, the generated list of DEG is written on disk file(s).}
  \item{Verbose}{\code{logical}; if TRUE, comments are printed out while running.}
}
\details{
The covariateNumb covariate (the 1st one by default) of the phenoData of the ExpressionSet \sQuote{data} is
expected to contain the necessary information about te experimental design. The values of this covariate
must be sample labels, that have to be identical for samples to be treated as replicates.
In particular, the ExpressionSet \sQuote{data} must have at least two conditions in the \sQuote{covariateNumb} covariate;
by default the first one is considered the baseline.

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
performed comparisons, i.e. the number of conditions in the starting ExpressionSet minus the baseline. In each list-item the values are
the observed PLGEM-STN and the names are the DEG probeset ids.
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
  \code{\link{plgem.write.summary}}
}
\examples{
  data(LPSeset)
  set.seed(123)
  LPSdegList <- run.plgem(esdata=LPSeset)
}
\keyword{models}
