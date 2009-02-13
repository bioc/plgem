\name{plgem.fit}
\alias{plgem.fit}
\title{PLGEM Fitting and Evaluation}
\description{
  Function for fitting and evaluating goodness of fit of a \bold{PLGEM} on a set
  of replicated samples defined in an \code{ExpressionSet}.
}
\usage{
  plgem.fit(data, covariateNumb=1, fit.condition=1, p=10, q=0.5,
    trimAllZeroRows=FALSE, zeroMeanOrSD=c("replace", "trim"), fittingEval=FALSE, 
    plot.file=FALSE, verbose=FALSE)
}
\arguments{
  \item{data}{an object of class \code{ExpressionSet}; see Details for important
    information on how the \code{phenoData} slot of this object will be
    interpreted by the function.}
  \item{covariateNumb}{\code{integer} (or coercible to \code{integer}); the
    covariate used to determine on which samples to fit the \bold{PLGEM}.}
  \item{fit.condition}{\code{integer} (or coercible to \code{integer}); the
    condition used for \bold{PLGEM} fitting, according to the order of unique
    values of the \sQuote{covariateNumb} covariate.}
  \item{p}{\code{integer} (or coercible to \code{integer}); number of intervals
    used to partition the expression value range.}
  \item{q}{\code{numeric} in [0,1]; the quantile of standard deviation used for
  \bold{PLGEM} fitting.}
  \item{trimAllZeroRows}{\code{logical}; if \code{TRUE}, rows in the data set
    containing only zero values are trimmed before fitting \bold{PLGEM}.}
  \item{zeroMeanOrSD}{either \code{NULL} or \code{character}; what should be
    done if a row with non-positive mean or zero standard deviation is
    encountered before fitting \bold{PLGEM}? Current options are one of
    \code{"replace"} or \code{"trim"}. Partial matching is used to switch
    between the options and setting the value to \code{NULL} will cause the
    default behaviour to be enforced, i.e. to \code{"replace"} (see Details).}
  \item{fittingEval}{\code{logical}; if \code{TRUE}, the fitting is evaluated
    generating a diagnostic plot.}
  \item{plot.file}{\code{logical}; if \code{TRUE}, a png file is written on the
    current working directory.}
  \item{verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  \code{plgem.fit} fits a Power Law Global Error Model (\bold{PLGEM}) to an
  expression set and optionally evaluates the quality of the fit. This
  \bold{PLGEM} aims to find the mathematical relationship between standard
  deviation and mean gene expression values (or protein abundance levels) in a
  set of replicated microarray (or proteomics) samples, according to the
  following power law:

  ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept

  It has been demonstrated that this model fits to Affymetrix GeneChip datasets,
  as well as to datasets of normalized spectral counts obtained by mass
  spectrometry-based proteomics (see References for details).
  
  The \code{covariateNumb} covariate (the first one by default) defined in the
  input \code{ExpressionSet} is expected to contain the necessary information
  to distinguish the various experimental conditions from one another. The
  values of this covariate must be condition names, which have to be identical
  for samples that are to be treated as replicates.
  
  Setting \code{trimAllZeroRows = TRUE} is especially useful in proteomics data
  sets, where there is no guarantee of identifying a protein across all
  experimental conditions. Since \bold{PLGEM} is fitted only to the data
  corresponding to a single experimental condition (as defined by
  \code{fit.condition}), it is possible to generate a non-negligible number of
  rows containing only zero values, even if there were no such rows in the
  original (complete) data set containing all experimental conditions.

  Setting \code{zeroMeanOrSD = "replace"} (the current default, for backward
  compatibility) will cause the function to replace zero or negative means with
  the smallest positive mean found in the data set and to replace zero standard
  deviations with the smallest non-zero standard deviation found in the data
  set. Setting \code{zeroMeanOrSD = "trim"} is the current recommended option,
  especially for spectral counting proteomics data sets that are typically
  characterized by a high data granularity or for microarray data sets with a
  small number of replicates. In both cases, there are chances for data values
  for a same gene or protein to be identical across replicates (and therefore
  with zero standard deviation) by chance alone. Note that setting
  \code{trimAllZeroRows = TRUE} does not guarantee that there will be no rows
  with zero mean or zero standard deviation.

  If argument \sQuote{fittingEval} is \code{TRUE}, a graphical control of the
  goodness of the \bold{PLGEM} fitting is produced and a plot containing four
  panels is generated. The top-left panel shows the power law, characterized by
  \sQuote{SLOPE} and \sQuote{INTERCEPT}. The top-right panel represents the
  distribution of model residuals. The bottom-left reports the contour plot of
  ranked residuals. The bottom-right panel finally shows the relationship
  between the distribution of observed residuals and the normal distribution.
  A good fit normally gives a horizontal symmetric rank-plot and a near normal
  distribution of residuals.
}
\value{
\code{plgem.fit} returns a list of five numbers (see Details):
  \item{SLOPE}{the slope of the fitted PLGEM.}
  \item{INTERCEPT}{the intercept of the fitted PLGEM.}
  \item{DATA.PEARSON}{the Pearson correlation coefficient of the linear model
  fitted on the original data.}
  \item{ADJ.R2.MP}{the adjusted R squared of PLGEM fitted on the modelling
  points.}
  \item{FIT.CONDITION}{the condition used for fitting PLGEM.}
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
  \code{\link{plgem.obsStn}}, \code{\link{plgem.resampledStn}},
  \code{\link{plgem.pValue}}, \code{\link{plgem.deg}}, \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset, fittingEval=FALSE)
  sapply(LPSfit, function(x) return(as.vector(x)))
}
\keyword{models}
