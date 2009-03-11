\name{plgem.fit}
\alias{plgem.fit}
\title{PLGEM Fitting and Evaluation}
\description{
  Function for fitting and evaluating goodness of fit of a \bold{PLGEM} on a set
  of replicated samples defined in an \code{ExpressionSet}.
}
\usage{
  plgem.fit(data, covariate=1, fitCondition=1, p=10, q=0.5,
    trimAllZeroRows=FALSE, zeroMeanOrSD=c("replace", "trim"), fittingEval=FALSE, 
    plot.file=FALSE, verbose=FALSE)
}
\arguments{
  \item{data}{an object of class \code{ExpressionSet}; see Details for important
    information on how the \code{phenoData} slot of this object will be
    interpreted by the function.}
  \item{covariate}{\code{integer}, \code{numeric} or \code{character}; specifies
    the covariate to be used to fit the \bold{PLGEM}. See Details for how to
    specify the \code{covariate}.}
  \item{fitCondition}{\code{integer}, \code{numeric} or \code{character};
    specifies the condition to be used to fit the \bold{PLGEM}. See Details for
    how to specify the \code{fitCondition}.}
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
  \code{ExpressionSet} and optionally evaluates the quality of the fit. This
  \bold{PLGEM} aims to find the mathematical relationship between standard
  deviation and mean gene expression values (or protein abundance levels) in a
  set of replicated microarray (or proteomics) samples, according to the
  following power law:
  \deqn{\log{(modeledSpread)}=PLGEMslope*\log{(mean)}+PLGEMintercept}{%
  ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept}
  It has been demonstrated that this model fits to Affymetrix GeneChip datasets,
  as well as to datasets of normalized spectral counts obtained by mass
  spectrometry-based proteomics. Technically, two replicates are required and
  sufficient to fit a \bold{PLGEM}. Having 3 or more replicates, of course,
  improves the fitting and is recommended (see References for details).
  
  The \code{phenoData} slot of the \code{ExpressionSet} given as input is
  expected to contain the necessary information to distinguish the various
  experimental conditions from one another. The columns of the \code{pData} are
  referred to as \sQuote{covariates}. There has to be at least one covariate
  defined in the input \code{ExpressionSet}. The sample attributes according to
  this covariate must be distinct for samples that are to be treated as distinct
  experimental conditions and identical for samples that are to be treated as
  replicates.
  
  There is a couple different ways to specify the \code{covariate}: If an
  \code{integer} or a \code{numeric} is given, it will be taken as the covariate
  number (in the same order in which the covariates appear in the
  \code{colnames} of the \code{pData}). If a \code{character} is given, it will
  be taken as the covariate name itself (in the same way the covariates are
  specified in the \code{colnames} of the \code{pData}). By default, the first
  covariate appearing in the \code{colnames} of the \code{pData} is used.
  
  Similarly, there is a couple different ways to specify on which experimental
  condition to fit the model. The available \sQuote{condition names} are taken
  from \code{unique(as.character(pData(data)[, covariate]))}. If
  \code{fitCondition} is given as a \code{character}, it will be taken as the
  condition name itself. If \code{fitCondition} is given as an \code{integer}
  or a \code{numeric} value, it will be taken as the condition number (in the
  same order of appearance as in the \sQuote{condition names}). By default, the
  first condition name is used.
  
  Setting \code{trimAllZeroRows=TRUE} is especially useful in proteomics data
  sets, where there is no guarantee of identifying a protein across all
  experimental conditions. Since \bold{PLGEM} is fitted only to the data
  corresponding to a single experimental condition (as defined by
  \code{fitCondition}), it is possible to generate a non-negligible number of
  rows containing only zero values, even if there were no such rows in the
  original (complete) data set containing all experimental conditions.

  Setting \code{zeroMeanOrSD="replace"} (the current default, for backward
  compatibility) will cause the function to replace zero or negative means with
  the smallest positive mean found in the data set and to replace zero standard
  deviations with the smallest non-zero standard deviation found in the data
  set. Setting \code{zeroMeanOrSD="trim"} is the current recommended option,
  especially for spectral counting proteomics data sets that are typically
  characterized by a high data granularity or for microarray data sets with a
  small number of replicates. In both cases, there are chances for data values
  for a same gene or protein to be identical across replicates (and therefore
  with zero standard deviation) by chance alone. Note that setting
  \code{trimAllZeroRows=TRUE} does not guarantee that there will be no rows with
  zero mean or zero standard deviation.

  If argument \code{fittingEval} is set to \code{TRUE}, a graphical control of
  the goodness of the \bold{PLGEM} fitting is produced and a plot containing
  four panels is generated. The top-left panel shows the power law,
  characterized by a \sQuote{SLOPE} and an \sQuote{INTERCEPT}. The top-right
  panel represents the distribution of model residuals. The bottom-left reports
  the contour plot of ranked residuals. The bottom-right panel finally shows the
  relationship between the distribution of observed residuals and the normal
  distribution. A good fit normally gives a horizontal symmetric rank-plot and a
  near normal distribution of residuals.
}
\value{
  A \code{list} of six elements (see Details):
  \item{SLOPE}{the slope of the fitted PLGEM.}
  \item{INTERCEPT}{the intercept of the fitted PLGEM.}
  \item{DATA.PEARSON}{the Pearson correlation coefficient between the
    \eqn{\log{(sd)}}{log(sd)} and the \eqn{\log{(mean)}}{log(mean)} in the
    original data.}
  \item{ADJ.R2.MP}{the adjusted \eqn{r^{2}}{r^2} of PLGEM fitted on the
    modelling points.}
  \item{COVARIATE}{a \code{character} indicating the covariate used for
    fitting.}
  \item{FIT.CONDITION}{a \code{character} indicating the condition used for
    fitting.}
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
  \code{\link{plgem.obsStn}}, \code{\link{plgem.resampledStn}},
  \code{\link{plgem.pValue}}, \code{\link{plgem.deg}}, \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset, fittingEval=TRUE)
  as.data.frame(LPSfit)
}
\keyword{models}
