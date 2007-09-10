\name{plgem.fit}
\alias{plgem.fit}
\title{PLGEM Fitting and Evaluation}
\description{
Function for fitting and evaluating goodness of fit of \bold{PLGEM} on a 
\sQuote{data} ExpressionSet, using the set of replicated samples identified by the \sQuote{fit.condition} condition of the
\sQuote{covariateNumb} covariate;
partitioning the range of expression values in \sQuote{p} intervals,
using the \sQuote{q}-th quantile of expression value standard deviations.  
}
\usage{plgem.fit(data, covariateNumb = 1, fit.condition = 1, p = 10, q = 0.5, fittingEval = FALSE,
plot.file = FALSE, verbose = FALSE)
}
\arguments{
  \item{data}{an object of class \sQuote{ExpressionSet}.}
  \item{covariateNumb}{number; the covariate used to determine on which samples to fit plgem.}
  \item{fit.condition}{number; the condition used for plgem fitting, according to the order of unique values of the covariateNumb covariate.}
  \item{p}{number of intervals used to partition the expression value range.}
  \item{q}{number in [0,1]; the quantile of standard deviation used for PLGEM fitting.}
  \item{fittingEval}{logical; if TRUE, the fitting is evaluated generating a diagnostic plot.}
  \item{plot.file}{logical; if TRUE, a png file is written on the current working directory.}
  \item{verbose}{logical; if TRUE, comments are printed out while running.}
}
\details{
\sQuote{plgem.fit} fits PLGEM on an expression set and eventually evaluates goodness of fit.
This Power Law Global Error Model aims to find the
mathematical relationship between standard deviation and mean expression values in a set of replicated
microarray samples, according to a power law:

ln(modeledSpread) = PLGEMslope * ln(mean) + PLGEMintercept

The covariateNumb covariate (the 1st one by default) of the phenoData of the ExpressionSet \sQuote{data} is
expected to contain the necessary information about te experimental design. The values of this covariate
must be sample labels, that have to be identical for samples to be treated as replicates. 
This function returns \sQuote{SLOPE} and \sQuote{INTERCEPT} of this power law; moreover it returns the
Pearson's coefficient of correlation \sQuote{DATA.PEARSON} of the linear model fitted on the original data,
as well as the adjusted R squared \sQuote{ADJ.R2.MP} of the linear model fitted on the modelling points.

If argument \sQuote{fittingEval} is TRUE, a graphical control of the goodness of the plgem fitting is produced and 
a plot containing four panels is generated. The top-left panel shows the power law, characterized by \sQuote{SLOPE}
and \sQuote{INTERCEPT}. The top-right panel represents the distribution of model residuals. The bottom-left
reports the contour plot of ranked residuals. The bottom-right panel finally shows the relationship between
the distribution of observed residuals and the normal distribution.
The goodness of the fit is principally judged by an horizontal symmetric rank-plot and a near normal
distribution of residuals.
}
\value{
\sQuote{plgem.fit} returns a list of five numbers (see details): 
  \item{SLOPE}{the slope of the fitted PLGEM.}
  \item{INTERCEPT}{the intercept of the fitted PLGEM.}
  \item{DATA.PEARSON}{the Pearson correlation coefficient of the linear model fitted on the original data.}
  \item{ADJ.R2.MP}{the adjusted R squared of PLGEM fitted on the modelling points.}
  \item{FIT.CONDITION}{the condition used for fitting PLGEM.}
}
\references{N. Pavelka et al., BMC Bioinformatics, 2004 Dec 17;5(1):203; \url{http://www.genopolis.it} }
\author{Mattia Pelizzola \email{mattia.pelizzola@gmail.com} and Norman Pavelka \email{NXP@stowers-institute.org}}
\seealso{\code{\link{plgem.obsStn}},\code{\link{plgem.resampledStn}},\code{\link{plgem.deg}},\code{\link{run.plgem}} }
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
}
\keyword{models}
