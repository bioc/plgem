\name{plgem.fit}
\alias{plgem.fit}
\title{PLGEM Fitting and Evaluation}
\description{
Function for fitting and evaluating goodness of fit of \bold{PLGEM} on a 
\sQuote{data} exprSet, using the condition \sQuote{fit.condition} containing replicates;
partitioning the range of expression values in \sQuote{p} intervals,
using the \sQuote{q}-th quantile of expression value standard deviations.  
}
\usage{plgem.fit(data, fit.condition, p = 10, q = 0.5, fittingEval = FALSE,
plot.file = FALSE, verbose = FALSE)
}
\arguments{
  \item{data}{an object of class \sQuote{exprSet} with a \sQuote{conditionName} covariate, see details.}
  \item{fit.condition}{number; the condition used for plgem fitting, according to the order of unique values of conditionName covariate.}
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

The exprSet \sQuote{data} must have a phenoData slot with a covariate called \sQuote{conditionName}. The values of this covariate
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
\author{Mattia Pelizzola \email{mattia.pelizzola@unimib.it} and Norman Pavelka \email{norman.pavelka@unimib.it}}
\seealso{\code{\link{plgem.obsStn}},\code{\link{plgem.resampledStn}},\code{\link{plgem.deg}},\code{\link{run.plgem}} }
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
}
\keyword{models}
