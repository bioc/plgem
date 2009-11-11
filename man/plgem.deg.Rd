\name{plgem.deg}
\alias{plgem.deg}
\title{
  Selection of Differentially Expressed Genes/Proteins With PLGEM
}
\description{
  This function selects differentially expressed genes/proteins (DEG) at a given
  significance level, based on observed \bold{PLGEM} signal-to-noise ratio (STN)
  values (typically obtained via a call to \code{\link{plgem.obsStn}}) and
  pre-computed p-values (typically obtained via a call to
  \code{\link{plgem.pValue}}).
}
\usage{
  plgem.deg(observedStn, plgemPval, delta=0.001, verbose=FALSE)
}
\arguments{
  \item{observedStn}{\code{list} containing a \code{matrix} of observed
    PLGEM-STN values; output of function \code{\link{plgem.obsStn}}.}
  \item{plgemPval}{\code{matrix} of p-values; output of function
    \code{\link{plgem.pValue}}.}
  \item{delta}{numeric vector; the significance level(s) to be used for the
    selection of DEG; value(s) must be between 0 and 1 (excluded).}
  \item{verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  This function allows for the selection of DEG by setting a significance
  cut-off on pre-calculated p-values. The significance level \code{delta}
  roughly represents the false positive rate of the DEG selection, e.g. if a
  \code{delta} of 0.001 is chosen in a microarray dataset with 10,000 genes
  (none of which is truly differentially expressed), on average 10
  genes/proteins are expected to be selected by chance alone.
}
\value{
  A \code{list} of four elements:
  \item{fit}{the input \code{plgemFit}.}
  \item{PLGEM.STN}{the input \code{matrix} of observed PLGEM-STN values (see
    \code{\link{plgem.obsStn}} for details).}
  \item{p-value}{the input \code{matrix} of p-values (see
    \code{\link{plgem.pValue}} for details).}
  \item{significant}{a \code{list} with a number of elements equal to the number
    of different significance levels (\code{delta}) used as input. Each element
    of this list is again a list, whose number of elements correspond to the
    number of performed comparisons (i.e. the number of conditions in the
    starting \code{ExpressionSet} minus the baseline). Each of these second
    level elements is a \code{character} vector of significant gene/protein
    names that passed the statistical test at the corresponding significance
    level.}
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
  \code{\link{plgem.fit}}, \code{\link{plgem.obsStn}},
  \code{\link{plgem.resampledStn}}, \code{\link{plgem.pValue}},
  \code{\link{run.plgem}}
}
\examples{
  data(LPSeset)
  LPSfit <- plgem.fit(data=LPSeset, fittingEval=FALSE)
  LPSobsStn <- plgem.obsStn(data=LPSeset, plgemFit=LPSfit)
  set.seed(123)
  LPSresampledStn <- plgem.resampledStn(data=LPSeset, plgemFit=LPSfit)
  LPSpValues <- plgem.pValue(LPSobsStn, LPSresampledStn)
  LPSdegList <- plgem.deg(observedStn=LPSobsStn, plgemPval=LPSpValues,
    delta=0.001)
}
\keyword{models}
