\name{LPSeset}
\alias{LPSeset}
\docType{data}
\title{
  ExpressionSet for Testing PLGEM
}
\description{
  This \code{ExpressionSet} object contains a subset of the microarray data used
  in References for \bold{PLGEM} set-up and validation. Briefly, it contains
  normalized gene expression values from a total of 6 hybridizations on
  MG-U74Av2 Affymetrix GeneChip microarrays. Two experimental conditions are
  represented in this dataset: the baseline condition (\sQuote{C}) contains data
  of immature murine dendritic cells (4 replicates); the other condition
  (\sQuote{LPS}) contains data of the same cells stimulated for 4 hours with LPS
  (2 replicates).
}
\usage{
  data(LPSeset)
}
\format{
  An object of class \code{\link[Biobase:class.ExpressionSet]{ExpressionSet}}.
}
\references{
  Pavelka N, Pelizzola M, Vizzardelli C, Capozzoli M, Splendiani A, Granucci F,
  Ricciardi-Castagnoli P. A power law global error model for the identification
  of differentially expressed genes in microarray data. BMC Bioinformatics. 2004
  Dec 17;5:203.; \url{http://www.biomedcentral.com/1471-2105/5/203}
}
\author{
  Mattia Pelizzola \email{mattia.pelizzola@gmail.com}

  Norman Pavelka \email{normanpavelka@gmail.com}
}
\examples{
  data(LPSeset)
  library(Biobase)
  head(exprs(LPSeset))
}
\keyword{models}
