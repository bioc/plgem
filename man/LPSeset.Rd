\name{LPSeset}
\alias{LPSeset}
\docType{data}
\title{ExprSet for Testing PLGEM}
\description{
This exprSet contains the same microarray data used in References for PLGEM set up and validation. Briefly, it
contains the result of a total of 6 hybridizations on MG-U74Av2 Affymetrix GeneChips.
MAS5 scaled data represent two experimental conditions: the baseline condition (C) contains data of immature murine
dendritic cells (4 replicates), the other condition (LPS) contains data of the same cells stimulated
4 hours with LPS (2 replicates).
This dataset is part of the data used in references.
}
\usage{data(LPSeset)}
\format{expression Set}
\references{N. Pavelka et al., BMC Bioinformatics, 2004 Dec 17;5(1):203; \url{http://www.genopolis.it}}
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
}
\keyword{datasets}
