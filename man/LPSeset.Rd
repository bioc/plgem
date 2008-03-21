\name{LPSeset}
\alias{LPSeset}
\docType{data}
\title{ExpressionSet for Testing PLGEM}
\description{
This ExpressionSet object contains the same microarray data used in References for PLGEM set up and validation.
Briefly, it contains the result of a total of 6 hybridizations on MG-U74Av2 Affymetrix GeneChips.
MAS5 scaled data represent two experimental conditions: the baseline condition (C) contains data of immature murine
dendritic cells (4 replicates), the other condition (LPS) contains data of the same cells stimulated
4 hours with LPS (2 replicates).
This dataset is part of the data used in references.
}
\usage{data(LPSeset)}
\format{An object of class \code{\linkS4class{ExpressionSet}}.}
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
\examples{
data(LPSeset)
LPSfit<-plgem.fit(data = LPSeset, fittingEval = TRUE)
}
\keyword{models}
