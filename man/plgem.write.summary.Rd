\name{plgem.write.summary}
\alias{plgem.write.summary}
\title{
  Write a list of differentially expressed genes/proteins to the disk
}
\description{
  This function writes the list of differentially expressed genes/proteins
  obtained via a call to either \code{\link{plgem.deg}} or
  \code{\link{run.plgem}} to a series of files in the current working directory.
}
\usage{
  plgem.write.summary(x, verbose=FALSE)
}
\arguments{
  \item{x}{\code{list}; the output of either \code{\link{plgem.deg}} or
    \code{\link{run.plgem}}, i.e. a list of list(s) of named vectors.}
  \item{verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  The gene or protein lists are written to the current working directory, using
  conveniently chosen filenames that reflect the specific comparisons that were
  performed (i.e. which experimental condition was compared to the baseline) and
  the specific significance threshold that were used in the DEG selection step.
}
\value{
  The function returns no value. It is called for its side effect to write files
  to the working directory.
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
  \code{\link{plgem.deg}}, \code{\link{run.plgem}}
}
\examples{
  ## Not run:
  data(LPSeset)
  LPSdegList <- run.plgem(LPSeset)
  plgem.write.summary(LPSdegList, verbose=TRUE)
  ## End(Not run)
}
\keyword{models}
