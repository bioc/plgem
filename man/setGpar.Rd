\name{setGpar}
\alias{setGpar}
\title{Set graphical parameters for PLGEM fitting evaluation plots}
\description{
  Function to set graphical parameters for \bold{PLGEM} fitting evaluation
  plots produced by function \code{plgem.fit}.
}
\usage{
  setGpar(minLnM=NULL, maxLnM=NULL, minLnSD=NULL, maxLnSD=NULL,
    minRes=NULL, maxRes=NULL)
}
\arguments{
  \item{minLnM}{minimum 'ln(mean)' value in upper left plot.}
  \item{maxLnM}{maximum 'ln(mean)' value in upper left plot.}
  \item{minLnSD}{minimum 'ln(sd)' value in upper left plot.}
  \item{maxLnSD}{maximum 'ln(sd)' value in upper left plot.}
  \item{minRes}{minimum 'residual' value in upper right plot.}
  \item{maxRes}{maximum 'residual' value in upper right plot.}
}
\details{
  A call to \code{setGpar()} is the recommended way to set graphical parameters
  in \code{plgem.fit}. If parameters are left unspecified, suitable defaults
  will be found by \code{plgem.fit}.
}
\value{
  A \code{list} of six elements:
  \item{minLnM}{minimum 'ln(mean)' value in upper left plot.}
  \item{maxLnM}{maximum 'ln(mean)' value in upper left plot.}
  \item{minLnSD}{minimum 'ln(sd)' value in upper left plot.}
  \item{maxLnSD}{maximum 'ln(sd)' value in upper left plot.}
  \item{minRes}{minimum 'residual' value in upper right plot.}
  \item{maxRes}{maximum 'residual' value in upper right plot.}
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
  Norman Pavelka \email{normanpavelka@gmail.com}
}
\seealso{
  \code{\link{plgem.fit}}, \code{\link{run.plgem}}
}
\examples{
  setGpar(minLnM=-1, maxLnM=8)
}
\keyword{models}
