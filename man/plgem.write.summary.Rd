\name{plgem.write.summary}
\alias{plgem.write.summary}
\title{
  Write the Result of a PLGEM Analysis to the Working Directory
}
\description{
  This function writes the output of function \code{\link{plgem.deg}} or
  \code{\link{run.plgem}} to a series of files in the current working directory.
}
\usage{
  plgem.write.summary(x, prefix=NULL, verbose=FALSE)
}
\arguments{
  \item{x}{\code{list}; the output of either \code{\link{plgem.deg}} or
    \code{\link{run.plgem}} (see corresponding help papes for details).}
  \item{prefix}{optional \code{character} to use as a prefix of the file names
    to be written.}
  \item{verbose}{\code{logical}; if \code{TRUE}, comments are printed out while
    running.}
}
\details{
  This function writes three types of files to the current working directory:
  
  1) A comma-separated text file containing the PLGEM fitting parameters;
  
  2) A comma-separated text file containing the observed STN values and their
  associated p-values for all performed comparisons; (STN values and p-values
  from different comparisons appear in different columns, with column headers
  reflecting the underlying comparison)
  
  3) One or more plain text files containing the identifiers of the significant
  genes or proteins, with filenames reflecting the specific comparisons that
  were performed (i.e. which experimental condition was compared to the
  baseline) and the specific significance threshold that were used in the DEG
  selection step.
  
  Before files are written, the function checks for existence of files with
  identical names in the working directory and prompts the user to decide
  whether to abort the writing process or to overwrite the existing files.
}
\value{
  The function returns no value. It is called for its side effect to write files
  to the working directory.
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
  \code{\link{plgem.deg}}, \code{\link{run.plgem}}
}
\examples{
\dontrun{
  data(LPSeset)
  LPSdegList <- run.plgem(LPSeset, fitting.eval=FALSE)
  plgem.write.summary(LPSdegList, prefix="test", verbose=TRUE)}
}
\keyword{models}
