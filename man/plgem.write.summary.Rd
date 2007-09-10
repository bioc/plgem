\name{plgem.write.summary}
\alias{plgem.write.summary}
\title{Writing Differentially Expressed Genes on the Disk}
\description{
This function writes the differentially expressed genes contained in the list obtained by either \sQuote{plgem.deg} or \sQuote{run.plgem} functions.
}
\usage{plgem.write.summary(x, verbose=FALSE)}
\arguments{
  \item{x}{list; the output of either the \sQuote{plgem.deg} or the \sQuote{run.plgem} functions, i.e. a list of list(s) of named vectors}
  \item{verbose}{logical; if TRUE, comments are printed out while running.}
}
\details{
The gene lists are written on the current working directory.
}
\value{}
\references{N. Pavelka et al., BMC Bioinformatics, 2004 Dec 17;5(1):203; \url{http://www.genopolis.it} }
\author{Mattia Pelizzola \email{mattia.pelizzola@gmail.com} and Norman Pavelka \email{NXP@stowers-institute.org}}
\seealso{\code{\link{plgem.deg}},\code{\link{run.plgem}} }
\examples{}
\keyword{models}
