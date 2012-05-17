\name{protAnnotations}
\alias{protAnnotations}
\title{Protein Annotations}
\description{

}
\usage{
protAnnotations(msfFile, minConf = "High")
}
\arguments{
\item{msfFile}{
    MSF file from a Thermo Proteome Discoverer Run
  }
  \item{minConf}{
    Minimum confidence level for identifications. Can be either \sQuote{Low},
    \sQuote{Medium}, or \sQuote{High}.
  }
}
\value{

}
\author{
Ashoka Polpitiya
}
\seealso{
  \code{\link{protModsTable}}, \code{\link{scoreTable}}, \code{\link{mtzTable}}, \code{\link{protModsTable}}, \code{\link{aucTable}}
}
\examples{
\dontrun{M <- protAnnotations(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{misc}
