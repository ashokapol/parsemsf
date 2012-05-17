\name{mtzTable}
\alias{mtzTable}
\title{Mass, Retention time, Charge information Table}
\description{
  Returns a data frame with identified peptide sequences, their modifications and 
  adds mass, retention time/scan and charge information.
}
\usage{
mtzTable(msfFile, minConf = "High")
}
\arguments{
\item{msfFile}{
    MSF file from a Thermo Proteome Discoverer Run
  }
  \item{minConf}{
    Minimum confidence level for identifications. Can be either \sQuote{Low},
    \sQuote{Medium}, or \sQuote{High}.
  }
  \item{proteins}{
    Specify if you want the corresponding proteins listed. (TRUE or FALSE)
  }
}
\value{
}
\author{
Ashoka Polpitiya
}
\seealso{
\code{\link{protModsTable}}, \code{\link{scoreTable}}, \code{\link{aucTable}}, \code{\link{protAnnotations}}, \code{\link{protModsTable}}
}
\examples{
\dontrun{M <- mtzTable(msfFile="D:/Research/R/msfRead/test_dimethyl.msf", minConf="High")}
}
\keyword{misc}
