\name{protModsTable}
\alias{protModsTable}
\title{Aminoacid Modifications Table}
\description{
  Returns a data frame with identified peptide sequences, their modifications and 
  corresponding proteins.
}
\usage{
protModsTable(msfFile, minConf = "High", proteins = TRUE)
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
  A data frame with peptides, their modifications and the proteins.
}
\author{
  Ashoka Polpitiya
}
\seealso{
  \code{\link{scoreTable}}, \code{\link{mtzTable}}, \code{\link{protAnnotations}}, \code{\link{protModsTable}}, \code{\link{aucTable}}
}
\examples{
  \dontrun{M <- protModsTable(msfFile="test_msf.msf", minConf="High", proteins=TRUE)}
}
\keyword{misc}
