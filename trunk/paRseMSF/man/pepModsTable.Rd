\name{pepModsTable}
\alias{pepModsTable}
\title{Internal Function to get Peptide Modifications}
\usage{
pepModsTable(msfFile, minConf = "High")
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
  \code{\link{protModsTable}}, \code{\link{scoreTable}}, \code{\link{mtzTable}}, \code{\link{protAnnotations}}, \code{\link{protModsTable}}, \code{\link{aucTable}}
}
\keyword{misc}
