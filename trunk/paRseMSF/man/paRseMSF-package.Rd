\name{paRseMSF-package}
\alias{paRseMSF-package}
\alias{paRseMSF}
\docType{package}
\title{Retrieving data from a Proteome Discoverer Results File (MSF extension)}
\description{
  Retrieving data from a Proteome Discoverer Results File (MSF extension)
}
\details{
\tabular{ll}{
Package: \tab paRseMSF\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-05-16\cr
License: \tab Apache License 2.0\cr
LazyLoad: \tab yes\cr
}
}
\author{
Ashoka Polpitiya
Maintainer: Ashoka Polpitiya <ashoka@tgen.org>
}
\keyword{misc}
\seealso{
\code{\link{protModsTable}}, \code{\link{scoreTable}}, \code{\link{mtzTable}}, \code{\link{protAnnotations}}, \code{\link{protModsTable}}, \code{\link{aucTable}}
}
\examples{
\dontrun{M <- protModsTable(msfFile="test_msf.msf", minConf="High", proteins=TRUE)}
}
