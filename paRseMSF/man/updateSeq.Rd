\name{updateSeq}
\alias{updateSeq}
\title{Readable Peptide Modification info}
\description{
  This is an internal function that updates each peptide modification information
  based on a vector of values in 
  c(PeptideID, Sequence, Position, ModificationName, TerminalModificationName, ConfidenceLevel)
}
\usage{
  updateSeq(x)
}
\arguments{
  \item{x}{
    Vector of values 
    c(PeptideID, Sequence, Position, ModificationName, TerminalModificationName, ConfidenceLevel)
  }
}
\value{
  Returns a vector of values:
  c(PeptideID, Sequence, ModifiedSequence, Modifications, Confidence)
}
\author{
  Ashoka Polpitiya
}
\seealso{
    \code{\link{protModsTable}}, \code{\link{scoreTable}}, \code{\link{mtzTable}}, \code{\link{protAnnotations}}, \code{\link{protModsTable}}, \code{\link{aucTable}}
}
\keyword{misc}

