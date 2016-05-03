\name{uniprot}
\alias{uniprot}
\title{ The UniProt Identifier in the PSI-MI 2.5 XML file }
\description{
  A character string representing the UniProt identifier in primaryRef
  or secondaryRef attributes
}
\usage{
uniprot(x)
}
\arguments{
  \item{x}{ An object of \code{\link{psimi25Source-class}} }
}
\value{
  A character string of the identifier
}
\author{ Tony Chiang <tchiang@ebi.ac.uk>, Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{ \code{\link{psimi25Source-class}} }
\examples{
uniprot(HPRD.PSIMI25)
}
\keyword{ models }
