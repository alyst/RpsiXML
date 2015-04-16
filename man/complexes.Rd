\name{complexes}
\alias{complexes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Extract complexes from psimi25ComplexEntry }
\description{
  Extract complexes from an object of \code{\link{psimi25ComplexEntry-class}}
}
\usage{
complexes(x)
}
\arguments{
  \item{x}{ An object of \code{\link{psimi25ComplexEntry-class}} }
}
\value{
  A list of \code{psimi25Complex}
}
\references{ PSI-MI 2.5 XML data exchange format \url{http://www.psidev.info/index.php?q=node/60} }
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk> }
\seealso{ \code{\link{psimi25Complex-class}} }
\examples{
## Not run:
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexSet <- parsePsimi25Complex(intactComplexxml,
INTACT.PSIMI25)

complexes(intactComplexSet)
## End(Not run)
}
\keyword{ models }
