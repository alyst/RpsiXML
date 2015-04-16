\name{interactors}
\alias{interactors}
\alias{interactors-methods}
\alias{interactors<-}
\alias{interactors<--methods}
\alias{interactors,interactorListBase-method}
\alias{interactors<-,interactorListBase,list-method}
\alias{interactors<-,interactorListBase,NULL-method}
\alias{numInteractors}
\alias{numInteractors-methods}
\alias{numInteractors,interactorListBase-method}
\title{ Get information of interactors of the given object}
\description{
  \code{interactors} gets the list of interactors stored in the given
  object. One can also assign a list to replace old interactors.

  \code{numInteractors} returns the length of interactors.
}
\usage{
interactors(x)
interactors(x) <- value
numInteractors(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of the subclass of \code{interactorListBase}}
  \item{value}{A list of interactors}
}
\value{
  The getting method returns a list of interactors. The setting method
  mutates the object. The count method returns the list length as integer.
}
\author{Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk>}
\seealso{ \code{\link{psimi25InteractionEntry-class}},
  \code{\link{psimi25ComplexEntry-class}}
}
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)
interactors(hprdSet)

## Not run:
intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexSet <- parsePsimi25Complex(intactComplexxml,
INTACT.PSIMI25)

interactors(intactComplexSet)
numInteractors(intactComplexSet)
## End(Not run)
}

