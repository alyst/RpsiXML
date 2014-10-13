\name{availableXrefs}
\alias{availableXrefs}
\alias{xref}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Get cross reference(s) of interactors from PSI-MI 2.5 files  }
\description{
  PSI-MI 2.5 data exchange format encourages data providers to encode
  cross referenes of the interactors into the file, which prove valuable
  for annotating the data. These functions help finding out the available cross references of
  interactors. 
}
\usage{
availableXrefs(x,...)
xref(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    for \code{availableXrefs}, 'x' can either be an object of
    \code{\link{psimi25Interactor-class}}, 
    \code{link{psimi25InteractionEntry}} or a list of
    \code{\link{psimi25Interactor-class}} objects.

    for \code{xref}, 'object' should be an object of \code{\link{psimi25Interactor-class}}
  }
  \item{object}{
    New interface of RpsiXML uses \sQuote{object} uniformly. 
    }
  \item{...}{intersect, logical option. If set to \code{TRUE} then only
    the cross references available for every interactor are returned
    (intersected), otherwise all the cross references are unioned and
    returned. See the examples}
}
\details{
  If \code{psimi25InteractionEntry} is provided as the parameter of
  \code{availableXrefs}, an option named "intersect" can be set to
  extract only those cross references intersected among all the
  interactors, namely the ones assigned in every interactor. Please see
  the examples.
}
\value{
  \item{xref }{Returns a matrix with two columns: \code{db} (external database
    name) and \code{id} (external database index) }
  \item{availableXrefs }{A vector of characters, the names of external
    databases referenced in the file}
}
\author{ Jitao David Zhang }
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)
hprdInteractors <- interactors(hprdSet)

availableXrefs(hprdSet)
xref(hprdInteractors[[1]])
availableXrefs(hprdInteractors[[1]])
}
\keyword{ models }
