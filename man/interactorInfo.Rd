\name{interactorInfo}
\alias{interactorInfo}
\title{ Interactor info in a matrix }
\description{
  The function returns the essential information of interactors in a
  matrix. Xrefs are left out since they have arbitrary numbers of
  annotation and cannot be summarized into a matrix.
}
\usage{
interactorInfo(x)
}
\arguments{
  \item{x}{ An object which contains psimi25Interactor information, for
    example objects of \code{\link{psimi25InteractionEntry-class}},
    \code{\link{psimi25Graph-class}}, \code{\link{psimi25Hypergraph-class}} or \code{\link{psimi25ComplexEntry-class}}
    }.
}
\value{
  A matrix of interactor information, each row represents an interactor
  The columns are
  \item{sourceDb }{source database}
  \item{sourceId }{source database index}
  \item{shortLabel}{short label assigned by the source database}
  \item{uniprotId}{UniProt ID, NA if not available}
  \item{organismName}{the organism of the interactor protein}
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{\code{\link{xref}}, \code{\link{availableXrefs}}}
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)

hprdInteractorInfo <- interactorInfo(hprdSet)
}
\keyword{ models }

