\name{translateID}
\alias{translateID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Finds identifiers of a given object}
\description{
  see \code{\link{translateID-methods}}
}
\usage{
translateID(r, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{An object of psimi25Graph, psimi25Hypergraph,
    psimi25Interactor or a list of \code{psimi25Interactor}}
  \item{...}{the symbol of the ID to translate}
}
\details{
  see \code{\link{translateID-methods}}
}
\value{
  The object of the same class as the parameter
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{see \code{\link{translateID-methods}} }
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)

it <- interactors(hprdSet)[[1]]
translateID(it, "uniprot")
translateID(it, "entrezgene")

##Not run
intactxml <- file.path(xmlDir, "intact_2008_test.xml")
intactSet <- parsePsimi25Interaction(intactxml, INTACT.PSIMI25)
intactGraph <- psimi25XML2Graph(intactxml, INTACT.PSIMI25)
intactGraphNew <- translateID(intactGraph,"sourceId")## translate the nodes of the graph to another identifier
                              
intactSetInteractors <- interactors(intactSet)
intactXrefExample <- xref(intactSetInteractors[[1]])
translateID(intactSetInteractors,"intact")
translateID(intactSetInteractors[[1]],"intact")

intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexSet <- parsePsimi25Complex(intactComplexxml, INTACT.PSIMI25)
intactComplexGraph <- psimi25XML2Graph(intactComplexxml, INTACT.PSIMI25, type="complex")
translateID(intactComplexGraph, "intact", "P49432")
translateID(intactComplexGraph, "intact", NA)
## End(Not run)

}
\keyword{ models }