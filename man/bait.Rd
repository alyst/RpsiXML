\name{bait}
\alias{bait}
\alias{prey}
\alias{participant}
\alias{inhibitor}
\alias{pubmedID}
\alias{confidenceValue}
\alias{neutralComponent}
\title{ Extract bait, prey, participant, inhibitor, pubmed, confidence value, interaction type, 
or neutral component information from an object of psimi25Interaction-class }
\description{
  The functions return bait/prey UniProt identifier of the given
  \code{psimi25Interaction} object.
}
\usage{
bait(x,...)
prey(x,...)
participant(x,...)
inhibitor(x,...)
pubmedID(x,...)
confidenceValue(x,...)
neutralComponent(x,...)
}
\arguments{
  \item{x}{ An object of \code{\link{psimi25Interaction-class}}, see
    example }
  \item{...}{ Other parameters to control the identifier returned, not
    implemented yet}
}
\value{
  The source database identifier is returned.
}
\references{The UniProt database \url{http://www.expasy.uniprot.org/}}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk>}
\seealso{ \code{\link{psimi25Interaction-class}} }
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

gridxml <- file.path(xmlDir, "biogrid_200804_test.xml")
gridSet <- parsePsimi25Interaction(gridxml, BIOGRID.PSIMI25)

interExp <- interactions(gridSet)[[1]]
bait(interExp)
prey(interExp)
}
\keyword{ models }

