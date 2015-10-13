\name{interactions}
\alias{interactions}
\title{ List interactions in the given psimi25InteractionEntry object }
\description{
  List interactions from an object of the
  \code{\link{psimi25InteractionEntry-class}} object
}
\usage{
interactions(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of \code{\link{psimi25InteractionEntry-class}},
    see example }
}
\value{
  A list of interactions
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk> }
\seealso{ \code{\link{psimi25InteractionEntry-class}} }
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)

interactions(hprdSet)
}

\keyword{ models}
