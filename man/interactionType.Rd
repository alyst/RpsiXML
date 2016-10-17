\name{interactionType}
\alias{interactionType}
\alias{interactionType,psimi25InteractionElementType-method}
\title{ Type of the interaction }
\description{
  Return the interaction type of the \code{psimi25Interaction}
}
\usage{
interactionType(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ An object of \code{\link{psimi25Interaction-class}} }
}
\value{
  A character string representing the interaction type
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk> }
\seealso{ \code{\link{psimi25Interaction-class}} }
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

hprdxml <- file.path(xmlDir, "hprd_200709_test.xml")
hprdSet <- parsePsimi25Interaction(hprdxml, HPRD.PSIMI25)

interExp <- interactions(hprdSet)[[1]]
interExpTyp <- interactionType(interExp)
}
\keyword{ models }

