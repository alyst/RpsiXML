\name{complexName}
\alias{complexName}
\alias{members}
\alias{attributesList}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Accessor functions for complex }
\description{
  These functions are used to extract useful information of complex in
  the form of \code{\link{psimi25Complex-class}} object.
}
\usage{
complexName(x,...)
members(x)
attributesList(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ An object of \code{\link{psimi25Complex-class}} }
  \item{...}{Not implemented yet}
}
\details{
  See examples
}
\value{
  \item{complexName}{Returns the name of the complex in characters}
  \item{members}{A data frame of protein members building the
    complex and their information}
  \item{attributesList}{A list of \code{psimi25Attribute} objects,
    recording the attribute name, name accession and value.}
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk> }
\seealso{\code{\link{psimi25Complex-class}}}
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexSet <- parsePsimi25Complex(intactComplexxml, INTACT.PSIMI25)
complexSample <- complexes(intactComplexSet)[[2]]

complexName(complexSample)
attributesList(complexSample)
members(complexSample)
}
\keyword{ models }
