\name{parsePsimi25Interaction}
\alias{parsePsimi25Interaction}
\alias{parsePsimi25Complex}
\title{ Parsing PSI-MI 2.5 XML documents into interactions }
\description{
  The PSI-MI 2.5 XML format is used widely by many repositories to
  record protein-protein interaction data as well as protein complex
  data. This functions parse such files into interactions or complexes.

  \code{parsePsimi25Interaction} is the parser for interaction data and
  \code{parsePsimi25Complex} is the parser for complex data.
}
\usage{
parsePsimi25Interaction(psimi25file, psimi25source, verbose=TRUE)
parsePsimi25Complex(psimi25file, psimi25source, verbose=FALSE)
}

\arguments{
  \item{psimi25file}{ character, file name or URL of the XML document }
  \item{psimi25source}{ A supported data repository source, see also
    \code{\link{psimi25Source-class}} }
  \item{verbose}{ logical, whether the parsing state should be displayed verbosely.}
}
\value{
  \code{psimi25Interaction} returns a list of \code{psimi25InteractionEntry}
  objects, each represents one \code{entry} in the XML document
  \code{psimi25Complex} returns a \code{psimi25ComplexEntry} objects,
  representing the complex data from one XML document.
}
\references{
 \describe{
  \item{PSI-MI XML v2.5 data exchange
    format}{\url{http://www.psidev.info/index.php?q=node/60}}
  \item{BioGRID}{\url{http://www.thebiogrid.com}}
  \item{DIP}{\url{http://dip.doe-mbi.ucla.edu/}}
  \item{HPRD}{\url{http://www.hprd.org/}}
  \item{IntAct}{\url{http://www.ebi.ac.uk/intact}}
  \item{MINT}{\url{http://mint.bio.uniroma2.it/mint}}
  \item{MIPS/CORUM}{\url{http://mips.gsf.de/}}
 }
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang
  <tchiang@ebi.ac.uk> }
\seealso{
  \code{\link{psimi25Interaction-class}},
  \code{\link{psimi25InteractionEntry-class}},
  \code{\link{psimi25Complex-class}}
  \code{\link{psimi25ComplexEntry-class}},
}
\examples{
# parse interaction data
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

gridxml <- file.path(xmlDir, "biogrid_200804_test.xml")
gridSet <- parsePsimi25Interaction(gridxml, BIOGRID.PSIMI25)

intactxml <- file.path(xmlDir, "intact_2008_test.xml")
intactSet <- parsePsimi25Interaction(intactxml, INTACT.PSIMI25, verbose=TRUE)

# parse complex data
intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexSet <- parsePsimi25Complex(intactComplexxml, INTACT.PSIMI25)

}
\keyword{ models }
