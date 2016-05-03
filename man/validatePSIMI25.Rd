\name{validatePSIMI25}
\alias{validatePSIMI25}
\title{ Validating PSI-MI 2.5 file with MIF25 XML schema }
\description{
  The function validates given PSI-MI 2.5 file with MIF25 XML schema
  provided by the Molecular Interactions Workgroup of HUPO Proteomics Standards Initiative 
}
\usage{
validatePSIMI25(file,
schema = system.file("extdata/schemas/MIF25.xsd", package = "RpsiXML"),
ignore.stderr = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{ The name or the URL of the file to be validated }
  \item{schema}{ The schema file of PSI-MI 2.5 file format by default }
  \item{ignore.stderr}{ Whether to print out errors in the console }
}
\details{
  We advice to set "ignore.stderr=FALSE" for the first time of
  validating. The error number of the file will be printed in the
  console. When the number is not 0, one can set "ignore.stderr=TRUE"
  and see the errors produced during the validation.

}
\value{
  The number of errors in validating the file
}
\references{
  libxml2 library \url{http://xmlsoft.org/}
  PSI-MI 2.5 XML schema \url{http://psidev.sourceforge.net/mi/rel25/src/MIF25.xsd}
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk> }
\note{
  This function sofar depends on the tool "xmllint" bundled with libxml2
  library. We plan to compile it later to run independent of the tool.
}
\examples{
okFile <- system.file("extdata/psi25files/intact_2008_test.xml",
                      package="RpsiXML")
validatePSIMI25(okFile)
}
\keyword{ models }

