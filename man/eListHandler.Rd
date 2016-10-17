\name{eListHandler}
\alias{eListHandler}
\alias{iListHandler}
\title{ xmlEventParse handlers for PSI-MI XML documents}
\description{
xmlEventParse handlers for PSI-MI XML documents
}
\usage{
eListHandler()
iListHandler()
}
%- maybe also 'usage' for other objects documented here.
\details{
A dump() method is supplied to deliver the list.
}
\value{
 list; see examples for structure
}
\author{ Vince Carey <stvjc@channing.harvard.edu> }
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")
fn <- file.path(xmlDir, "hprd_200709_test.xml")
XML::xmlEventParse(fn, eListHandler())$dump()
}
\keyword{ models }
