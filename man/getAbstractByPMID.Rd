\name{getAbstractByPMID}
\alias{getAbstractByPMID}
\title{A function to obtain the abstract information via a Pubmed ID}
\description{
  This function takes a character vector of pubmed IDs and returns a list of 
  pubMedAbst objects indexed by each ID.
}
\usage{
getAbstractByPMID(pmID)
}
\arguments{
  \item{pmID}{ A chacater vector of pubmed IDs}
}
\value{
 A list of pubMedAbst objects.
}
\author{Tony Chiang}

\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")
intactxml <- file.path(xmlDir, "intact_2008_test.xml")
x <- parsePsimi25Interaction(intactxml, INTACT.PSIMI25)
y <- interactions(x)[[1]]
getAbstractByPMID(pubmedID(y))
}
\keyword{ models }
