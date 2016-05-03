\name{separateXMLDataByExpt}
\alias{separateXMLDataByExpt}
\title{Convert a vector of PSI-MI 2.5 XML files into graph objects based on pubmedID}
\description{
  The function \code{psimi25XML2Graph} take a vector of XML 2.5 files
  from te same data source and generates a graph object based on the
  type of the files.
}
\usage{
separateXMLDataByExpt(xmlFiles, psimi25source, type = "direct", directed = TRUE, abstract = FALSE,...)
}
\arguments{
  \item{xmlFiles}{ Single file name or a vector of PSI-MI 2.5 XML file names or
    URLs.}
  \item{type}{A character. Currently the user can specify to cull either "direct"
    interactions or "indirect" interactions.}
  \item{psimi25source}{ Source of the PSI-MI 2.5 XML file, see \code{\link{psimi25Source-class}}}
  \item{directed}{ Logical, whehter the returned graph object should be
    directed or undirected.}
  \item{abstract}{Logical; if TRUE, the abstract information will be
    appended to the graph object.}
  \item{...}{ Other parameters passed to
    \code{\link{parsePsimi25Interaction}}, for example \code{verbose=TRUE}}
}
\value{
 A list of \code{\link{psimi25Graph-class}} are generated indexed by the pubmedID of
 each bait-prey interaction. WARNING - the abstract information is obtained using the 
 pubmed and buildPubMedAbst functions from the 
 annotate package which warns the user that NCBI may block access to
 their site. The default is to not obtain the abstract for this reason. 
}
\author{ Jitao David Zhang, Tony Chiang }
\seealso{
  \code{\link{psimi25Source-class}},
  \code{\link{psimi25Graph-class}},
  \code{\link{psimi25Hypergraph-class}}
}
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")

intactxml <- file.path(xmlDir, "intact_2008_test.xml")
intactGraph <- separateXMLDataByExpt(intactxml, INTACT.PSIMI25, type="indirect")
}
\keyword{ models }
