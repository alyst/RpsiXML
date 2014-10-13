\name{psimi25XML2Graph}
\alias{psimi25XML2Graph}
\alias{interactionEntry2graph}
\alias{complexEntry2graph}
\title{ Convert a vector of PSI-MI 2.5 XML files into graph objects }
\description{
  The function \code{psimi25XML2Graph} take a vector of XML 2.5 files
  from te same data source and generates a graph object based on the
  type of the files.
  \code{psimi25XML2Graph} is a wrapper for \code{interactionEntry2graph}
  and \code{complexEntry2graph}, which transform interactionEntry list
  and complexEntry list into graphs respectively.
}
\usage{
psimi25XML2Graph(psimi25files, psimi25source, type = "interaction",
directed = TRUE, ...)
}

\arguments{
  \item{psimi25files}{ Single file name or a vector of PSI-MI 2.5 XML file names or
    URLs. In case of splitted data the latter form is preferred. Different datasets or datasets from different sources should
    not be put into the same vector. }
  \item{psimi25source}{ Source of the PSI-MI 2.5 XML file, see \code{\link{psimi25Source-class}} }
  \item{type}{ A character string which is either "interaction" or
    "complex". As the value suggests, use "interaction" if the XML
    file contains experimental physical data, and "complex" if the file
    contains curated protein complex membership data.}
  \item{directed}{ Logical, whehter the returned graph object should be
    directed or undirected }
  \item{...}{ Other parameters passed to
    \code{\link{parsePsimi25Interaction}} }
}
\value{
 If type is "interaction",  then a resulting psimi25Graph object is genertated on the
 aggregation of the XML files. Otherwise if type is "complex,", a resulting psimi25HyperGraph
 object is generated on the aggregate of the XML files.
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
intactGraph <- psimi25XML2Graph(intactxml, INTACT.PSIMI25, type="interaction")

intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
intactComplexGraph <- psimi25XML2Graph(intactComplexxml, INTACT.PSIMI25, type="complex")
}
\keyword{ models }
