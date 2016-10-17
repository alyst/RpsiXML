\name{buildPCHypergraph}
\alias{buildPCHypergraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Build protein complex hypergraph from PSI-MI 2.5 files}
\description{
  The protein complexes can be represented by hypergraph models, with
  proteins as nodes and complexes as hypergraphs. This function builds
  protein complex hypergraph from one or more PSI-MI 2.5 files (complex
  mode), with the option to split the dataset by organism name or
  taxonomy ID of the complexes.
}
\usage{
buildPCHypergraph(xmlFiles, psimi25source, split.by = c("none", "organismName", "taxId"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmlFiles}{ PSI-MI 2.5 files, must be of complex mode. }
  \item{psimi25source}{ PSI-MI 2.5 source indicator, for example
    INTACT.PSIMI25 for files from IntAct database}
  \item{split.by}{ The qualifier to split the dataset, \code{none}
    specifies not to split the dataset, \code{organismName} and
    \code{taxId} splits the dataset according to organism name or
    taxonomy ID respectively. In the case of splitting, the results are
    a list of hypergraphs indexed by respective split qualifier}
  \item{\dots}{ other parameters passed to \code{parsePsimi25Complex}}
}
\details{
  See \code{\link{psimi25Hypergraph-class}} for the use of resulting
  hypergraphs and examples.
}
\value{
  In case the dataset is not split (by setting the option
  \code{split.by} as \code{none}, the result is a
  \code{\link{psimi25Hypergraph-class}} object. If the dataset was
  split, a list of \code{\link{psimi25Hypergraph-class}} is returned,
  which is indexed by either the organism name or the taxonomy ID.
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{ \code{\link{psimi25Hypergraph-class}} for the use of
  hypergraph objects, \code{\link{separateXMLDataByExpt}} for similar
  functionality but for interaction mode files}
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")
intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")

pc2 <- buildPCHypergraph(intactComplexxml, INTACT.PSIMI25,split.by="organismName")

pc2[2]
complexes(pc2[[2]])[1:3]
}

