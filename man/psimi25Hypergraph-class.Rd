\name{psimi25Hypergraph-class}
\docType{class}
\alias{psimi25Hypergraph-class}
\alias{hyperedgeNodes}
\alias{edgeLabel}
\alias{revInciMat}
\alias{initialize,Hypergraph-method}
\alias{initialize,psimi25Hypergraph-method}
\alias{show,psimi25Hypergraph-method}
\alias{complexes,psimi25Hypergraph-method}
\alias{edgeLabel,psimi25Hypergraph-method}
\alias{hyperedgeNodes,psimi25Hypergraph-method}
\alias{organismName,psimi25Hypergraph-method}
\alias{numEdges,psimi25Hypergraph-method}
\alias{interactorInfo,psimi25Hypergraph-method}
\alias{hyperedgeNodes,Hypergraph-method}
\alias{revInciMat,matrix-method}

\title{Class "psimi25Hypergraph"}
\description{ Class to present PSI-MI 2.5 XML data as
  hypergraph. Proteins are projected as hypergraph nodes and complex
  composition as hyperedges.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{psimi25XML2Graph}
}
\section{Slots}{
  \describe{
    \item{\code{interactors}:}{Object of class \code{"matrix"}, Object of class \code{"matrix"},
      interactor information in a matrix, Each row represents one
      interactor. Source IDs are used as row names. Each column represents one annotation. 
      Annotations include: UniProt ID, short label, organism name, and NCBI taxonomy ID. Only those interactors which are 
      the nodes of the psimi25Graph is given.}
  }
}
\section{Extends}{
Class \code{Hypergraph}, directly
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "psimi25Hypergraph")} }
    \item{show}{\code{signature(object = "psimi25Hypergraph")}: show
      method,print complex and protein number}
    \item{interactors}{\code{signature(object = "psimi25Hypergraph")}:
      list of \code{psimi25Interactor} objects, providing full
      information of complex members}
    \item{edgeLabel}{\code{signature(object = "psimi25Hypergraph")}:
      returns complex names (as hyperedge label), as a character vector}
    \item{hyperedgeNodes}{\code{signature(object = "Hypergraph")}:
      returns a list of characters: names of the list are complex names and character vector in each
      list item are the members of that complex}
    \item{complexes}{\code{signature(object = "psimi25Hypergraph")}:
      a wrapper of \code{hyperedgeNodes}}
    \item{translateSourceID}{\code{signature(r = "psimi25Hypergraph")}:
      translate source ID into other IDs }
    \item{numInteractors}{\code{signature(r = "psimi25Hypergraph")}:
      returns the number of proteins}
    \item{interactorInfo}{\code{signature(r = "psimi25Hypergraph")}:
      returns a data frame containing essential information of the interactors}
    \item{numEdges}{\code{signature(r = "psimi25Hypergraph")}:
      returns the number of complexes}
    \item{revInciMat}{\code{signature(r = "matrix")}: returns the
      hypergraph built from the incidence matrix}
    
  }
}
\author{ Tony Chiang <tchiang@ebi.ac.uk> , Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{
  \code{\link{psimi25XML2Graph}},\code{S4class\link{psimi25Graph-class}}
}
\examples{
xmlDir <- system.file("/extdata/psi25files",package="RpsiXML")
intactComplexxml <- file.path(xmlDir,"intact_complexSample.xml")
pc1 <- buildPCHypergraph(intactComplexxml, INTACT.PSIMI25)

## print number of proteins and complexes (edges)
numNodes(pc1)
## the same as numInteractors(pc1)
numEdges(pc1)

## print proteins (nodes)
nodes(pc1)[1:3]

## print complex names
edgeLabel(pc1)[1:3]

## print complexes (not so informative with 'hyperedges')
hyperedges(pc1)[1:3]
## better with 'complexes' or 'hyperedgeNodes'
complexes(pc1)[1:3]

## get interactor detailed information
interactors(pc1)[[1]]
}
\keyword{classes}
