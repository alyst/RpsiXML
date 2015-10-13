\name{translateID-methods}
\docType{methods}
\alias{translateID-methods}
\alias{translateID,psimi25Graph-method}
\alias{translateID,psimi25Hypergraph-method}
\alias{translateID,psimi25Interactor-method}
\alias{translateID,list-method}
\title{ Translate interactors into other identifiers }
\description{
  The method finds any annotation IDs of of
  \code{\link{psimi25Graph-class}} or
  \code{\link{psimi25Hypergraph-class}} or
  \code{\link{psimi25Interactor-class}}
}
\section{Methods}{
\describe{

\item{r = "psimi25Graph"}{ An object of
  \code{\link{psimi25Graph-class}}. It replaces the nodes with
  translated IDs and returns a new psimi25Graph object back.}

\item{r = "psimi25Hypergraph"}{ An object of
  \code{\link{psimi25Hypergraph-class}}. Since the nodes of
  psimi25Hypergraph are not mutable, this method accepts a UniProt ID of
  interactors and returns the translated ID back. NA is also acceptable,
  which can be used to check interactors without mapping to UniProt}

\item{r = "list" }{ A list of \code{\link{psimi25Interactor-class}},
  returns a vector of translated identifiers back}

\item{r = "psimi25Interactor"}{An object of
  \code{\link{psimi25Interactor-class}}, returns translated
  identifier. Before using the method it is advisable to used xref(x)
  method to check which cross references are provided.}
}}
\keyword{methods}
\keyword{models}
