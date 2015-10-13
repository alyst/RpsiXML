\name{psimi25Graph-class}
\docType{class}
\alias{psimi25Graph-class}
\alias{show,psimi25Graph-method}
\alias{abstract,psimi25GraphBase-method}

\title{Class "psimi25Graph" ~~~ }
\description{ A graph object representing data extracted from PSI-MI 2.5
files}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{psimi25XML2Graph}}
}
\section{Slots}{
	 \describe{
    \item{\code{interactors}:}{Object of class \code{"matrix"},
      interactor information in a matrix, Each row represents one
      interactor. Source IDs are used as row names. Each column represents one annotation. 
      Annotations include: UniProt ID, short label, organism name, and NCBI taxonomy ID. Only those interactors which are 
      the nodes of the psimi25Graph is given. }
    \item{\code{abstract}:}{Object of class \code{pubMedAbst}}
 
  }
}
\section{Extends}{
Class \code{\link[graph:graphNEL-class]{graphNEL}}, directly.
Class \code{\link[graph:graphNEL-class]{graphNEL}}, by class "graphNEL", distance 2.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "psimi25Graph")}: show method }
    \item{translateSourceID}{\code{signature(r = "psimi25Graph")}:
      translate the source ID into other IDs }
    \item{abstract}{\code{signature(object="psimi25Graph")}: get the
      abstract information for the dataset from NCBI}
  }
}
\author{ Tony Chiang <tchiang@ebi.ac.uk> , Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{
  \code{\link{psimi25XML2Graph}},\code{S4class\link{psimi25Hypergraph-class}}
}
\examples{
showClass("psimi25Graph")
}
\keyword{classes}
