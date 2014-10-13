\name{psimi25InteractionEntry-class}
\docType{class}
\alias{psimi25InteractionEntry-class}
\alias{interactions,psimi25InteractionEntry-method}
\alias{interactors,psimi25InteractionEntry-method}
\alias{organismName,psimi25InteractionEntry-method}
\alias{releaseDate}
\alias{releaseDate-methods}
\alias{releaseDate,psimi25Entry-method}
\alias{show,psimi25InteractionEntry-method}
\alias{taxId,psimi25InteractionEntry-method}
\alias{pubmedID,psimi25InteractionEntry-method}
\title{Class "psimi25InteractionEntry" }
\description{ A class representing interaction data parsed from PSI-MI
  2.5 XML files}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{parsePsimi25Interaction}
}
\section{Slots}{
  \describe{
    \item{\code{organismName}:}{Object of class \code{"character"}, the unique organism name(s) of the interactors }
    \item{\code{taxId}:}{Object of class \code{"character"}, the unique
      NCBI taxonomy ID(s) of the interactors }
    \item{\code{releaseDate}:}{Object of class \code{"character"}, character, release date of the data entry, 
	recorded in the \code{entry} element of the XML file}
    \item{\code{interactors}:}{Object of class \code{"list"}, a list of
      \code{psimi25Interactor-class} objects, each represents one interactor.}
    \item{\code{interactions}:}{Object of class \code{"list"}, list, a list of \code{psimi25Interaction-class} objects, 
	each represents one interaction. }
  }
}
\section{Methods}{
  \describe{
    \item{interactions}{\code{signature(x = "psimi25InteractionEntry")}:
    return a list of \code{psimi25Intearaction} objects, each
    representing one interaction}
    \item{interactors}{\code{signature(x = "psimi25InteractionEntry")}: a method to find 
      all the interactors within the psimi25InteractionEntry and print them }
    \item{organismName}{\code{signature(x = "psimi25InteractionEntry")}:
      returns organism names of the interactors }
    \item{releaseDate}{\code{signature(x = "psimi25InteractionEntry")}:
      returns the release date, serving as a proxy for versioning }
    \item{show}{\code{signature(object = "psimi25InteractionEntry")}: a
      print method}
    \item{taxId}{\code{signature(x = "psimi25InteractionEntry")}:
      returns NCBI taxonomy ID }
    \item{numInteractors}{\code{signature(x = "psimi25InteractionEntry")}:
      returns the number of interactors.}
    \item{numInteractions}{\code{signature(x = "psimi25InteractionEntry")}:
      returns the number of interactions.}
    \item{pubmedID}{\code{signature(x = "psimi25InteractionEntry")}:
      returns the (unique) PubMed IDs of the papers reporting the
      interactions in the entry.}
  }
  
}
\author{ Tony Chiang <tchiang@ebi.ac.uk> , Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{
  \code{\link{parsePsimi25Interaction}}, \code{S4class\link{psimi25Interaction-class}}
}
\examples{
showClass("psimi25InteractionEntry")
}
\keyword{classes}
