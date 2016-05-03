\name{psimi25Interaction-class}
\docType{class}
\alias{psimi25Interaction-class}
\alias{interactionType,psimi25Interaction-method}
\alias{show,psimi25Interaction-method}
\alias{bait,psimi25Interaction-method}
\alias{prey,psimi25Interaction-method}
\alias{confidenceValue,psimi25Interaction-method}
\alias{inhibitor,psimi25Interaction-method}
\alias{neutralComponent,psimi25Interaction-method}
\alias{participant,psimi25Interaction-method}
\alias{pubmedID,psimi25Interaction-method}
\title{Class "psimi25Interaction"}
\description{ A class representing interaction data }
\section{Objects from the Class}{
Objects can be created by calls of the form \code{parsePsimi25Interaction}
}
\section{Slots}{
	 \describe{
    \item{\code{sourceDb}:}{Object of class \code{"character"}, source database}
    \item{\code{sourceId}:}{Object of class \code{"character"}, source
      database ID}
    \item{\code{interactionType}:}{Object of class \code{"character"},
      character, the method used for detecting the interaction ,such as
      "pull down"}
    \item{\code{expPubMed}:}{Object of class \code{"character"}, PubMed
      ID of the publication that describes the experiment}
    \item{\code{sourceId}:}{Object of class \code{"character"},
      source database ID of the experiment}
    \item{\code{confidenceValue}:}{Object of class \code{"character"},
      confidence value of the experimental interaction}
    \item{\code{participant}:}{Object of class \code{"character"},
      UniProt IDs of the participants. Important when no bait/prey
      information is available}
    \item{\code{bait}:}{Object of class \code{"character"}, UniProt ID
      of the bait }
    \item{\code{prey}:}{Object of class \code{"character"}, UniProt
      ID(s) of the prey(s)}
    \item{\code{inhibitor}:}{Object of class \code{"character"}, UniProt
    ID of the inhibitor, NA when missing}
    \item{\code{neutralComponent}:}{Object of class \code{"character"},
      UniProt ID of the neutral components, NA when missing}
    \item{\code{baitUniProt}:}{ An object of class \code{"character"}
     }
    \item{\code{preyUniProt}:}{ An object of class \code{"character"}
     }
 
   }
}
\section{Methods}{
  \describe{
    \item{interactionType}{\code{signature(x = "psimi25Interaction")}:
      finds out interaction type }
    \item{show}{\code{signature(object = "psimi25Interaction")}: a print
    method}
    \item{sourceDb}{\code{signature(x = "psimi25Interaction")}: finds
      out source database}
    \item{sourceId}{\code{signature(x = "psimi25Interaction")}: returns
      ID of the source database }
    \item{bait}{\code{signature(x = "psimi25Interaction")}: returns the
      UniProt ID of the bait}
    \item{prey}{\code{signature(x = "psimi25Interaction")}: returns the
      UniProt ID(s) of the prey(s)}
    \item{confidenceValue}{\code{signature(x = "psimi25Interaction")}:
    returns the confidence value of the interaction}
  }
}
\author{ Tony Chiang <tchiang@ebi.ac.uk> , Jitao David Zhang <jitao_david.zhang@roche.com> }
\seealso{
  \code{\link{parsePsimi25Interaction}},\code{S4CLASS\link{psimi25Interactor-class}},
  \code{S4CLASS\link{psimi25Experiment-class}}, \code{S4class\link{psimi25InteractionEntry-class}}
}
\examples{
showClass("psimi25Interaction")
}
\keyword{classes}
