\name{psimi25ComplexEntry-class}
\docType{class}
\alias{psimi25ComplexEntry-class}
\alias{complexes,psimi25ComplexEntry-method}
\alias{interactors,psimi25ComplexEntry-method}
\alias{show,psimi25ComplexEntry-method}

\title{Class "psimi25ComplexEntry"}
\description{ A class representing the complex data extracted from
  file in PSI-MI 2.5 format.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{parsePsimi25Complex}
}
\section{Slots}{
	 \describe{
    \item{\code{releaseDate}:}{Object of class \code{"character"},
      release date}
    \item{\code{interactors}:}{Object of class \code{"matrix"},
      interactors involved }
    \item{\code{complexes}:}{Object of class \code{"list"}, complexes}
  }
}
\section{Methods}{
  \describe{
    \item{complexes}{\code{signature(x = "psimi25ComplexEntry")}:
      returns all complexes }
    \item{interactors}{\code{signature(x = "psimi25ComplexEntry")}:
      returns information of interactors}
    \item{show}{\code{signature(object = "psimi25ComplexEntry")}: shows
      information of all the complexes}
	 }
}
\author{ Jitao David Zhang <jitao_david.zhang@roche.com>, Tony Chiang <tchiang@ebi.ac.uk>}
\seealso{
  \code{\link{parsePsimi25Complex}}
  \code{\linkS4class{psimi25ComplexEntry}}
}
\examples{
showClass("psimi25ComplexEntry")
}
\keyword{classes}
