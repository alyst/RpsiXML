\name{psimi25Interactor-class}
\docType{class}
\alias{psimi25Interactor-class}
\alias{show,psimi25Interactor-method}
\alias{xref,psimi25Interactor-method}

\title{A class representing interactor from PSI-MI 2.5 XML files}
\description{A class representing interactor from PSI-MI 2.5 XML files}
\section{Slots}{
	 \describe{
    \item{\code{sourceDb}:}{Object of class \code{"character"}, source database}
    \item{\code{sourceId}:}{Object of class \code{"character"}, ID of
      the source database }
    \item{\code{shortLabel}:}{Object of class \code{"character"}, short
      label of the interactor, if not available the slot will be filled
      with "fullName" node of the file}
    \item{\code{uniprotId}:}{Object of class \code{"character"}, UniProt
    ID of the interactor}
    \item{\code{organismName}:}{Object of class \code{"character"}, the
      unique name of the organism name }
    \item{\code{taxId}:}{Object of class \code{"character"}, the NCBI
      taxonomy ID}
    \item{\code{xref}:}{ Object of class \code{"environment"}
      }
  }
}
\section{Methods}{
 \describe{
  \item{show}{\code{signature(x = "psimi25Interactor")}: a print method}
  \item{sourceDb}{\code{signature(x = "psimi25Interactor")}: finds
    out source database}
  \item{sourceId}{\code{signature(x = "psimi25Interactor")}: returns
    ID of the source database }
  \item{xref}{\code{sigature(x = "psimi25Interactor")}: returns
    cross-references of the given interactor}
 }
}
\author{ Jitao David Zhang \email{jitao_david.zhang@roche.com} }

%\seealso{
%  \code{\link{parseInteractor}}
%}
\examples{
showClass("psimi25Interactor")
}
\keyword{classes}
