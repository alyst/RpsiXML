\name{psimi25Source-class}
\docType{class}
\alias{psimi25Source-class}
\alias{uniprot,psimi25Source-method}

\title{Class "psimi25Source"}
\description{ A class representing data source compatible with PSI-MI
  2.5 XML data exchange format}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("psimi25Source", ...)}
}
\section{Slots}{
	 \describe{
    \item{\code{label}:}{Object of class \code{"character"}, a
      meaningful unique label of the database, does not necessarily have to be
      the identifier used in the PSI-MI files. The value is advised to
      be overwritten by the objects}
    \item{\code{sourceDb}:}{Object of class \code{"character"}, the
      character string used in the PSI-MI files to identify the source
      database. The value must be overwritten by the object, depending
      on the database identifier used in the PSI-MI file.}
    \item{\code{uniprotSymbol}:}{Object of class \code{"character"}, the
    identifier of UniProf reference in "primaryRef" or "secondaryRef"
    attributes of the interactors, telling the parser where to find
    UniProt symbol of the interactors. The value must be overwritten by
    the objects.}
  }
}
\section{Methods}{
  \describe{
    \item{sourceDb}{\code{signature(x = "psimi25Source")}: returns the
      source database}
    \item{uniprot}{\code{signature(x = "psimi25Source")}: returns the
      uniprot identifier of the PSI-MI 2.5 XML file }
  } 
}
       
\author{ Jitao David Zhang <jitao_david.zhang@roche.com> }
\note{
  In case of a new data repository, you can first try creating a
  new object of the class with \code{label}, \code{sourceDb} and
  \code{uniprotSymbol} adjusted as needed. Try parsing a function with
  this new psimi25Source object. In case of warnings or errors, it may
  mean that the file provided is not strictly compatible with PSI-MI 2.5
  data exchange format standard, please then contact the data provider.
}

\seealso{
  \code{\link{parsePsimi25Interaction}},
  \code{\link{parsePsimi25Complex}},
  \code{\link{psimi25XML2Graph}}
}
\examples{
showClass("psimi25Source")
}
\keyword{classes}
