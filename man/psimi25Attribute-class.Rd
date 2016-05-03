\name{psimi25Attribute-class}
\Rdversion{1.1}
\docType{class}
\alias{psimi25Attribute-class}
\alias{iValue,psimi25Attribute-method}
\alias{iValue<-,psimi25Attribute-method}
\alias{name,psimi25Attribute-method}
\alias{name<-,psimi25Attribute,ANY-method}
\alias{name<-,psimi25Attribute,character-method}
\alias{nameAc,psimi25Attribute-method}
\alias{nameAc<-,psimi25Attribute-method}
\alias{show,psimi25Attribute-method}
\alias{show,psimi25AvailabilityType-method}
\alias{show,psimi25DbReferenceType-method}
\alias{show,psimi25NamesType-method}
\alias{value<-,psimi25Attribute-method}
\alias{xref,psimi25CommonNameRef-method}

\title{Class "psimi25Attribute" }
\description{ Persistence of the data structures specified by the PSI-MI
2.5 standard}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{new("psimi25Attribute", ...)}, or by the constructors.
}
\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"character"}: value}
    \item{\code{name}:}{Object of class \code{"character"}: name attribute}
    \item{\code{nameAc}:}{Object of class \code{"character"}: nameAc attribute }
  }
}

\section{Extends}{
Class \code{"\linkS4class{character}"}, from data part.
Class \code{"\linkS4class{vector}"}, by class "character", distance 2.
}
\section{Methods}{
  \describe{
    \item{iValue}{\code{signature(object = "psimi25Attribute")}: ... }
    \item{iValue<-}{\code{signature(object = "psimi25Attribute")}: ... }
    \item{name}{\code{signature(object = "psimi25Attribute")}: ... }
    \item{name<-}{\code{signature(object = "psimi25Attribute", value = "ANY")}: ... }
    \item{name<-}{\code{signature(object = "psimi25Attribute", value = "character")}: ... }
    \item{nameAc}{\code{signature(object = "psimi25Attribute")}: ... }
    \item{nameAc<-}{\code{signature(object = "psimi25Attribute")}: ... }
    \item{show}{\code{signature(object = "psimi25Attribute")}: ... }
    \item{value<-}{\code{signature(object = "psimi25Attribute")}: ... }
  }
}
\references{http://psidev.sourceforge.net/mi/rel25/doc/}
\author{ Jitao David Zhang }
\keyword{classes}
