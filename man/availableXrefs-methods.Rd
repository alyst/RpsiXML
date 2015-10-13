\name{availableXrefs-methods}
\docType{methods}
\alias{availableXrefs-methods}
\alias{availableXrefs,list-method}
\alias{availableXrefs,psimi25InteractionEntry-method}
\alias{availableXrefs,psimi25Interactor-method}
\title{ Methods for function availableXrefs}
\description{
 Methods for function \code{availableXrefs}. See
 \code{\link{availableXrefs}} for more details
}
\section{Methods}{
\describe{

\item{x = "list"}{ A list of \code{\link{psimi25Interactor-class}} objects, an
logical option "intersect" can be set to tell the function whether to
return the intersected or the unioned cross ref databases.}

\item{x = "psimi25InteractionEntry"}{ Similar to that of "list",
  "intersect" option is also provided}

\item{x = "psimi25Interactor"}{ returns the avaiable cross refences of
  the given interactor, "intersect" option not available }

}}
\keyword{methods}
\keyword{models}
