\name{list2Matrix}
\alias{list2Matrix}
\title{ Converts list into matrix }
\description{
  Converts a named list representation of a graph into matrix
}
\usage{
list2Matrix(namedList, type="interaction")
}

\arguments{
  \item{namedList}{Named list.}
  \item{type}{Character: either interaction or complex}
}
\value{
  A matrix. A bait-prey matrix if type is "interaction" with baits indexing
  the columns and prey the rows or a protein complex incidence graph with
  complexes indexing the columns and proteins, the rows.

  If type is interaction, the names of the list are the baits and the entries
  of the list correspond the prey found by each bait.

  If type is complex, the names of the list are the names of the protein complex
  and the entries are the members of each corresponding complex.
}
\author{ Tony Chiang <tchiang@ebi.ac.uk> }
\keyword{ models }