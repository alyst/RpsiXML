\name{null2na}
\alias{null2na}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Turns null or NA into character "NA"}
\description{
  The functions turns NULL or NA into character "NA"
}
\usage{
null2na(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A vector }
}
\value{
  Either the original vector (if not NULL or logical NA) or character NA
}
\author{ Tony Chiang <tchiang@ebi.ac.uk> }
\examples{
null2na(NA)
null2na(NULL)
}
\keyword{ models}