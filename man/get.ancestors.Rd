\name{get.ancestors}
\alias{get.ancestors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get ancestor terms of term set
}
\description{
Recursive function to find all the terms from which a term in the given set descends.
}
\usage{
get.ancestors(hpo.terms, target.terms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{target.terms}{
	Character vector of HPO terms.
}
}
\value{
Character vector of HPO terms which are ancestors to at least one of the `target.terms'
}

\examples{
	data(hpo)
	example.ancestors <- get.ancestors(hpo.terms, "HP:0002076")
	print(hpo.terms$name[example.ancestors])
}
\keyword{ ancestors }
