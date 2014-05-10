\name{clean.terms}
\alias{clean.terms}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Leaves Of HPO Term Set
}
\description{
Get a set of leaf nodes which represent a given set of terms.
}
\usage{
clean.terms(hpo.terms, terms)
}
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{terms}{
	Character vector of terms to `clean'.
}
}
\value{
Character vector of leaf terms
}

\examples{
	data(hpo)
	headache.terms <- hpo.terms$id[
		hpo.terms$name %in% c("Headache", "Migraine")
	]
	leaf.terms <- clean.terms(hpo.terms, headache.terms)
}
\keyword{ redundancy }
