\name{get.descendants}
\alias{get.descendants}
\title{
Get Descendants
}
\description{
Get all terms which descend from a given term in the HPO.
}
\usage{
get.descendants(hpo.terms, ancestor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{ancestor}{
	Single HPO term
}
}
\value{
Character vector of HPO terms which descend from `ancestor'.
}

\examples{
	data(hpo)

	types.of.headache <- get.descendants(
		hpo.terms,
		hpo.terms$id[hpo.terms$name == "Headache"]
	)

	print(hpo.terms$name[types.of.headache])
}
\keyword{ descendants }
