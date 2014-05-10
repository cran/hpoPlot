\name{term.set.list.from.character}
\alias{term.set.list.from.character}
\title{
Get Term Set List From Character Vector
}
\description{
Function which returns a list of HPO term sets (typically to represent patients and used by functions in this package) given a character vector of comma separated HPO terms.
}
\usage{
term.set.list.from.character(hpo.terms, character.vector, remove.empty = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{character.vector}{
	A character vector of comma separated lists of HPO terms
}
  \item{remove.empty}{
	Boolean to dictate whether to remove sets with no HPO terms.
}
}
\value{
	List of character vectors of HPO terms.
}
\examples{
	data(hpo)

	patients.character.vector <- c(
		A="HP:0001382,HP:0004272,HP:0007917,HP:0004912,HP:0001596",
		B="HP:0001382,HP:0004272,HP:0002165,HP:0004800,HP:0004912",
		C="HP:0004800,HP:0001382,HP:0004912,HP:0007917,HP:0008743",
		D="HP:0001257,HP:0001382,HP:0007917,HP:0012623,HP:0002165",
		E="HP:0007917,HP:0004800,HP:0004272,HP:0001596,HP:0002165" 
	)

	patients.of.interest <- term.set.list.from.character(hpo.terms, patients.character.vector)

}
