\name{swap.out.alt.ids}
\alias{swap.out.alt.ids}
\title{
Swap Out Alternative/Previous IDs
}
\description{
Function to replace the terms in a set with their updated/replacement IDs (if they exist) in the given HPO object.
}
\usage{
swap.out.alt.ids(hpo.terms, hpo.term.set, remove.dead = FALSE)
}
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{hpo.term.set}{
	Character vector of HPO terms.
}
  \item{remove.dead}{
	Boolean to dictate whether to remove terms which are not in the HPO object given. Defaults to FALSE.
}
}
\value{
Character vector of HPO terms.
}
\examples{
	data(hpo)
	old.code.for.migraine <- "HP:0007194"
	current.code.for.migraine <- swap.out.alt.ids(hpo.terms, old.code.for.migraine)
}
