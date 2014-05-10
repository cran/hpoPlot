\name{remove.duplicate.term.sets}
\alias{remove.duplicate.term.sets}
\title{
	Remove Duplicate Term Sets
}
\description{
	Filters a list of character vectors of HPO terms so that there are no two equivalent.
}
\usage{
remove.duplicate.term.sets(hpo.terms, hpo.terms.set)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{hpo.terms.set}{
	List of character vectors of HPO terms.
}
}
\value{
	List of character vectors of HPO terms.
}
\examples{
	data(hpo)

	patients.character.vector <- c(
		Term.Set.A=
			"HP:0001382,HP:0004272,HP:0004912,HP:0007917",
		Term.Set.B=
			"HP:0001382",
		Duplicate.Of.Term.Set.A=
			"HP:0001382,HP:0004272,HP:0007917,HP:0004912"
	)

	patient.hpo.terms <- term.set.list.from.character(hpo.terms, patients.character.vector)

	no.dupes <- remove.duplicate.term.sets(hpo.terms, patient.hpo.terms)
}
\keyword{ duplicates }
