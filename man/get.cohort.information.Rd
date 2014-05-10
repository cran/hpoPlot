\name{get.cohort.information}
\alias{get.cohort.information}
\title{
Get Information Content Of Terms
}
\description{
Get information content (i.e. - log frequency) of term with respect to a given cohort of patients/term-sets.
}
\usage{
get.cohort.information(hpo.terms, patient.hpo.terms, patch.missing = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{patient.hpo.terms}{
	List of HPO term sets.
}
  \item{patch.missing}{
	Return an information content for all HPO terms (i.e. including those not present in the cohort) and set it to the highest information content.
}
}
\value{
Numeric vector of information contents, named by term.
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

	patient.hpo.terms <- term.set.list.from.character(hpo.terms, patients.character.vector)

	cohort.information <- get.cohort.information(hpo.terms, patient.hpo.terms)
}
\keyword{ information }
\keyword{ frequency }
