\name{remove.uninformative.terms}
\alias{remove.uninformative.terms}
\title{
Remove Uninformative Terms From Set
}
\description{
Removes terms which don't add information about which terms people represented in the plot are annotated with.
}
\usage{
remove.uninformative.terms(hpo.terms, patient.hpo.terms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{patient.hpo.terms}{
	A named list of HPO term character vectors representing the phenotypes of a set of patients.
}
}
\value{
	A character vector of the HPO terms which are informative in the context of the patients' phenotypes, in the sense that the frequency can't be deduced from the frequencies of its child terms alone.
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

	patients <- term.set.list.from.character(hpo.terms, patients.character.vector)

	informative.terms <- remove.uninformative.terms(hpo.terms, patients)
}
\keyword{ uninformative }

