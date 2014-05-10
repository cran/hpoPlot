\name{get.term.patient.matrix}
\alias{get.term.patient.matrix}
\title{
Get Term By Patient Matrix
}
\description{
Create a matrix of logicals where the patients are represented by the rows, and the HPO terms present in at least 1 patient are represented by the columns. The values in the cells are then TRUE if the patient of the row has the phenotype of the column and FALSE otherwise.
}
\usage{
get.term.patient.matrix(patient.hpo.terms.with.ancs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{patient.hpo.terms.with.ancs}{
	List of HPO term character vectors, most usually containing ancestors. Such a list can be obtained from a list containing only the leaf terms using `lapply(list.of.term.sets, function(x) get.ancestors(hpo.terms, x))'.
}
}
\value{
Logical matrix 
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
	
	term.patient.matrix <- get.term.patient.matrix(
		lapply(
			patient.hpo.terms, 
			function(patient.terms) get.ancestors(hpo.terms, patient.terms)
		)
	)
}
