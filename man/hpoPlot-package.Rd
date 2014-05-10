\name{hpoPlot-package}
\alias{hpoPlot-package}
\alias{hpoPlot}
\docType{package}
\title{
Functions to plot graphviz style graphs of sets of HPO terms
}
\description{
Functions for visualising a set of HPO-encoded phenotypes in a variety of ways. 
The most important function in the package is hpo.graph, which expects to be given a list of HPO term character vectors (each typically the HPO encoded phenotype of a single patient), and some functions telling it how to display the terms in the resulting plot. 
}
\details{
\tabular{ll}{
Package: \tab hpoPlot\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-14\cr
License: \tab Public Domain\cr
}
}
\author{
Daniel Greene

Maintainer: Daniel Greene <dg333@cam.ac.uk>
}
\keyword{ package }
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

	#simple example
	hpo.graph(
		hpo.terms=hpo.terms,
		patients=patients.of.interest,
		main.title="Simple Example",
		filter.out.uninformative=TRUE
	)

	#complex example
	hpo.graph(
		hpo.terms=hpo.terms,
		patients=patients.of.interest,
		main.title="Complex Example",
		filter.out.uninformative=FALSE,
		colouring.function=colouring.functions$get.patient.based.colours,
		labelling.function=labelling.functions$get.patient.based.labels,
		size.function=size.functions$get.frequency.based.sizes,
		border.function=border.functions$get.no.borders
	)

}
