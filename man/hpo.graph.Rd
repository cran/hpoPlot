\name{hpo.graph}
\alias{hpo.graph}
\title{
Plot HPO Terms
}
\description{
Function to display, or save to a PDF, a plot of sets of HPO terms.
}
\usage{
hpo.graph(
	hpo.terms, 
	patients, 
	term.information.content = NULL, 
	term.population.frequency = NULL, 
	terms=NULL, 
	colouring.function = colouring.functions$get.frequency.based.colours, 
	labelling.function = labelling.functions$get.frequency.based.labels, 
	border.function = border.functions$get.no.borders, 
	size.function = size.functions$get.standard.sizes, 
	main.title = NULL, 
	draw.legend = FALSE, 
	filter.out.uninformative = FALSE, 
	filter.out.non.pa.terms = FALSE, 
	min.occs = 1, 
	nodeAttrs = NULL, 
	pdf.file.name = NULL, 
	pdf.height = 10, 
	pdf.width = NULL
)
}

\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{patients}{
	A list of character vectors of HPO terms, where typically each list item represents a patient. e.g. patients=list(A="HP:00000118")
}
  \item{term.information.content}{
	The `Information-Content' (i.e. -log frequency) of the terms relevant to the plot as a named (by the HPO ID) vector in the context of a particular cohort of patients (as returned by \code{\link{get.cohort.information}}).
}
  \item{term.population.frequency}{
	An alternative to the `term.information.content', the frequency of the terms relevant to the plot in the population, as opposed to in the sample to be plotted, as a named vector of values between 0 and 1. If both this and `term.information.content' are given, `term.information.content' is used.
}
  \item{terms}{
	Character vector of HPO term codes for terms which are to appear in the plot. Defaults to NULL whereafter the terms in the plot are chosen automatically depending on the filtering options (for example `filter.out.uninformative').
}
  \item{colouring.function}{
	A named vector of colours for the terms in the plot or a function which returns one (see \code{\link{colouring.functions}}).
}
  \item{labelling.function}{
	A named vector of labels for the terms in the plot or a function which returns one (see \code{\link{labelling.functions}}).
}
  \item{border.function}{
	A named vector of border colours for the terms in the plot or a function which returns one (see \code{\link{border.functions}}).
}
  \item{size.function}{
	A named vector of node radii (in inches) for the terms in the plot or a function which returns one (see \code{\link{size.functions}}).
}
  \item{main.title}{
	Title for the plot
}
  \item{draw.legend}{
	Boolean to decide whether legend is drawn. NOTE: The legend drawn may intersect the plot.
}
  \item{filter.out.uninformative}{
	Filter out uninformative nodes (i.e. nodes add no new information to the plot, given the presence of their children).
}
  \item{filter.out.non.pa.terms}{
	Filter out nodes which are not descendants of the `Phenotypic abnormality' HPO-term.
}
  \item{min.occs}{
	Filter out nodes which occur less than `min.occs' times in `patients'.
}
  \item{nodeAttrs}{
	A custom Rgraphviz node attributes object which overrides any other graphical parameters passed to this function. Use this to set font size, font colour and more (for a full list see the Rgraphviz documentation).
}
  \item{pdf.file.name}{
	File name. If null will plot to the current graphical device.
}
  \item{pdf.height}{
	Height in inches of PDF created.
}
  \item{pdf.width}{
	Width in inches of PDF created (if left NULL, the function will try to set this based on the pdf.height parameter and the number of leaves in the plot, so that the graph looks legible).
}
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

	hpo.graph(
		hpo.terms=hpo.terms,
		patients=patients.of.interest,
		main.title="Test Plot",
		filter.out.uninformative=TRUE
	)

}
