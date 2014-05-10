\name{suggest.aspect.ratio.width.to.height}
\alias{suggest.aspect.ratio.width.to.height}
\title{
Suggest Aspect Ratio For Plot
}
\description{
Suggests a width to height aspect ratio for an HPO term plot based on the number of terms which are to be plotted.
}
\usage{
suggest.aspect.ratio.width.to.height(hpo.terms, list.of.patients.or.vector.of.terms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{list.of.patients.or.vector.of.terms}{
	List of vectors of terms or vector of terms which are to be plotted.
}
}
\value{
	Scalar value of suggested aspect ratio.
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

	height=10 

#	pdf(
#		"test.pdf", 
#		height=height,
#		width=height * suggest.aspect.ratio.width.to.height(
#			hpo.terms,
#			patients.of.interest
#		)
#	)
#
#	hpo.graph(
#		hpo.terms=hpo.terms,
#		patients=patients.of.interest,
#		main.title="Simple Example",
#		filter.out.uninformative=TRUE
#	)
#
#	dev.off()

}
\keyword{ aspect-ratio }
