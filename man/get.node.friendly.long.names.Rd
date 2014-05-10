\name{get.node.friendly.long.names}
\alias{get.node.friendly.long.names}
\title{
Get Node-Friendly Long Names
}
\description{
Get the clinical names of HPO terms, shortened, where possible, to exclude phrases like `Abnormality of' and split onto several lines where necessary so they fit neatly onto nodes
}
\usage{
get.node.friendly.long.names(hpo.terms, terms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{terms}{
	Character vector of HPO terms.
}
}
\value{
	Character vector of descriptions
}
\examples{
	data(hpo)

	body.height <- hpo.terms$id[hpo.terms$name == "Abnormality of body height"]

	print(
		paste(
			hpo.terms$name[hpo.terms$id == body.height],
			" -> ",
			get.node.friendly.long.names(
				hpo.terms,
				body.height
			),
			sep=""
		)
	)
}
