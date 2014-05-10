\name{get.term.adjacency.matrix}
\alias{get.term.adjacency.matrix}
\alias{get.term.pseudo.adjacency.matrix}
\title{
	Get HPO Term Adjacency Matrix
}
\description{
	Get adjacency matrix based on the HPO for a set of terms.
}
\usage{
get.term.adjacency.matrix(hpo.terms, terms)
get.term.pseudo.adjacency.matrix(hpo.terms, terms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hpo.terms}{
	An R-Object representing the HPO, as returned by \code{\link{get.hpo.terms}}.
}
  \item{terms}{
	Character vector of hpo terms.
}
}
\details{
get.term.adjacency.matrix only creates edges between parents and children, whereas get.term.pseudo.adjacency.matrix also creates edges from a node to its most recent ancestors in the set of terms if its parent is absent. 
}
\value{
	Matrix of logicals x, where x[row, col] is TRUE if row is a parent of column and FALSE otherwise.
}
\examples{
	data(hpo)

	migraine.terms <- get.ancestors(hpo.terms, hpo.terms$id[hpo.terms$name == "Migraine"])

	migraine.adj.mat <- get.term.adjacency.matrix(hpo.terms, migraine.terms)

	#plot simple graph of migraine terms
	plot(
		new(
			"graphAM", 
			adjMat=migraine.adj.mat, 
			edgemode="directed"
		)
	)

}
\keyword{ adjacency }
