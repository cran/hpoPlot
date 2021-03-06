% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sets.R
\name{get.descendants}
\alias{get.descendants}
\title{Get set of all descendants of single term}
\usage{
get.descendants(hpo.terms, ancestor, remove.ancestor = FALSE)
}
\arguments{
\item{hpo.terms}{R-Object representation of HPO}

\item{ancestor}{Character vector of length 1 - the HPO code of the term whose descendants you wish to retrieve}

\item{remove.ancestor}{Boolean indicating whether to remove the given ancestor or not}
}
\value{
Character vector of terms
}
\description{
Get set of all descendants of single term
}
\examples{
data(hpo.terms)
get.descendants(hpo.terms, ancestor=c("HP:0001873"))
}
\seealso{
\code{link{get.ancestors}}
}

