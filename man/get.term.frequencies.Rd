% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\name{get.term.frequencies}
\alias{get.term.frequencies}
\title{Get frequency of each term in a set of phenotypes}
\usage{
get.term.frequencies(hpo.terms, hpo.phenotypes, patch.missing = FALSE)
}
\arguments{
\item{hpo.terms}{R-Object representation of HPO}

\item{hpo.phenotypes}{List of HPO term character vectors}

\item{patch.missing}{Logical indicating whether to include all HPO even if they're not present in the \code{hpo.phenotypes} as if they had occurred once}
}
\value{
Numeric vector of information contents, named by corresponding HPO terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
}
\description{
Get frequency of each term in a set of phenotypes
}
\examples{
data(hpo.terms)
get.term.frequencies(hpo.terms, list("HP:0001873"))
}
\seealso{
\code{\link{get.term.info.content}}
}

