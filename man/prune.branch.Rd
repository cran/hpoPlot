% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sets.R
\name{prune.branch}
\alias{prune.branch}
\title{Prune all terms descending from given term down to that term and ensure no degeneracy}
\usage{
prune.branch(hpo.terms, prune.to.point, terms)
}
\arguments{
\item{hpo.terms}{R-Object representation of HPO}

\item{prune.to.point}{HPO term which can be included, but whose descendants should be excluded}

\item{terms}{Character vector of HPO terms}
}
\value{
Character vector of terms
}
\description{
Prune all terms descending from given term down to that term and ensure no degeneracy
}

