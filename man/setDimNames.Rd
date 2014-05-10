\name{setDimNames}
\alias{setDimNames}
\title{
Set Names For Dimensions Of Array
}
\description{
A higher dimensional version of \code{\link{setNames}}.
}
\usage{
setDimNames(array.object, list.of.dimension.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{array.object}{
	Array to be named
}
  \item{list.of.dimension.names}{
	A list of character vectors whose lengths are the same as dim(array.object).
}
}
\value{
Named array
}
\examples{
	setDimNames(	
		matrix(1,2,2),
		list(
			c("Row1","Row2"),
			c("Col1","Col2")
		)
	)
}
\keyword{ arrays }
