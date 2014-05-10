\name{get.hpo.terms}
\alias{get.hpo.terms}
\title{
Get HPO As R-Object
}
\description{
Parses a .obo file containing the HPO and returns an R-Object representing the HPO.
}
\usage{
get.hpo.terms( file )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
	The location of the .obo file (leave blank to download the current official version from "http://compbio.charite.de/hudson/job/hpo/lastStableBuild/artifact/ontology/release/hp.obo".
}
}
\value{
	\item{id}{Character vector of HPO term IDs}
	\item{name}{Character vector of HPO clinical descriptions names}
	\item{children}{Per-term list of character vectors of direct descedants}
	\item{parents}{Per-term list of character vectors of parent terms}
	\item{ancestors}{Per-term list of character vectors of ancestors}
	\item{siblings}{Per-term list of character vectors of sibling nodes}
	\item{alt.id}{Per-term list of character vectors of alternate/previous IDs}
	\item{date.downloaded}{Date that the object was constructed}
}
\references{
	Sebastian Kohler,Sandra C Doelken,Christopher J. Mungall,Sebastian Bauer,Helen V. Firth,Isabelle Bailleul-Forestier,Graeme C. M. Black,Danielle L. Brown,Michael Brudno,Jennifer Campbell,David R. FitzPatrick,Janan T. Eppig,Andrew P. Jackson,Kathleen Freson,Marta Girdea,Ingo Helbig,Jane A. Hurst,Johanna Jahn,Laird G. Jackson,Anne M. Kelly,David H. Ledbetter,Sahar Mansour,Christa L. Martin,Celia Moss,Andrew Mumford,Willem H. Ouwehand,Soo-Mi Park,Erin Rooney Riggs,Richard H. Scott,Sanjay Sisodiya,Steven Van Vooren,Ronald J. Wapner,Andrew O. M. Wilkie,Caroline F. Wright,Anneke T. Vulto-van Silfhout,Nicole de Leeuw,Bert B. A. de Vries,Nicole L. Washingthon,Cynthia L. Smith,Monte Westerfield,Paul Schofield,Barbara J. Ruef,Georgios V. Gkoutos,Melissa Haendel,Damian Smedley,Suzanna E. Lewis,and Peter N. Robinson - The Human Phenotype Ontology project: linking molecular biology and disease through phenotype data - Nucl. Acids Res. (1 January 2014) 42 (D1): D966-D974 doi:10.1093/nar/gkt1026
}

\examples{
	#to create ontology object from official .obo file...
	#hpo.terms <- get.hpo.terms()
}
\keyword{ HPO }
\keyword{ ontology }
