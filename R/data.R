#' Get human readable, shortened (where possible) HPO term names
#'
#' @template hpo.terms
#' @template terms
#' @return Character vector
#' @examples
#' data(hpo.terms)
#' get.shortened.names(hpo.terms, c("HP:0001873", "HP:0011877"))
#' @export
#' @import magrittr
get.shortened.names <- function(hpo.terms, terms) gsub(
	"Impaired |(Abnormality of (the )?)|(Abnormal )", 
	"", 
	hpo.terms$name[terms]
) %>% 
(function (x) gsub("^\\s+|\\s+$", "", x)) %>%
sapply(simpleCap)

#' Get R-Object representation of ontology from obo file
#'
#' @param file File path of obo file
#' @param species Character vector - "H" for HPO, "M" for MPO
#' @return R-Object (list) representing ontology
#' @export
#' @import magrittr
get.ontology <- function(file, species="H") {
	hpo.obo.lines <- readLines(file)

	hpo.term.id.pattern <- paste("^id: (", species, "P:\\d+)", sep="")
	hpo.term.name.pattern <- paste("^name: (\\.*)", sep="")
	hpo.term.parent.pattern <- paste("^is_a: (", species, "P:\\d+)", sep="")
	hpo.term.pattern <- paste("", species, "P:\\d+", sep="")
	hpo.term.alt.id.pattern <- paste("^alt_id: (", species, "P:\\d+)", sep="")

	hpo.term.id.lines <- grep(hpo.term.id.pattern, hpo.obo.lines)
	hpo.term.name.lines <- grep(hpo.term.name.pattern, hpo.obo.lines)
	hpo.term.parent.lines <- grep(hpo.term.parent.pattern, hpo.obo.lines)
	hpo.term.alt.id.lines <- grep(hpo.term.alt.id.pattern, hpo.obo.lines)

	hpo.terms <- NULL

	hpo.terms$id <- sub(
		hpo.term.id.pattern,
		"\\1",
		hpo.obo.lines[hpo.term.id.lines]
	)

	hpo.terms$name <- sub(
		hpo.term.name.pattern, 
		"\\1", 
		hpo.obo.lines[hpo.term.name.lines]
	)

	names(hpo.terms$name) <- hpo.terms$id

	Encoding(hpo.terms$name) <- "latin1"
	hpo.terms$name <- iconv(
		hpo.terms$name,
		"latin1",
		"ASCII",
		sub=""
	)

	hpo.parent.term.matches <- regexpr(
		hpo.term.pattern,
		hpo.obo.lines[hpo.term.parent.lines]
	)

	hpo.parent.terms <- substr(
		hpo.obo.lines[hpo.term.parent.lines],
		hpo.parent.term.matches,
		hpo.parent.term.matches + attr(hpo.parent.term.matches, "match.length")-1
	)

	hpo.terms$parents <- split(
		hpo.parent.terms,
		cut(
			hpo.term.parent.lines,
			breaks=c(hpo.term.id.lines,length(hpo.obo.lines)+1),
			labels=hpo.terms$id
		)
	)
	
	hpo.alt.id.matches <- regexec(
		hpo.term.pattern,
		hpo.obo.lines[hpo.term.alt.id.lines]
	)
	
	hpo.alt.ids <- regmatches(
		hpo.obo.lines[hpo.term.alt.id.lines],
		hpo.alt.id.matches
	)

	hpo.alt.ids <- sapply(hpo.alt.ids, function(x) x[1])
	
	hpo.terms$alt.id <- as.character( 
		cut(
			hpo.term.alt.id.lines,
			breaks=c(hpo.term.id.lines,length(hpo.obo.lines)+1),
			labels=hpo.terms$id
		)
	)	
	names(hpo.terms$alt.id) <- hpo.alt.ids

	get.ancestors <- function(terms) {
		result <- unique(terms)
		for (term in terms)
			for (parent in hpo.terms$parents[[term]])
				result <- union(result, get.ancestors(parent))
		result
	}

	names(hpo.terms$id) <- hpo.terms$id

	hpo.terms$children <- lapply(hpo.terms$id, function(x) c())
	for (hpo.term in hpo.terms$id)
		for (parent.term in hpo.terms$parents[[hpo.term]])
			hpo.terms$children[[parent.term]] <- c(
				hpo.terms$children[[parent.term]],
				hpo.term
			)

	hpo.terms$ancestors <- lapply(
		hpo.terms$id,
		get.ancestors
	)
 	
	hpo.terms$siblings <- lapply(
		hpo.terms$id,
		function(x) unlist(hpo.terms$children[hpo.terms$parents[[x]]])
	)

	hpo.terms$date.downloaded <- Sys.Date()

	hpo.terms
}

#' Gets HPO R-Object
#'
#' @param file File path of obo file
#' @return R-Object (list) representing ontology
#' @export
#' @import magrittr
get.hpo.terms <- function(file) get.ontology(
	file=file,
	species="H"
)

#' Remove alternate/deprecated HPO term IDs and swap for new ones
#'
#' @template hpo.terms
#' @template terms
#' @param remove.dead Boolean to indicate whether to strip out terms which can't be found in the given hpo.terms database argument
#' @return A directed adjacency matrix of \code{terms} based on DAG structure of HPO, whereby each term is considered adjacent to it's MRCA in \code{terms}
#' @examples
#' data(hpo.terms)
#' swap.out.alt.ids(hpo.terms, c("HP:0001873"))
#' @export
#' @import magrittr
swap.out.alt.ids <- function(hpo.terms, terms, remove.dead=FALSE) {
	need.swap <- terms %in% names(hpo.terms$alt.id)
	terms[need.swap] <- hpo.terms$alt.id[terms[need.swap]]
	if (remove.dead) 
		terms <- terms[terms %in% hpo.terms$id]
	terms
}

#' Get list of character vector of HPO terms, given character vector of comma separated terms
#'
#' @param character.vector Character vector of comma separated terms
#' @return List of character vectors of HPO terms
#' @examples
#' term.set.list.from.character(c("HP:0001873", "HP:0001902"))
#' @export
#' @import magrittr
term.set.list.from.character <- function(character.vector) strsplit(
	gsub(" ", "", character.vector), 
	split=","
)

#' Get frequency of each term in a set of phenotypes
#'
#' @template hpo.terms
#' @template hpo.phenotypes
#' @param patch.missing Logical indicating whether to include all HPO even if they're not present in the \code{hpo.phenotypes} as if they had occurred once
#' @return Numeric vector of information contents, named by corresponding HPO terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
#' @seealso \code{\link{get.term.info.content}}
#' @examples
#' data(hpo.terms)
#' get.term.frequencies(hpo.terms, list("HP:0001873"))
#' @export
#' @import magrittr
get.term.frequencies <- function(
	hpo.terms, 
	hpo.phenotypes,
	patch.missing=FALSE
) {
	exp(-get.term.info.content(hpo.terms, hpo.phenotypes, patch.missing=FALSE))
}

#' Get information content of each term in a set of phenotypes
#'
#' @template hpo.terms
#' @template hpo.phenotypes
#' @param patch.missing Logical indicating whether to include all HPO even if they're not present in the \code{hpo.phenotypes} as if they had occurred once
#' @return Numeric vector of information contents, named by corresponding HPO terms. Takes into account ancestors, in the sense that all ancestor terms implied by the phenotypes are considered 'on'
#' @examples
#' data(hpo.terms)
#' get.term.info.content(hpo.terms, list("HP:0001873"))
#' @export
#' @import magrittr
get.term.info.content <- function(
	hpo.terms, 
	hpo.phenotypes,
	patch.missing=FALSE
) {
	terms.tab <- table(unlist(lapply(hpo.phenotypes, function(x) get.ancestors(hpo.terms, x))))
	total.patients <- length(hpo.phenotypes)
	terms.numeric <- as.numeric(terms.tab)
	names(terms.numeric) <- names(terms.tab)

	result <- log(total.patients) - ifelse(terms.numeric==0, log(total.patients), log(terms.numeric))
	
	if (patch.missing) {
		#include missing terms and give max information content...
		missing.terms <- setdiff(hpo.terms$id, names(result))
		missing.infos <- rep(max(result), length(missing.terms))
		names(missing.infos) <- missing.terms
		result <- c(result, missing.infos) 
	}

	result
}

#' @name hpo.terms
#' @title HPO Terms object (based on version 887 of the HPO)
#' @docType data
#' @format List of indices containing metadata and structure of HPO 
NULL

#' @name mpo.terms
#' @title MPO Terms object
#' @docType data
#' @format List of indices containing metadata and structure of MPO 
NULL

#' @name mpo.to.hpo
#' @title Object containing data for mapping between MPO and HPO
#' @docType data
#' @format List of HPO terms per MPO term
NULL

