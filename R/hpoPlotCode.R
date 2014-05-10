library(Hmisc)
library(functional)
library(plotrix)
library(Rgraphviz)

setDimNames <- function(array.object, list.of.dimension.names) {
	dimnames(array.object) <- list.of.dimension.names
	array.object
}

#function to remove implied ancestors from a set of hpo terms
clean.terms <- function(hpo.terms, terms) {
	redundant <- unlist(
		lapply(
			terms,
			function(x) setdiff(hpo.terms$ancestors[[x]], x)
		)
	)
	setdiff(terms,redundant)
}

#get ancestors for a set of terms...
get.ancestors <- function(hpo.terms, target.terms) {
	unique(
		unlist(
			lapply(target.terms, function(x) hpo.terms$ancestors[[x]])
		)
	)
}

get.descendants <- function(hpo.terms, ancestor) unique(
	c(
		ancestor,
		do.call(
			c,
			lapply(hpo.terms$children[[ancestor]], function(child.term) get.descendants(hpo.terms, child.term))
		)
	)
)

#get a matrix with columns of hpo terms and rows of patients,
#entry for a patient/hpo term = 1 if the patient has the term and 0 otherwise. 
get.term.patient.matrix <- function(patient.hpo.terms.with.ancs) {
	all.terms <- unique(unlist(patient.hpo.terms.with.ancs))
	result <- t(sapply(patient.hpo.terms.with.ancs, function(x) all.terms %in% x))
	colnames(result) <- all.terms
	return(result)
}


get.hpo.terms <- function(
	file="http://compbio.charite.de/hudson/job/hpo/lastStableBuild/artifact/ontology/release/hp.obo"
) {
	cat("Creating hpo.terms database...\n")
	hpo.obo.lines <- readLines(file)

	hpo.term.id.pattern <- "^id: (HP:\\d+)"
	hpo.term.name.pattern <- "^name: (\\.*)"
	hpo.term.parent.pattern <- "^is_a: (HP:\\d+)"
	hpo.term.pattern <- "HP:\\d+"
	hpo.term.alt.id.pattern <- "^alt_id: (HP:\\d+)"

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

swap.out.alt.ids <- function(hpo.terms, hpo.term.set, remove.dead=FALSE) {
	need.swap <- hpo.term.set %in% names(hpo.terms$alt.id)
	hpo.term.set[need.swap] <- hpo.terms$alt.id[hpo.term.set[need.swap]]
	if (remove.dead) 
		hpo.term.set <- hpo.term.set[hpo.term.set %in% hpo.terms$id]
	return(hpo.term.set)
}

term.set.list.from.character <- function(
	hpo.terms, 
	character.vector,
	remove.empty=FALSE
) {
	processed <- lapply(
		strsplit(
			gsub(" ", "", character.vector), 
			split=","
		),
		Compose(
			function(x) swap.out.alt.ids(hpo.terms, x, remove.dead=FALSE),
			function(x) clean.terms(hpo.terms, x)
		)
	)

	if (remove.empty)
		return(Filter(x=processed, f=function(x) length(x) > 0))
	else
		return(processed)
}

get.cohort.information <- function(
	hpo.terms, 
	patient.hpo.terms,
	patch.missing=FALSE
) {
	terms.tab <- table(unlist(lapply(patient.hpo.terms, function(x) get.ancestors(hpo.terms, x))))
	total.patients <- length(patient.hpo.terms)
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

	return(result)
}

remove.duplicate.term.sets <- function(hpo.terms, hpo.terms.set) {
	term.patient.matrix <- get.term.patient.matrix(
		lapply(hpo.terms.set, function(x) get.ancestors(hpo.terms, x))
	)
	return(
		apply(
			unique(data.frame(term.patient.matrix)),
			1,
			function(x) clean.terms(
				hpo.terms,
				colnames(term.patient.matrix)[x]
			)
		)
	)
}

get.term.pseudo.adjacency.matrix <- function(hpo.terms, terms) setDimNames(
	t(
		sapply(
			setNames(terms, terms),
			function(term) "%in%"(
				setNames(terms, terms),
				clean.terms(	
					hpo.terms, 
					setdiff(
						intersect(terms, hpo.terms$ancestors[[term]]),
						term
					)
				)
			)
		)
	),
	rep(
		list(terms),
		2
	)
)

get.term.adjacency.matrix <- function(hpo.terms, terms) {
	names(terms) <- terms
	adj.mat <- sapply(
		terms,
		function(term) terms %in% hpo.terms$parents[[term]]
	)
	rownames(adj.mat) <- colnames(adj.mat)
	adj.mat <- t(adj.mat)
	return(adj.mat)
}

suggest.aspect.ratio.width.to.height <- function(hpo.terms, list.of.patients.or.vector.of.terms) {
	leaves <- clean.terms(hpo.terms, unlist(list.of.patients.or.vector.of.terms))
	(120 + 7 * length(leaves))/120
}

remove.uninformative.terms <- function(hpo.terms, patient.hpo.terms) {
	#loop through terms, removing parents where their children are possessed by the same people - then clean up and link to implied nodes...
	with.ancs <- lapply(patient.hpo.terms, function(x) get.ancestors(hpo.terms, x))
	all.terms <- unique(unlist(with.ancs))
	terms <- Filter(
		x=unique(unlist(with.ancs)),
		f=function(term) {
			if (
				"=="(
					length(
						intersect(
							hpo.terms$children[[term]],
							all.terms
						)
					),
					0
				)
			) return(TRUE)

			patients.of.each.child <- lapply(
				intersect(hpo.terms$children[[term]], all.terms),
				function(child) sapply(
					with.ancs, 
					function(patient.terms) child %in% patient.terms
				)
			)
			patients.of.term <- sapply(
				with.ancs,
				function(patient.terms) term %in% patient.terms
			)
			!all(
				sapply(
					patients.of.each.child,
					function(x) identical(
						x,
						patients.of.term
					)
				)
			)
		}
	)
	
	terms
}

	
get.node.friendly.long.names <- function(hpo.terms, terms) {
	#get rid of the "Abnormality of"s because they take up too much room...
	simpleCap <- function(x) {
		s <- strsplit(x, " ")[[1]]
		paste(
			toupper(substring(s, 1,1)), 
			substring(s, 2),
			sep="", collapse=" "
		)
	}

	reorglabs <- sapply(
		gsub(
			"(Abnormality of (the )?)|(Abnormal)", 
			"", 
			hpo.terms$name[terms]
		),
		simpleCap
	)

	reorglabs <- sapply(
		reorglabs, 
		function(x) {
			words <- strsplit(x, split=" |-")[[1]]
			if (length(words) == 1)
				return(words)
			
			lines <- list(words[1])
			for (word.no in 2:length(words))
				if (nchar(paste(c(words[word.no], lines[[length(lines)]]), collapse=" ")) > 17)
					lines <- c(lines, words[word.no])
				else
					lines[[length(lines)]] <- c(lines[[length(lines)]], words[word.no])

			desc.lines <- paste(
				lapply(lines, function(line) paste(line, collapse=" ")),
				collapse="\\\n"
			)

			desc.lines
		}
	)

	reorglabs
}

colouring.functions <- list(
	get.patient.based.colours=function(hpo.terms, plotting.context) {
		patient.list <- plotting.context$patient.list
		ancestors.by.patient <- lapply(patient.list, function(x) get.ancestors(hpo.terms, x))
		terms <- plotting.context$terms
		
		term.pat.mat <- data.frame(t(get.term.patient.matrix(ancestors.by.patient)))

		patient.combos <- unique(term.pat.mat)

		colours <- rainbow(nrow(patient.combos), alpha=0.5)
		
		node.colors <- setNames(
			colours[
				apply(
					term.pat.mat[terms,], 
					1, 
					function(term.patients) which(
						apply(
							patient.combos, 
							1, 
							function(combo) identical(
								as.logical(term.patients), 
								as.logical(combo)
							)
						)
					)
				)
			],
			terms
		)

		list(
		colours=node.colors,
			legend.callback=Curry(
				legend,
				legend=apply(patient.combos, 1, function(x) paste(
					colnames(patient.combos)[x],
					collapse=" & "
				)),
				col=colours,
				x="bottomright",
				title=NULL,
				lwd=2
			)
		)
	}, 
	get.frequency.based.colours=function(
		hpo.terms, 
		plotting.context, 
		colour.func=NULL
	) {
		if (is.null(colour.func))
			colour.func <- colorRampPalette(c("#0099FF", "Green", "Yellow"))

		ancestors.by.patient <- lapply(plotting.context$patient.list, function(x) get.ancestors(hpo.terms, x))
		
		patients.with.term.count <- sapply(
			plotting.context$terms,
			function(term) sum(
				sapply(ancestors.by.patient, function(ancs) term %in% ancs)
			)
		)

		node.colours <- colour.func(1+diff(range(patients.with.term.count)))[
			patients.with.term.count-min(patients.with.term.count)+1
		]
		names(node.colours) <- plotting.context$terms
		
		#strip out the terms not relevant to this plot...
		node.colours <- node.colours[which(names(node.colours) %in% plotting.context$terms)]

		n <- length(plotting.context$patient.list)	
		legend.colours <- colour.func(1+diff(range(patients.with.term.count)))

		freq.range <- seq(
			min(patients.with.term.count), 
			max(patients.with.term.count)
		)	
		names(legend.colours) <- paste(
			freq.range,	
			ifelse(freq.range==1,"patient","patients")
		)

		list(
			colours=node.colours,
			legend.callback=NULL
		)
	}, 
	get.information.based.colours=function(
		hpo.terms, 
		plotting.context, 
		colourPalette=colorRampPalette(c("Yellow", "Green", "#0099FF"))(10)
	) {
		stopifnot(!is.null(plotting.context$information))

		terms.freq <- exp(-plotting.context$information[plotting.context$terms])

		freq.groups <- cut(terms.freq, length(colourPalette))

		if (diff(range(terms.freq)) == 0)
			node.colors <- rep(colourPalette[1], length(terms.freq))
		else
			node.colors <- colourPalette[as.integer(freq.groups)]

		names(node.colors) <- plotting.context$terms
		list(
			colours=node.colors,
			legend.callback=Curry(
				function(xl,yb,xr,yt,rect.col,legend,gradient) {
					color.legend(
						xl=xl,
						yb=yb,
						xr=xr,
						yt=yt,
						rect.col=rect.col,
						legend=legend,
						gradient=gradient
					)
					text(
						x=as.integer(xr/2), 
						y=65, 
						labels="Proportion Of Cohort"
					)
				},
				xl=0, 
				yb=0, 
				xr=as.integer(
					dev.size(units="px")[1]/3
				), 
				yt=30, 
				rect.col=colourPalette,
				legend=paste(
					round(
						((0:5) * (max(terms.freq)-min(terms.freq))/5 + min(terms.freq)) * 100
					),
					"%",
					sep=""
				),
				gradient="x"
			)
		)
	},
	get.white.nodes=function(hpo.terms, plotting.context) list(
		colours=setNames(rep("White", length(plotting.context$terms)), plotting.context$terms),
		legend.colours=NULL
	)
)

size.functions <- list(
	get.significance.based.sizes=function(hpo.terms, plotting.context) {
		stopifnot(!is.null(plotting.context$information))
		num.in.group <- sapply(
			plotting.context$terms,
			function(term) sum(
				sapply(
					lapply(
						plotting.context$patient.list, 
						function(x) get.ancestors(hpo.terms, x)
					), 
					function(ancs) term %in% ancs
				)
			)
		)

		cohort.freq <- exp(-plotting.context$information[plotting.context$terms])
		
		setNames(
			(
				function(x, high, low) "+"(
					low,
					"*"(
						"/"(
							x-min(x),
							diff(range(x))
						),
						high-low
					)
				)		
			)(
				-pbinom(
					q=num.in.group-1,
					size=length(plotting.context$patient.list),
					prob=cohort.freq,
					lower.tail=FALSE,
					log.p=TRUE
				),
				3,
				1
			),
			plotting.context$terms
		)
	}, 
	get.frequency.based.sizes=function(hpo.terms, plotting.context) {
		group.freq <- sapply(
			plotting.context$terms,
			function(term) sum(
				sapply(
					lapply(
						plotting.context$patient.list, 
						function(x) get.ancestors(hpo.terms, x)
					), 
					function(ancs) term %in% ancs
				)
			)
		)

		setNames(
			(
				function(x, high, low) "+"(
					low,
					"*"(
						"/"(
							x-min(x),
							diff(range(x))
						),
						high-low
					)
				)		
			)(group.freq, 3, 1),
			plotting.context$terms
		)
	},
	get.standard.sizes=function(hpo.terms, plotting.context) setNames(rep(0.75, length(plotting.context$terms)), plotting.context$terms)
)

labelling.functions <- list(
	get.frequency.based.labels=function(hpo.terms, plotting.context) {
		node.friendly <- get.node.friendly.long.names(hpo.terms, plotting.context$terms)
		freqs <- sapply(
			plotting.context$terms,
			function(term) sum(
				sapply(lapply(plotting.context$patient.list, function(x) get.ancestors(hpo.terms, x)), function(ancs) term %in% ancs)
			)
		)

		result <- "if"(
			!is.null(plotting.context$information),
			paste(
				plotting.context$terms,
				node.friendly,
				paste(
					round(100 * exp(-plotting.context$information[plotting.context$terms])),	
					"% Of Cohort", 
					sep=""
				),
				paste(
					freqs,
					" / ",
					length(plotting.context$patient.list),
					sep=""
				),
				sep="\\\n"
			),
			paste(
				plotting.context$terms,
				node.friendly,
				paste(
					freqs,
					" / ",
					length(plotting.context$patient.list),
					sep=""
				),
				sep="\\\n"
			)
		)
		
		setNames(
			result,
			plotting.context$terms
		)
	}, 
	get.simple.node.labels=function(hpo.terms, plotting.context) setNames(
		get.node.friendly.long.names(hpo.terms, plotting.context$terms),
		plotting.context$terms
	), 
	get.informative.node.labels=function(hpo.terms, plotting.context) {
		stopifnot(!is.null(plotting.context$information))
		setNames(
			paste(
				plotting.context$terms,
				get.node.friendly.long.names(hpo.terms, plotting.context$terms),
				paste(
					round(100 * exp(-plotting.context$information[plotting.context$terms])),	
					"%", 
					sep=""
				),
				sep="\\\n"
			), 
			plotting.context$terms
		)
	},
	get.patient.based.labels=function(hpo.terms, plotting.context) setNames(
		paste(
			plotting.context$terms,
			get.node.friendly.long.names(hpo.terms, plotting.context$terms),
			sapply(
				plotting.context$terms,
				function(term) paste(
					names(
						Filter(
							x=plotting.context$patient.list,
							f=function(patient.terms) term %in% get.ancestors(
								hpo.terms,
								patient.terms
							)
						)
					),
					collapse=","
				)
			),
			sep="\\\n"
		),
		plotting.context$terms
	),
	get.code.node.labels=function(hpo.terms, plotting.context) setNames(plotting.context$terms, plotting.context$terms)
)

border.functions <- list(
	get.grey.borders=function(hpo.terms, plotting.context) setNames(rep("#00000080", length(plotting.context$terms)), plotting.context$terms),
	get.no.borders=function(hpo.terms, plotting.context) setNames(rep("#00000000", length(plotting.context$terms)), plotting.context$terms)
)

hpo.graph <- function(
	hpo.terms,
	patients,
	term.information.content=NULL,
	term.population.frequency=NULL,
	terms=NULL,
	colouring.function=colouring.functions$get.frequency.based.colours,
	labelling.function=labelling.functions$get.frequency.based.labels,
	border.function=border.functions$get.no.borders,
	size.function=size.functions$get.standard.sizes,
	main.title=NULL,
	draw.legend=FALSE,
	filter.out.uninformative=FALSE,
	filter.out.non.pa.terms=FALSE,
	min.occs=1,
	nodeAttrs=NULL,
	pdf.file.name=NULL,
	pdf.height=10,
	pdf.width=NULL
) {
	select.n.most.frequent.terms <- function(terms, hpo.terms, plotting.context, n) intersect(
		terms,
		Compose(	
			unlist,
			table, 
			function(x) sort(x, decreasing=TRUE), 
			names
		)(
			lapply(
				plotting.context$patient.list,
				function(x) get.ancestors(hpo.terms, x)
			)
		)[1:min(n, length(terms))]
	)

	remove.terms.with.less.than.n.occurrences <- function(terms, hpo.terms, plotting.context, n) intersect(
		terms,
		Compose(unlist, table, function(x) which(x >= n), names)(
			lapply(
				plotting.context$patient.list,
				function(x) get.ancestors(hpo.terms, x)
			)
		)
	)

	remove.uninformative.for.plot <- function(terms, hpo.terms, plotting.context) {
		intersect(
			remove.uninformative.terms(
				hpo.terms,
				lapply(
					plotting.context$patient.list,
					function(x) clean.terms(
						hpo.terms, 
						intersect(
							get.ancestors(hpo.terms, x), 
							terms
						)
					)
				)
			),
			terms
		)
	}

	remove.non.pa.terms <- function(terms, hpo.terms, plotting.context) Filter(
		x=terms,
		f=Curry(
			function(is.this, a.descendant.of.this) "%in%"(
				a.descendant.of.this,
				get.ancestors(
					hpo.terms,
					is.this	
				)
			),
			a.descendant.of.this=hpo.terms$id[hpo.terms$name == "Phenotypic abnormality"]
		)
	)

	get.hpo.graph <- function(
		hpo.terms,
		plotting.context,
		colouring.function,
		labelling.function,
		border.function,
		size.function,
		nodeAttrs=NULL
	) {
		adj.mat <- get.term.pseudo.adjacency.matrix(
			hpo.terms,
			plotting.context$terms
		)

		get.hpo.graph.rgraphviz <- function(adj.mat) new(
			"graphAM", 
			adjMat=adj.mat, 
			edgemode="directed"
		)

		hpo.graph <- get.hpo.graph.rgraphviz(adj.mat)

		node.colouring <- colouring.function(hpo.terms, plotting.context)
		
		if (is.null(nodeAttrs)) {
			attrs <- list(
				fontsize=setNames(
					rep(
						"if"(
							is.null(plotting.context$font.size),
							30,
							plotting.context$font.size
						),
						length(plotting.context$terms)
					),
					plotting.context$terms
				),
				shape=setNames(
					rep(
						"if"(
							is.null(plotting.context$shape),
							"circle",
							plotting.context$shape
						),
						length(plotting.context$terms)
					),
					plotting.context$terms
				),
				width=size.function(hpo.terms, plotting.context),
				color=border.function(hpo.terms, plotting.context),
				fillcolor=node.colouring$colours,
				label=labelling.function(hpo.terms, plotting.context)
			)
		} else {
			attrs <- nodeAttrs
		}

		result <- agopen(graph=hpo.graph, nodeAttrs=attrs, name="HPO.term.plot") 
		if (length(result@AgEdge) > 0)
			for (i in 1:length(result@AgEdge)) {
				result@AgEdge[[i]]@lty <- ifelse(
					result@AgEdge[[i]]@head %in% hpo.terms$parents[[result@AgEdge[[i]]@tail]],
					"solid",
					"dashed"
				)
			}

		result
	}

	if (class(colouring.function) != "function") {
		stopifnot(!is.null(names(colouring.function)))
		colouring.function <- Curry(
			function(hpo.terms, plotting.context, result) list(colours=result),
			result=colouring.function
		)
	}
	if (class(labelling.function) != "function") {
		stopifnot(!is.null(names(labelling.function)))
		labelling.function <- Curry(
			function(hpo.terms, plotting.context, result) result,
			result=labelling.function
		)
	}
	if (class(size.function) != "function") {
		stopifnot(!is.null(names(size.function)))
		size.function <- Curry(
			function(hpo.terms, plotting.context, result) result,
			result=size.function
		)
	}
	if (class(border.function) != "function") {
		stopifnot(!is.null(names(border.function)))
		border.function <- Curry(
			function(hpo.terms, plotting.context, result) result,
			result=border.function
		)
	}
	
	plotting.context <- list(
		patient.list=patients,
		information="if"(
			!is.null(term.information.content),
			term.information.content,
			"if"(
				!is.null(term.population.frequency),
				-log(term.information.content),
				NULL
			)
		)
	)

	if (!is.null(terms))
		plotting.context$terms <- terms
	else {
		term.filters <- list()

		if (filter.out.non.pa.terms)
			term.filters <- c(term.filters, remove.non.pa.terms)

		term.filters <- c(
			term.filters, 
			Curry(remove.terms.with.less.than.n.occurrences, n=min.occs)
		)

		if (filter.out.uninformative) 
			term.filters <- c(term.filters, remove.uninformative.for.plot)

		plotting.context$terms <- Reduce(
			f=function(terms, filter.func) filter.func(terms=terms, hpo.terms=hpo.terms, plotting.context=plotting.context),
			init="if"(
				is.null(plotting.context$terms),
				Compose(unlist, unique)(
					lapply(
						plotting.context$patient.list,
						function(x) get.ancestors(hpo.terms, x)
					)
				),
				plotting.context$terms
			),
			x=term.filters
		)

	}

	if (!is.null(pdf.file.name)) {
		pdf.width <- "if"(
			is.null(pdf.width),
			pdf.height * suggest.aspect.ratio.width.to.height(hpo.terms, plotting.context$terms),
			pdf.width
		)
		pdf(file=pdf.file.name, height=pdf.height, width=pdf.width)
	}

	graph <- get.hpo.graph(
		hpo.terms=hpo.terms,
		plotting.context=plotting.context,
		colouring.function=colouring.function,
		labelling.function=labelling.function,
		border.function=border.function,
		size.function=size.function,
		nodeAttrs=nodeAttrs
	)
	
	plot(
		graph,
		main=main.title
	)

	node.colouring <- colouring.function(hpo.terms, plotting.context)
	if ((!is.null(node.colouring$legend.callback)) & (draw.legend==TRUE))
		node.colouring$legend.callback()

	if (!is.null(pdf.file.name)) 
		dev.off()
}
