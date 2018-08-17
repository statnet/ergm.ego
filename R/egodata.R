#  File R/egodata.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################


#' Convert to or Construct \code{\link{egodata}} Objects
#' 
#' \code{\link{as.egodata}} is a generic function to construct
#' \code{\link{egodata}} objects from a variety of sources.
#' \code{\link{egodata}} function is the standard constructor, taking two data
#' frames. For other methods for this class, see the Miscellaneous Methods
#' section.
#' 
#' 
#' @aliases as.egodata as.egodata.data.frame egodata.object egodata
#' na.omit.egodata dim.egodata dimnames.egodata sample.egodata head.egodata
#' sample sample.default
#' @param object,egos The object from which the egocentric data should be
#' constructed. For the \code{\link{data.frame}} methods and
#' \code{\link{egodata}} itself, a data frame containing at least the column
#' named in \code{egoIDcol}, whose values must all be unique. Other columns
#' contain information about the egos.
#' @param alters A data frame containing at least the column named in
#' \code{egoIDcol}, whose values do not have to be unique, and not every ego
#' must be represented. Other columns contain information about the alters.
#' @param egoWt A vector of the same length as number of rows in \code{egos} or
#' \code{object}, containing the relative sampling weight of each ego.
#' @param \dots Additional arguments; currently unused.
#' @param egoIDcol Name of the column in the ego table containing the unique
#' ego identifier.
#' @return An \code{\link{egodata}} object. The object is a list containing the
#' following elements:
#' 
#' \item{egos}{ A data frame with one row for each ego, containing at least the
#' column named in \code{egoIDcol}, and other columns containing attributes of
#' the egos.  }
#' 
#' \item{alters}{ A data frame containing at least the column named in
#' \code{egoIDcol}, and other columns containing attributes of the alters.  }
#' 
#' \item{egoWt}{A vector of the same length as the number of egos, containing
#' the relative sampling weight of each ego.}
#' 
#' \item{egoIDcol}{ Name of the column in the ego table containing the unique
#' ego identifier.  }
#' @section Miscellaneous Methods: The following \dQuote{standard} methods have
#' also been implemented for \code{\link{egodata}}: \describe{
#' \item{"dim.egodata"}{A vector with three elements containing the
#' \dQuote{dimensions} of the \code{\link{egodata}} object: number of egos,
#' number of columns in the \code{egos} table, and number of columns in the
#' \code{alters} table, inclsive of the ego identifier column. As a corollary,
#' \code{\link{nrow}} returns the number of egos in the dataset.  }
#' 
#' \item{"dimnames.egodata"}{A list with three elements containing
#' the \dQuote{dimension names} of the \code{\link{egodata}} object: ego IDs,
#' column names of the \code{egos} table, and column names of the \code{alters}
#' table, inclsive of the ego identifier column.  }
#' 
#' \item{"sample.egodata"}{ As \code{\link{sample}}, but takes and
#' returns a simulated \code{\link{egodata}} dataset by resampling egos,
#' adjusting ego weights as necessary, if weighted sampling was used.  }
#' 
#' \item{"head.egodata"}{ As \code{\link{head}}, but returns the
#' first \code{n} rows of egos, alters, and weights.  }
#' 
#' \item{"na.omit.egodata"}{ As \code{\link{na.omit.data.frame}},
#' but takes and returns an \code{\link{egodata}} dataset, with egos with
#' \code{NA} in their rows or in their alters' rows. An optional argument
#' \code{relevant}, defaulting to all columns, can be used to select (by index
#' or name) based on which columns an ego may be dropped. (I.e., \code{NA}s in
#' those not \dQuote{relevant} are ignored.)  } }
#' @author Pavel N. Krivitsky
#' @seealso \code{\link{ergm.ego}} for examples,
#' \code{\link{as.network.egodata}}, \code{\link{as.egodata.network}},
#' \code{\link{subset.egodata}}, \code{\link{[.egodata}}
#' @keywords manip methods
#' @export
egodata <- function(egos, alters, egoWt=1, ..., egoIDcol="egoID"){
  egoWt <- rep(egoWt, length.out=nrow(egos))
  out <- list(egos=egos, alters=.prune.alters(egos, alters, egoIDcol), egoWt = egoWt, egoIDcol=egoIDcol)  
  class(out) <- "egodata"
  out
}

#' @rdname egodata
#' @export
as.egodata <- function(object, ..., egoIDcol="egoID"){
  UseMethod("as.egodata")
}

#' @rdname egodata
#' @export
as.egodata.data.frame <- function(object, alters, egoWt = 1, ..., egoIDcol="egoID"){
  egodata(egos=object, alters=alters, egoWt=egoWt, ..., egoIDcol=egoIDcol)
}

# Conduct an egocentric census from the undirected network y=,
# returning an egodata object. The corresponding vertex attributes of
# y= are copied into columns in these data frames, excluding
# attributes listed in special.cols=.


#' Construct an Egocentric View of a Network
#' 
#' Given a \code{\link[network]{network}} object, construct an
#' \code{\link{egodata}} object representing a census of all the actors in the
#' network. Used mainly for testing.
#' 
#' 
#' @param object A \code{\link[network]{network}} object.
#' @param special.cols Vertex attributes that should not be copied to the
#' \code{egos} and \code{alters} tables. Defaults to attributes special to the
#' \code{\link[network]{network}} objects.
#' @param \dots Additional arguments, currently unused.
#' @param egoIDcol The name of the vertex attribute containg unique ego IDs.
#' Defaults to \code{"vertex.names"}.
#' @return An \code{\link{egodata}} object.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link{as.network.egodata}}, which performs the inverse
#' operation (though drops the ties).
#' @keywords datagen manip
#' @examples
#' 
#' # See example(ergm.ego) and example(as.network.egodata).
#' 
#' @export
as.egodata.network<-function(object,special.cols=c("na","vertex.names"),...,egoIDcol="vertex.names"){
  N<-network.size(object)

  if(egoIDcol%in%list.vertex.attributes(object)){
    egoIDs<-object%v%egoIDcol
    if(any(duplicated(egoIDs))){
      warning("Non-unique ego IDs; using 1..N.")
      egoIDs <- seq_along(egoIDs)
    }
  }else{
    message("Network does not have vertex attribute ",sQuote(egoIDcol)," to use as ego ID; using 1..N.")
    egoIDs <- seq_len(N)
  }
  
  egos<-list()
  egos[[egoIDcol]]<-egoIDs
  
  for(a in list.vertex.attributes(object))
    if(!(a %in% special.cols)) egos[[a]]<-get.vertex.attribute(object,attrname=a)

  el<-as.edgelist(object)
  el<-rbind(el,el[,2:1])
  alterS<-unlist(tapply(el[,2],INDEX=el[,1],FUN=c,simplify=FALSE))
  alter.eID<-egoIDs[unlist(tapply(el[,1],INDEX=el[,1],FUN=c,simplify=FALSE))]
  
  alters<-list()

  alters[[egoIDcol]]<-alter.eID
    
  for(a in list.vertex.attributes(object))
    if(!(a %in% special.cols)) alters[[a]]<-get.vertex.attribute(object,attrname=a)[alterS]

  egodata(egos=as.data.frame(egos,stringsAsFactors=FALSE),alters=as.data.frame(alters,stringsAsFactors=FALSE), egoIDcol=egoIDcol)
}

.prune.alters <- function(egos, alters, egoIDcol){
  eis <- egos[[egoIDcol]]
  aeis <- alters[[egoIDcol]]

  todel <- !(aeis %in% eis)

  if(any(todel)) alters[!todel,,drop=FALSE]
  else alters
}



#' Construct an Empty ``Template'' Network Consistent with an Egocentric Sample
#' 
#' Taking a \code{\link{egodata}} object, constructs a
#' \code{\link[network]{network}} object with no edges whose vertices have the
#' attributes of the egos in the dataset, replicating the egos as needed, and
#' taking into accounts their sampling weights.
#' 
#' 
#' @param x A \code{\link{egodata}} object.
#' @param N The target number of vertices the output network should have.
#' @param scaling If \code{\link{egodata}} contains weights or \code{N} is not
#' a multiple of number of egos in the sample, it may not be possible, for a
#' finite \code{N} to represent each ego exactly according to its relative
#' weight, and \code{scaling} controls how the fractional egos are allocated:
#' \describe{ \item{"round"}{(the default) Rather than treating \code{N} as a hard
#' setting, calculate \eqn{N w_i / w_\cdot} for each ego \eqn{i} and round it
#' to the nearest integer. Then, the \code{N} actually used will be the sum of
#' these rounded freqencies.} \item{"sample"}{Resample in
#' proportion to \eqn{w_i}.} }
#' @param \dots Additional arguments, currently unused.
#' @return A \code{\link[network]{network}} object.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link{as.egodata.network}}, which performs the inverse
#' operation.
#' @keywords manip
#' @examples
#' 
#' 
#' data(faux.mesa.high)
#' summary(faux.mesa.high, print.adj = FALSE)
#' 
#' fmh.ego <- as.egodata(faux.mesa.high)
#' 
#' # Same actor attributes
#' fmh.template <- as.network(fmh.ego, N=network.size(faux.mesa.high))
#' summary(fmh.template, print.adj = FALSE)
#' 
#' # Twice the actors, same distribution
#' fmh2.template <- as.network(fmh.ego, N=2*network.size(faux.mesa.high))
#' summary(fmh2.template, print.adj = FALSE)
#'
#' @import network
#' @export
as.network.egodata<-function(x, N, scaling=c("round","sample"), ...){
  scaling <- match.arg(scaling)
  egoinds <- switch(scaling,
                    round={
                      .round.scaling(N,x$egoWt)
                    },
                    sample={
                      sample(length(x$egoWt),N,TRUE,x$egoWt)
                    })

  N <- length(egoinds) # round scaling may modify N.
  y0<-network.initialize(N,directed=FALSE)
  egos <- x$egos
  
  egos <- egos[egoinds,]
  
  for(ego.col in names(egos))
    if(is.factor(egos[[ego.col]]))
      y0 <- set.vertex.attribute(y0,ego.col,as.character(egos[[ego.col]]))
    else
      y0 <- set.vertex.attribute(y0,ego.col,egos[[ego.col]])
  y0 %v% ".ego.ind" <- egoinds
  y0
}

.greedy.scaling <- function(N, w){
  ideal<-N*w/sum(w)
  n<-floor(ideal) # "Guaranteed" assignments.
  r<-ideal-n
  leftover<-sum(r)
  if(leftover){
    best<-order(rank(r*w,ties.method="random"),decreasing=TRUE)[1:leftover]
    n[best]<-n[best]+1
  }
  rep(seq_along(w),n)
}

.round.scaling <- function(N, w){
  ideal<-N*w/sum(w)
  n<-round(ideal)
  rep(seq_along(w),n)
}


# Note: The following functions use parts of na.omit.data.frame() and
# subset.data.frame() from the R's stats package under the terms of
# the GNU GPL v3.

#' @rdname subset.egodata
#' @export
`[.egodata` <- function(x, i, j, ..., dup.action=c("make.unique", "fail", "number")){
  subset(x, i, j, ..., dup.action=dup.action)
}



#' Subsetting \code{\link{egodata}} Objects
#' 
#' Returns subsets of \code{\link{egodata}} objects that meet conditions.
#' 
#' 
#' @aliases subset.egodata [.egodata
#' @param x An \code{\link{egodata}} object.
#' @param subset,i An expression (evaluated in the context of the \code{egos}
#' table of \code{x} producing a logical, integer, or character vector
#' indicating which egos to select (and, for the latter two, how many times).
#' @param select,j A numeric or character vector specifying the columns of
#' \code{egos} and \code{alters} to select.
#' @param \dots Additional arguments, currently unused.
#' @param dup.action What to do when an ego is referenced multiple times:
#' \describe{ \item{"make.unique"}{Construct new unique ego IDs using
#' the \code{\link[base]{make.unique}} function } \item{"fail"}{Exit
#' with an error.} \item{"number"}{Number the egos consecutively in
#' the order they were selected } }
#' @return An \code{\link{egodata}} object.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link{sample.egodata}}
#' @keywords manip
#' @export
subset.egodata <- function(x, subset, select, ..., dup.action=c("make.unique", "fail", "number")){
  if(missing(subset)) subset <- TRUE
  
  if (missing(select)) 
    egovars <- altervars <- TRUE
  else {
    if((is.numeric(select) || is.logical(select) || is.integer(select))
       && !isTRUE(all.equal(names(x$egos),names(x$alters)))){
      stop("Logical or numeric column subsets cannot be used if column orderings for egos and alters are not identical. Use column names instead.")
    }
    
    nl <- as.list(seq_along(x$egos))
    names(nl) <- names(x$egos)
    vars <- eval(substitute(select), nl, parent.frame())
    eIDi <- switch(mode(vars),
                   numeric = which(names(x$egos)==x$egoIDcol),
                   character = x$egoIDcol)                   
    egovars <- union(eIDi,vars)

    nl <- as.list(seq_along(x$alters))
    names(nl) <- names(x$alters)
    vars <- eval(substitute(select), nl, parent.frame())
    eIDi <- switch(mode(vars),
                   numeric = which(names(x$alters)==x$egoIDcol),
                   character = x$egoIDcol)
    altervars <- union(eIDi,vars)
  }

  egoIDs <- switch(mode(subset),
                   numeric =,
                   logical =,
                   integer = x$egos[[x$egoIDcol]][subset],
                   character = x$egos[[x$egoIDcol]][match(subset, x$egos[[x$egoIDcol]])])

  dup.action <- match.arg(dup.action)
  unique.egoIDs <- switch(dup.action,
                          fail=stop("Selected subset calls for duplicating egos."),
                          numeric=seq_along(egoIDs),
                          make.unique=make.unique(as.character(egoIDs)))
  
  egos <- cbind(x$egos[match(egoIDs,x$egos[[x$egoIDcol]]),if(is.character(egovars)) intersect(egovars,names(x$egos)) else egovars,drop=FALSE], .unique.egoIDs = unique.egoIDs, stringsAsFactors=FALSE)
  alters <- merge(egos[c(x$egoIDcol,".unique.egoIDs")], x$alters[if(is.character(altervars)) intersect(altervars,names(x$alters)) else altervars], by=x$egoIDcol)
  egoWt <- x$egoWt[match(egoIDs,x$egos[[x$egoIDcol]])]

  # If we have duplicated egoIDs, we have to handle it per dup.action.
  
  if(any(duplicated(egoIDs))){
    egos[[x$egoIDcol]] <- egos[[".unique.egoIDs"]]
    alters[[x$egoIDcol]] <- alters[[".unique.egoIDs"]]
  }
  egos[[".unique.egoIDs"]] <- NULL
  alters[[".unique.egoIDs"]] <- NULL

  out <- x
  out$egos <- egos
  out$alters <- alters
  out$egoWt <- egoWt
  
  class(out) <- "egodata"
  out
}

# TODO: A more efficient implementation of this.
#' @importFrom stats na.omit
#' @export
na.omit.egodata <- function(object, relevant=TRUE, ...){
  # Create a subdataset containing only the relevant variables.
  obj <- subset(object,select=relevant)
  
  n <- length(obj$egos)
  omit <- FALSE
  vars <- seq_len(n)
  for (j in vars) {
    x <- obj$egos[[j]]
    if (!is.atomic(x)) 
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L) 
      omit <- omit | x
    else for (ii in 1L:d[2L]) omit <- omit | x[, ii]
  }
  ego.omit <- obj$egos[[obj$egoIDcol]][omit]

  n <- length(obj$alters)
  omit <- FALSE
  vars <- seq_len(n)
  for (j in vars) {
    x <- obj$alters[[j]]
    if (!is.atomic(x)) 
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L) 
      omit <- omit | x
    else for (ii in 1L:d[2L]) omit <- omit | x[, ii]
  }
  
  alter.omit <- obj$alters[[obj$egoIDcol]][omit]

  subset(object, -match(union(ego.omit,alter.omit),object$egos[[object$egoIDcol]]))  
}

#' @export
dimnames.egodata <- function(x){
  list(x$egos[[x$egoIDcol]], names(x$egos), names(x$alters))
}

#' @export
dim.egodata <- function(x){
  c(nrow(x$egos), ncol(x$egos), ncol(x$alters))
}

# Not really a generic function, but perhaps should be.

#' @export
sample <- function(x, size, replace=FALSE, prob=NULL, ...) UseMethod("sample")
#' @export
sample.default <- function(x, ...) base::sample(x, ...)

#' @export
sample.egodata <- function(x, size, replace=FALSE, prob=NULL, ...){
  if(missing(size)) size <- nrow(x)
  
  is <- sample.int(nrow(x), size, replace, prob)

  out <- subset(x, is)

  if(is.null(prob)) prob <- rep(1, nrow(x))
  
  out$egoWt <- x$egoWt[is]/prob[is]
  out
}

#' @importFrom utils head
#' @export
head.egodata <- function(x, n=6L, ...) lapply(x, head, n=n, ...)
