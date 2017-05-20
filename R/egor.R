#  File R/egodata.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################


#' Convert \code{\link{egodata}} Objects to \code{\link{egor}} Objects
#' 
#' @param object a \code{\link{egodata}} object
#'
#' @return An \code{\link{egor}} object.
#' 
#' @section Miscellaneous Methods: The following \dQuote{standard} methods have
#' @author Pavel N. Krivitsky
#' @keywords manip methods
#' @export
as.egor.egodata <- function(object, ...){
  ego.design <- list(~1, weights = rep(object$egoWt, length.out=nrow(object$egos)))
  egor(egos.df=object$egos, alters.df=object$alters, egoID = object$egoIDcol, ego.design=ego.design)
}

#' Construct an Egocentric View of a \code{\link{network}} Object
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
#' @return An \code{\link{egor}} object.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link{as.network.egor}}, which performs the inverse
#' operation (though drops the ties).
#' @keywords datagen manip
#' @examples
#' 
#' # See example(ergm.ego) and example(as.network.egor).
#' @import tibble
#' @export
as.egor.network<-function(object,special.cols=c("na")){
  N<-network.size(object)

  egos<-list()
  
  for(a in list.vertex.attributes(object))
    if(!(a %in% special.cols)) egos[[a]]<-get.vertex.attribute(object,attrname=a)

  egos <- tibble::as_tibble(egos)

  el<-as.edgelist(object)
  el<-rbind(el,el[,2:1])
  alterS<-tapply(el[,2],INDEX=el[,1],FUN=c,simplify=FALSE)

  alters <- lapply(seq_len(N), get.neighborhood, x=object) # so v gets the index variable

  alters <- lapply(alters, function(js) egos[js,,drop=FALSE])

  egor(egos.df=egos,alters.df=alters)
}


#' Construct an Empty ``Template'' Network Consistent with an Egocentric Sample
#' 
#' Taking a \code{\link{egor}} object, constructs a
#' \code{\link[network]{network}} object with no edges whose vertices have the
#' attributes of the egos in the dataset, replicating the egos as needed, and
#' taking into accounts their sampling weights.
#' 
#' 
#' @param x A \code{\link{egor}} object.
#' @param N The target number of vertices the output network should have.
#' @param scaling If \code{\link{egor}} contains weights or \code{N} is not
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
#' @seealso \code{\link{as.egor.network}}, which performs the inverse
#' operation.
#' @keywords manip
#' @examples
#' 
#' 
#' data(faux.mesa.high)
#' summary(faux.mesa.high, print.adj = FALSE)
#' 
#' fmh.ego <- as.egor(faux.mesa.high)
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
as.network.egor<-function(x, N, scaling=c("round","sample"), ...){
  scaling <- match.arg(scaling)
  w <- weights(x)
  egoinds <- switch(scaling,
                    greedy={
                      .greedy.scaling(N,w)
                    },
                    round={
                      .round.scaling(N,w)
                    },
                    sample={
                      sample.int(length(w),N,TRUE,w)
                    })

  N <- length(egoinds) # round scaling may modify N.
  y0<-network.initialize(N,directed=FALSE)

  x <- x[egoinds,]
  
  for(ego.col in setdiff(names(x),c(".alters",".alter_ties")))
    if(is.factor(x[[ego.col]]))
      y0 <- set.vertex.attribute(y0,ego.col,as.character(x[[ego.col]]))
    else
      y0 <- set.vertex.attribute(y0,ego.col,x[[ego.col]])
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

# Not really a generic function, but perhaps should be.

#' @export
sample <- function(x, size, replace=FALSE, prob=NULL, ...) UseMethod("sample")
#' @export
sample.default <- function(x, ...) base::sample(x, ...)

#' @export
sample.egor <- function(x, size, replace=FALSE, prob=NULL, ...){
  if(missing(size)) size <- nrow(x)
  
  is <- sample.int(nrow(x), size, replace, prob)

  x[is, ,aspect="egos"]
}

