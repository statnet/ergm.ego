#  File R/egor.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################


#' Convert (deprecated) \code{\link{egodata}} Objects to
#' \code{\link{egor}} Objects
#'
#' @aliases egodata
#' 
#' @param x a \code{\link{egodata}} object
#' @param ... additional arguments, currently unused.
#'
#' @return An \code{\link{egor}} object.
#' 
#' @author Pavel N. Krivitsky
#' @keywords manip methods
#' @import egor
#' @export
as.egor.egodata <- function(x, ...){
  x$egos$.egoWt <- rep(x$egoWt, length.out=nrow(x$egos))
  ego_design <- list(weights = ".egoWt")
  x$alters$.alterIDdummy <- seq_len(nrow(x$alters)) # This is a workaround around egor() requiring alter ID even without alter-alter ties.
  egor(egos=x$egos, alters=x$alters, ID.vars = list(ego = x$egoIDcol, alter=".alterIDdummy"), ego_design=ego_design)
}

#' @rdname as.egor.egodata
#' @export
as_egor.egodata <- as.egor.egodata

#' Construct an Egocentric View of a \code{\link{network}} Object
#' 
#' Given a \code{\link[network]{network}} object, construct an
#' \code{\link{egor}} object representing a census of all the actors in the
#' network. Used mainly for testing.
#' 
#' 
#' @param x A \code{\link[network]{network}} object.
#' @param special.cols Vertex attributes that should not be copied to the
#' \code{egos} and \code{alters} tables. Defaults to attributes special to the
#' \code{\link[network]{network}} objects.
#' @param ... Additional arguments, currently unused.
#' @return An \code{\link{egor}} object.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link{template_network}}, which performs the inverse
#' operation (though drops the ties).
#' @keywords datagen manip
#' @examples
#' 
#' # See example(ergm.ego) and example(template_network).
#' @import tibble
#' @importFrom dplyr bind_rows bind_cols
#' @export
as.egor.network<-function(x,special.cols=c("na"),...){
  N<-network.size(x)

  egos<-list()
  
  for(a in list.vertex.attributes(x))
    if(!(a %in% special.cols)) egos[[a]]<-get.vertex.attribute(x,attrname=a)

  egos <- as_tibble(egos)

  # FIXME: Save edge attributes as well.
  alters <- lapply(seq_len(N), get.neighborhood, x=x) # so v gets the index variable

  aaties <- lapply(alters, lapply, get.neighborhood, x=x) # list of lists of alters' nominations

  # Note: Alter ID API is subject to change.
  aaties <- mapply(function(i, a, aa){
    # Only keep alters' ties that are with another alter of this ego.
    aa <- lapply(aa, function(ks) ks[ks %in% a])
    # FIXME: Save edge attributes as well.
    tibble(..Source=rep(a, sapply(aa, length)), ..Target=as.vector(unlist(aa), mode=storage.mode(a)), ..EgoID=rep.int(i, length(unlist(aa))))
  }, seq_along(aaties), alters, aaties, SIMPLIFY=FALSE) %>% bind_rows

  alters <- mapply(function(i, js) bind_cols(egos[js,,drop=FALSE], ..AlterID=js, ..EgoID=rep.int(i,length(js))), seq_along(alters), alters, SIMPLIFY=FALSE) %>% bind_rows

  egos <- bind_cols(egos, ..EgoID=seq_len(N))
  
  egor(egos=egos, alters=alters, aaties=aaties,
       ID.vars=list(ego="..EgoID", alter="..AlterID", source="..Source", target="..Target"))
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
#' @param ... Additional arguments, currently unused.
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
#' fmh.template <- template_network(fmh.ego, N=network.size(faux.mesa.high))
#' summary(fmh.template, print.adj = FALSE)
#' 
#' # Twice the actors, same distribution
#' fmh2.template <- template_network(fmh.ego, N=2*network.size(faux.mesa.high))
#' summary(fmh2.template, print.adj = FALSE)
#'
#' @import network
#' @export
template_network<-function(x, N, scaling=c("round","sample"), ...){
  scaling <- match.arg(scaling)
  w <- NVL(weights(x), rep(1,nrow(x)))
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

  x <- as_tibble(x$ego)[egoinds,]
  
  for(ego.col in names(x))
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
na.omit.egor <- function(object, relevant=TRUE, ...){
  # Create a subdataset containing only the relevant variables.
  obj <- object[,relevant,unit="ego"][,relevant,unit="alter"]

  # What this deos: for each row, for each column, sees if any element
  # (including in the alters table) is an NA.
  na.egos <- apply(obj, 1, function(r) any(sapply(lapply(r, is.na),any)))
  object[!na.egos,]
}

# Not really a generic function, but perhaps should be.

#' Draw random egocentric subsamples
#'
#' Implementations of the [base::sample()] function for [egor::egor()] data.
#'
#' @param x,size,replace,prob see [base::sample()].
#' @param ... extra arguments, currently unused.
#'
#' @note A reimplementation of sample as a generic was necessary
#'   because [base::sample()] is not a generic and cannot take
#'   data-frame-alikes as arguments.
#'
#' @return An [egor::egor()] object whose egos have been resampled in
#'   accordance with the arguments. Note that its [egor::ego_design()]
#'   information is overwritten in favor of the selection
#'   probabilities used in the sampling.
#'
#' @examples
#'
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#'
#' # Create a tiny weighted sample:
#' (s3 <- sample(fmh.ego, 3, replace=TRUE, prob=1:nrow(fmh.ego$ego)))
#' # Resampling with prob=weights(egor) creates a self-weighted
#' # sample:
#' (sample(s3, 3, replace=TRUE, prob=weights(s3)))
#'
#' # Create a large weighted sample, oversampling 12th-graders:
#' p <- ifelse(as_tibble(fmh.ego$ego)$Grade==12, 2, 1)
#' s2000 <- sample(fmh.ego, 2000, replace=TRUE, prob=p)
#'
#' # Summary function adjusts for weights:
#' (summ.net <- summary(faux.mesa.high ~ edges + nodematch("Grade") +
#'                      nodefactor("Race") + gwesp(0,fix=TRUE)))
#' (summ.ego <- summary(s2000 ~ edges + nodematch("Grade") + 
#'                      nodefactor("Race") + gwesp(0,fix=TRUE),
#'                      scaleto=network.size(faux.mesa.high)))
#'
#' \dontshow{
#' stopifnot(isTRUE(all.equal(
#'   as.vector(summ.net),
#'   as.vector(summ.ego),
#'   tolerance=.05,
#'   check.attributes=FALSE
#' )))
#' }
#' 
#' @export
sample <- function(x, size, replace=FALSE, prob=NULL, ...) UseMethod("sample")
#' @rdname sample
#' @export
sample.default <- function(x, ...) base::sample(x, ...)

#' @rdname sample
#' @export
sample.egor <- function(x, size, replace=FALSE, prob=NULL, ...){
  N <- nrow(x$ego)
  if(missing(size)) size <- N

  w <- weights(x)
  if(is.null(prob)) prob <- rep(size/N, N)
  is <- sample.int(N, size, replace, prob)

  x <- x[is, ,unit="ego"]
  x$ego$.sample_weights <- (w/prob)[is]
  ego_design(x) <- list(weights=".sample_weights")
  x
}

