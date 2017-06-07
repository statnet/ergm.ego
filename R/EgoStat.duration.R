#  File R/EgoStat.duration.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
#' @export
#' @rdname ergm.ego-terms
EgoStat.mean.age <- function(egor, emptyval=0){
  startcol <- attr(egor,"alter.design")$startcol
  if(is.null(startcol)) stop("Egocentric dataset does not appear to contain durational information.")

  nedges <- summary(egor~edges, individual=FALSE, scaleto=1)
  if(nedges==0){
    out <- emptyval
  }else{
    d <- sapply(egor, nrow)
    w <- rep(weights(egor), d)
    starts <- .allAlters(egor)[[startcol]]
    out <- sum(w*starts)/sum(w)
  }    

  names(out) <- "mean.age"
  attr(out, "order") <- 0
  out
}

