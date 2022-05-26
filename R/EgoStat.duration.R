#  File R/EgoStat.duration.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################
EgoStat.mean.age <- function(egor, emptyval=0){
  startcol <- attr(egor,"alter_design")$startcol
  if(is.null(startcol)) stop("Egocentric dataset does not appear to contain durational information.")

  nedges <- nrow(egor$alter)
  if(nedges==0){
    out <- emptyval
  }else{
    w <- rep(weights(egor), .degreeseq(egor))
    starts <- egor$alter[[startcol]]
    out <- sum(w*starts)/sum(w)
  }    

  names(out) <- "mean.age"
  attr(out, "order") <- 0
  out
}

