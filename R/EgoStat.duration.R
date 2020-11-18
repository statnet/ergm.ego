#  File R/EgoStat.duration.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################
EgoStat.mean.age <- function(egodata, emptyval=0){
  startcol <- egodata$startcol
  if(is.null(startcol)) stop("Egocentric dataset does not appear to contain durational information.")

  nedges <- summary(egodata~edges, individual=FALSE, scaleto=1)
  if(nedges==0){
    out <- emptyval
  }else{
    egos <- egodata$egos
    alters <- egodata$alters
    egoIDcol <- egodata$egoIDcol
   
    ties<-merge(egos[c(egoIDcol,"egoWt")],alters[c(egoIDcol,startcol)],by=egoIDcol,suffixes=c(".ego",".alter"))
    names(ties) <- c(egoIDcol, "w", "a")
    out <- sum(ties$a*ties$w)/sum(ties$w)
  }    

  names(out) <- "mean.age"
  attr(out, "nonscaling") <- TRUE
  out
}

