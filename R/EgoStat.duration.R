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
   
    ties<-merge(egos[c(egoIDcol,egodata$egoWt)],alters[c(egoIDcol,startcol)],by=egoIDcol,suffixes=c(".ego",".alter"))
    names(ties) <- c(egoIDcol, "w", "a")
    total.ages <- sum(ties$a*ties$w)/sum(ties$w)
    out <- total.ages/nedges
  }    

  names(out) <- "mean.age"
  attr(out, "nonscalable") <- TRUE
  out
}

