summary.statistics.egodata <- function(object,..., basis=NULL, individual=FALSE, scaleto=NULL){
  egodata <-
    if(!is.null(basis)) basis
    else get(as.character(object[[2]]), envir=environment(object))

  stats <- NULL
  
  for(trm in term.list.formula(object[[length(object)]])){
    if(is.call(trm)){
      init.call <- list(as.name(paste("EgoStat.", trm[[1]],sep="")),egodata=egodata)
      init.call <- c(init.call,as.list(trm[-1]))
    }else{
      init.call <- list(as.name(paste("EgoStat.", trm,sep="")),egodata=egodata)
    }
    stat<-eval(as.call(init.call), environment(object))
    stats<-cbind(stats,stat)
  }
  rownames(stats) <- egodata$egos[[egodata$egoIDcol]]

  if(!individual){
    scaleto <- if(is.null(scaleto)) nrow(egodata$egos) else scaleto
    stats <- colSums(stats*egodata$egoWt)/sum(rep(egodata$egoWt,length.out=nrow(stats)))
    stats*scaleto
  }else stats
}
