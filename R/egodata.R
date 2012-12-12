as.egodata <- function(object, ..., egoIDcol="egoID"){
  UseMethod("as.egodata")
}

as.egodata.data.frame <- function(object, alters, egoWt = 1, ..., egoIDcol="egoID"){
  egoWt <- rep(egoWt, length.out=nrow(object))
  out <- list(egos=object[order(object[[egoIDcol]]),], alters=.prune.alters(object, alters, egoIDcol), egoWt = egoWt[order(object[[egoIDcol]])], egoIDcol=egoIDcol)  
  class(out) <- "egodata"
  out
}

# Conduct an egocentric census from the undirected network y=,
# returning an egodata object. The corresponding vertex attributes of
# y= are copied into columns in these data frames, excluding
# attributes listed in special.cols=.
as.egodata.network<-function(object,special.cols=c("na","vertex.names"),...,egoIDcol="vertex.names"){
  N<-network.size(object)

  egoIDs<-object%v%egoIDcol

  egos<-list()
  egos[[egoIDcol]]<-egoIDs
  
  for(a in list.vertex.attributes(object))
    if(!(a %in% special.cols)) egos[[a]]<-get.vertex.attribute(object,attrname=a)

  el<-as.edgelist(object)
  el<-rbind(el,el[,2:1])
  alterS<-unlist(tapply(el[,2],INDEX=el[,1],FUN=c,simplify=FALSE))
  alter.eID<-unlist(tapply(el[,1],INDEX=el[,1],FUN=c,simplify=FALSE))
  
  alters<-list()

  alters[[egoIDcol]]<-alter.eID
    
  for(a in list.vertex.attributes(object))
    if(!(a %in% special.cols)) alters[[a]]<-get.vertex.attribute(object,attrname=a)[alterS]

  as.egodata(as.data.frame(egos,stringsAsFactors=FALSE),alters=as.data.frame(alters,stringsAsFactors=FALSE), egoIDcol=egoIDcol)
}

.prune.alters <- function(egos, alters, egoIDcol){
  eis <- egos[[egoIDcol]]
  aeis <- alters[[egoIDcol]]

  todel <- !(aeis %in% eis)

  if(any(todel)) alters[!todel,]
  else alters
}

as.network.egodata<-function(x, N, scaling=c("greedy","sample"), ...){
  y0<-network.initialize(N,directed=FALSE)
  egos <- x$egos
  
  scaling <- match.arg(scaling)
  egoinds <- switch(scaling,
                    greedy={
                      .greedy.scaling(N,x$egoWt)
                    },
                    sample={
                      sample(length(x$egoWt),N,TRUE,x$egoWt)
                    })

  egos <- egos[egoinds,]
  
  for(ego.col in names(egos))
    if(is.factor(egos[[ego.col]]))
      y0 <- set.vertex.attribute(y0,ego.col,as.character(egos[[ego.col]]))
    else
      y0 <- set.vertex.attribute(y0,ego.col,egos[[ego.col]])
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
