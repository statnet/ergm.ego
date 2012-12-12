# An EgoStat.* function takes the data frame of egos, a data frame of
# alters, and the arguments passed to the corresponding ERGM terms,
# and returns a matrix of h(e[i]) values, with egos in rows and
# elements of h(e[i]) in columns.

EgoStat.edges <- function(egodata){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
   

  ties<-merge(egos[egoIDcol],alters[egoIDcol],by=egoIDcol,suffixes=c(".ego",".alter"))

  alterct <- as.data.frame(table(ties[[egoIDcol]]))
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[egoIDcol],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  h <- cbind(egos$.degree)
  colnames(h) <- "edges"
  rownames(h) <- egos[[egoIDcol]]
  
  h[order(rownames(h)),,drop=FALSE]/2
}

EgoStat.nodecov <- function(egodata, attrname){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
   
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],ties[[egoIDcol]],isolates),x=c(ties$.e,ties$.a,rep(0,length(isolates))))
  
  h <- cbind(sapply(tapply(ties$x,list(egoID=ties$egoID),FUN=sum),identity))
  colnames(h) <- paste("nodecov",attrname,sep=".")
  
  h[order(rownames(h)),,drop=FALSE]/2
}


EgoStat.nodefactor <- function(egodata, attrname, base=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],ties[[egoIDcol]],isolates),x=c(ties$.e,ties$.a,rep(0,length(isolates))))
  
  h <- t(sapply(tapply(ties$x, list(egoID=ties$egoID), FUN=tabulate, nbins=length(levs)),identity))
  colnames(h) <- paste("nodefactor",attrname,levs,sep=".")  

  if(length(base)==0 || base==0) h[order(rownames(h)),,drop=FALSE]/2
  else h[order(rownames(h)),-base,drop=FALSE]/2
}

EgoStat.nodematch <- function(egodata, attrname, diff=FALSE, keep=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  ties$match <- ifelse(ties$.e==ties$.a, as.integer(ties$.e), 0)
  if(!is.null(keep)) ties$match[!(ties$match%in%keep)] <- 0
  if(!diff) ties$match[ties$match!=0] <- 1
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 

  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), match=c(ties$match,rep(0,length(isolates))))

  h <- t(rbind(sapply(tapply(ties$match, list(egoID=ties$egoID), FUN=tabulate, nbins=if(diff) length(levs) else 1),identity)))

  colnames(h) <- if(diff) paste("nodematch",attrname,levs,sep=".") else paste("nodematch",attrname,sep=".")
  
  h[order(rownames(h)),if(!is.null(keep) && diff) keep else TRUE,drop=FALSE]/2
}

EgoStat.absdiff <- function(egodata, attrname, pow=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties)<-c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties$.absdiff <- abs(ties$.e-ties$.a)^pow 
  
  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), absdiff=c(ties$.absdiff,rep(0,length(isolates))))
  
  h <- cbind(sapply(tapply(ties$absdiff,list(egoID=ties$egoID),FUN=sum),identity))         
  colnames(h) <- if(pow==1) paste("absdiff",attrname,sep=".") else paste("absdiff",pow,".",attrname,sep="")
  
  h[order(rownames(h)),,drop=FALSE]/2
}

EgoStat.degree <- function(egodata, d, by=NULL, homophily=FALSE){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],alters[[by]])))
  }
  
  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",".a")
  if(!is.null(by) && homophily) ties <- ties[ties$.e==ties$.a,]
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]))
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by) && !homophily){
    bys <- rep(levs,each=length(d))
    degs <- rep(d,length(levs))
    
    h <- sapply(seq_along(bys), function(i) egos$.degree==degs[i] & egos[[by]]==bys[i])
    colnames(h) <- paste("deg",degs,".",by,bys,sep="")
  }else{
    h <- sapply(d, function(i) egos$.degree==i)
    colnames(h) <- if(homophily) paste("degree",d,".homophily.",by,sep="") else paste("degree",d,sep="")
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[order(rownames(h)),,drop=FALSE]
}
