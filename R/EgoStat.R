#  File R/EgoStat.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
# An EgoStat.* function takes the data frame of egos, a data frame of
# alters, and the arguments passed to the corresponding ERGM terms,
# and returns a matrix of h(e[i]) values, with egos in rows and
# elements of h(e[i]) in columns.

#' \code{\link[ergm]{ergm}} Terms Implemented for
#' \code{\link[=egodata.object]{egodata}}
#' 
#' This page describes the \code{\link[ergm]{ergm}} terms (and hence network
#' statistics) for which inference based on egocentrically sampled data is
#' implemented in \code{ergm.ego} package. Other packages may add their own
#' terms.
#' 
#' The current recommendation for any package implementing additional
#' egocentric calculator terms is to create a help file with a name or alias
#' \code{ergm.egodata-terms}, so that \code{help("ergm.egodata-terms")} will
#' list egocentric ERGM terms available from all loaded packages.
#' 
#' 
#' @name ergm.ego-terms
#' @aliases ergm.ego-terms terms-ergm.ego ergm.ego.terms terms.ergm.ego
#' ergm-terms ergm.terms terms-ergm terms.ergm EgoStat EgoStat.edges
#' EgoStat.nodecov EgoStat.nodefactor EgoStat.nodematch EgoStat.nodemix
#' EgoStat.absdiff EgoStat.degree EgoStat.degrange EgoStat.concurrent
#' EgoStat.concurrentties EgoStat.degreepopularity EgoStat.mean.age netsize.adj
#' InitErgmTerm.netsize.adj
#' @docType methods
#' @section Currently implemented egocentric statistics: For each of these,
#' please see their respective package's \code{ergm-terms} help for meaning and
#' parameters. The simplest way to do this is usually via \code{? TERM}.
#' 
#' \describe{ \item{Special-purpose terms:}{ \describe{
#' \item{netsize.adj}{A special-purpose term equivalent to
#' \code{\link[ergm]{edges}}, to house the network-size adjustment offset. This
#' term is added to the model automatically and should not be used in the model
#' formula directly.  } } }
#' 
#' \item{ergm:}{ \itemize{ \item \code{edges} \item \code{nodecov}
#' \item \code{nodefactor} \item \code{nodematch} \item \code{nodemix} \item
#' \code{absdiff} \item \code{degree} \item \code{degrange} \item
#' \code{concurrent} \item \code{concurrentties} \item \code{degreepopularity}
#' } }
#' 
#' \item{tergm:}{ \itemize{ \item \code{mean.age} } } }
#' @seealso \code{\link[ergm]{ergm-terms}}
#' @keywords models
NULL


#' @export
EgoStat.edges <- function(egodata){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
   

  ties<-merge(egos[egoIDcol],alters[egoIDcol],by=egoIDcol,suffixes=c(".ego",".alter"))

  alterct <- as.data.frame(table(ties[[egoIDcol]]), stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[egoIDcol],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  h <- cbind(egos$.degree)
  colnames(h) <- "edges"
  rownames(h) <- egos[[egoIDcol]]

  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
}


#' @export
EgoStat.nodecov <- function(egodata, attrname){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  alt <- !is.null(alters[[attrname]])
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,if(alt) attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],if(alt) ties[[egoIDcol]],isolates),x=c(ties$.e,if(alt) ties$.a,rep(0,length(isolates))),stringsAsFactors=FALSE)
  
  h <- cbind(sapply(tapply(ties$x,list(egoID=ties$egoID),FUN=sum),identity)) / if(alt) 2 else 1
  colnames(h) <- paste("nodecov",attrname,sep=".")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}


#' @export
EgoStat.nodefactor <- function(egodata, attrname, base=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  # If there are multiple attributes, concatenate their names with a
  # dot and concatenate their values with a dot.
  if(length(attrname)>1){
    attrnamename <- paste(attrname, collapse=".")
    egos[[attrnamename]] <- do.call(paste,c(as.list(egos[,attrname]),list(sep=".")))
    alters[[attrnamename]] <- do.call(paste,c(as.list(alters[,attrname]),list(sep=".")))
    attrname <- attrnamename
  }
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)

  alt <- !is.null(alters[[attrname]])
  
  if(alt) alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,if(alt) attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],if(alt) ties[[egoIDcol]],isolates),x=c(ties$.e,if(alt) ties$.a,rep(0,length(isolates))),stringsAsFactors=FALSE)

  h <- t(sapply(tapply(ties$x, list(egoID=ties$egoID), FUN=tabulate, nbins=length(levs)),identity)) / if(alt) 2 else 1
  colnames(h) <- paste("nodefactor",attrname,levs,sep=".")  

  if(length(base)==0 || base==0) h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
  else h[match(egodata$egos[[egoIDcol]],rownames(h)),-base,drop=FALSE]
}

#' @export
EgoStat.nodematch <- function(egodata, attrname, diff=FALSE, keep=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  # If there are multiple attributes, concatenate their names with a
  # dot and concatenate their values with a dot.
  if(length(attrname)>1){
    attrnamename <- paste(attrname, collapse=".")
    egos[[attrnamename]] <- do.call(paste,c(as.list(egos[,attrname]),list(sep=".")))
    alters[[attrnamename]] <- do.call(paste,c(as.list(alters[,attrname]),list(sep=".")))
    attrname <- attrnamename
  }
  
  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  ties$match <- ifelse(ties$.e==ties$.a, as.integer(ties$.e), 0)
  if(!is.null(keep)) ties$match[!(ties$match%in%keep)] <- 0
  if(!diff) ties$match[ties$match!=0] <- 1
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 

  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), match=c(ties$match,rep(0,length(isolates))),stringsAsFactors=FALSE)

  h <- t(rbind(sapply(tapply(ties$match, list(egoID=ties$egoID), FUN=tabulate, nbins=if(diff) length(levs) else 1),identity)))

  colnames(h) <- if(diff) paste("nodematch",attrname,levs,sep=".") else paste("nodematch",attrname,sep=".")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),if(!is.null(keep) && diff) keep else TRUE,drop=FALSE]/2
}


#' @export
EgoStat.nodemix <- function(egodata, attrname, base=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  # If there are multiple attributes, concatenate their names with a
  # dot and concatenate their values with a dot.
  if(length(attrname)>1){
    attrnamename <- paste(attrname, collapse=".")
    egos[[attrnamename]] <- do.call(paste,c(as.list(egos[,attrname]),list(sep=".")))
    alters[[attrnamename]] <- do.call(paste,c(as.list(alters[,attrname]),list(sep=".")))
    attrname <- attrnamename
  }

  levs <- sort(unique(c(egos[[attrname]],alters[[attrname]])))
  egos[[attrname]] <- match(egos[[attrname]], levs, 0)
  alters[[attrname]] <- match(alters[[attrname]], levs, 0)
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  
  
  namevec <- outer(levs,levs,paste,sep=".")
  namevec <- namevec[upper.tri(namevec,diag=TRUE)]
  
  if (!is.null(base) && !identical(base,0)) {
    namevec <- namevec[-base]
  }
  
  mat <- t(apply(cbind(ties[,2:3]),1,sort))
  h <- table(ties[,1],paste(levs[mat[,1]],levs[mat[,2]],sep="."))
  
  h<- t(apply(h,1,function(x)x[namevec,drop=FALSE]))
  
  if(length(isolates)){
    isolates.mat <- matrix(0,nrow=length(isolates),ncol=length(namevec))
    rownames(isolates.mat) <- isolates	
    h <- rbind(h,isolates.mat)
  }
  
  h <- h[order(as.numeric(rownames(h))),]
  h[is.na(h)] <- 0
  colnames(h) <- paste("mix",attrname,namevec,sep=".")
  h[match(egodata$egos[[egoIDcol]],rownames(h)),]/2
}

#' @export
EgoStat.absdiff <- function(egodata, attrname, pow=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties)<-c(egoIDcol,".e",".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties$.absdiff <- abs(ties$.e-ties$.a)^pow 
  
  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), absdiff=c(ties$.absdiff,rep(0,length(isolates))),stringsAsFactors=FALSE)
  
  h <- cbind(sapply(tapply(ties$absdiff,list(egoID=ties$egoID),FUN=sum),identity))         
  colnames(h) <- if(pow==1) paste("absdiff",attrname,sep=".") else paste("absdiff",pow,".",attrname,sep="")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
}

#' @export
EgoStat.degree <- function(egodata, d, by=NULL, homophily=FALSE){
  ## if(any(d==0)) warning("degree(0) (isolate) count statistic depends strongly on the specified population network size.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  alt <- !is.null(by) && !is.null(alters[[by]])
  if(homophily && !alt) stop("Attribute ", sQuote(by), " must be observed on alters if homophily=TRUE.")
  
  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],if(alt) alters[[by]])))
  }
  
  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,if(alt) by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e", if(alt) ".a")
  if(!is.null(by) && homophily) ties <- ties[ties$.e==ties$.a,]
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
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
    colnames(h) <- if(homophily) paste("deg",d,".homophily.",by,sep="") else paste("degree",d,sep="")
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

#' @export
EgoStat.degrange <- function(egodata, from=NULL, to=Inf, by=NULL, homophily=FALSE){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  to <- ifelse(to==Inf, .Machine$integer.max, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  alt <- !is.null(by) && !is.null(alters[[by]])
  if(homophily && !alt) stop("Attribute ", sQuote(by), " must be observed on alters if homophily=TRUE.")
  
  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],if(alt) alters[[by]])))
  }

  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,if(alt) by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  if(!is.null(by) && homophily) ties <- ties[ties$.e==ties$.a,]
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by) && !homophily){
    bys <- rep(levs,each=length(from))
    froms <- rep(from,length(levs))
    tos <- rep(to,length(levs))
    
    h <- sapply(seq_along(bys), function(i) egos$.degree>=froms[i] & egos$.degree<tos[i] & egos[[by]]==bys[i])
    colnames(h) <-  ifelse(tos>=.Machine$integer.max,
                           paste("deg", from, "+.",          by, bys, sep=""),
                           paste("deg", from, "to", to, ".", by, bys, sep=""))

  }else{
    h <- sapply(seq_along(from), function(i) egos$.degree>=from[i] & egos$.degree<to[i])
    colnames(h) <-
      if(homophily)
        ifelse(to>=.Machine$integer.max,
               paste("deg", from,  "+",     ".homophily.", by, sep=""),
               paste("deg", from, "to", to, ".homophily.", by, sep=""))
      else
        ifelse(to>=.Machine$integer.max,
               paste("deg", from,  "+", sep=""),
               paste("deg", from, "to", to, sep=""))

  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

#' @export
EgoStat.concurrent <- function(egodata, by=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  alt <- !is.null(by) && !is.null(alters[[by]])
   
  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],if(alt) alters[[by]])))
  }

  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,if(alt) by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by)){
    bys <- levs
    
    h <- sapply(seq_along(bys), function(i) egos$.degree>=2 & egos[[by]]==bys[i])
    colnames(h) <- paste("concurrent.", by, bys, sep="")

  }else{
    h <- cbind(egos$.degree>=2)
    colnames(h) <- "concurrent"
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

#' @export
EgoStat.concurrentties <- function(egodata, by=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  alt <- !is.null(by) && !is.null(alters[[by]])
  
  if(!is.null(by)){
    levs <- sort(unique(c(egos[[by]],if(alt) alters[[by]])))
  }

  ties<-merge(egos[c(egoIDcol,by)],alters[c(egoIDcol,if(alt) by)],by=egoIDcol,suffixes=c(".ego",".alter"))

  if(!is.null(by)) names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  ties$.a <- NULL

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol,by)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  if(!is.null(by)){
    bys <- levs
    
    h <- sapply(seq_along(bys), function(i) cbind(ifelse(egos[[by]]==bys[i], pmax(egos$.degree-1,0), 0)))
    colnames(h) <- paste("concurrentties.", by, bys, sep="")

  }else{
    h <- cbind(pmax(egos$.degree-1,0))
    colnames(h) <- "concurrentties"    
  }
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}


#' @export
EgoStat.degreepopularity <- function(egodata){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
    
  ties<-merge(egos[egoIDcol],alters[egoIDcol],by=egoIDcol,suffixes=c(".ego",".alter"))

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  h <- cbind(egos$.degree^(3/2))
  colnames(h) <- "degreepopularity"
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}
