#  File R/EgoStat.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2019 Statnet Commons
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
#' 
#' @aliases ergm.ego-terms terms-ergm.ego ergm.ego.terms
#'   terms.ergm.ego ergm-terms ergm.terms terms-ergm terms.ergm
#'   EgoStat EgoStat.offset EgoStat.edges EgoStat.nodecov
#'   EgoStat.nodefactor EgoStat.nodematch EgoStat.nodemix
#'   EgoStat.absdiff EgoStat.degree EgoStat.degrange
#'   EgoStat.concurrent EgoStat.concurrentties EgoStat.degree1.5
#'   EgoStat.mm EgoStat.mean.age netsize.adj InitErgmTerm.netsize.adj
#'
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
#' \item{ergm:}{ \itemize{ \item \code{offset} \item \code{edges}
#' \item \code{nodecov} \item \code{nodefactor} \item \code{nodematch}
#' \item \code{nodemix} \item \code{absdiff} \item \code{degree} \item
#' \code{degrange} \item \code{concurrent} \item \code{concurrentties}
#' \item \code{degree1.5} \item `mm` } }
#' 
#' \item{tergm:}{ \itemize{ \item \code{mean.age} } } }
#' @seealso \code{\link[ergm]{ergm-terms}}
#' @keywords models
NULL

# copied from ergm
LEVELS_BASE1 <- NULL

#' @export
EgoStat.offset <- function(egodata, trm){
  trm <- substitute(trm)
  if(is.call(trm)){
    init.call <- list(as.name(paste("EgoStat.", trm[[1]],sep="")),egodata=egodata)
    init.call <- c(init.call,as.list(trm[-1]))
  }else{
    init.call <- list(as.name(paste("EgoStat.", trm,sep="")),egodata=egodata)
  }
  h <- eval(as.call(init.call))
  
  # needed to match ergm's naming convention for offsets  
  colnames(h) <- paste0("offset(", colnames(h), ")")
  h
}

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
EgoStat.nodecov <- function(egodata, attr){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  xe <- ergm.ego_get_vattr(attr, egos, accept = "numeric", multiple = "matrix")
  xa <- ERRVL(try(ergm.ego_get_vattr(attr, alters, accept = "numeric", multiple = "matrix"), silent=TRUE), NULL)
  
  attrnames <- if(is.matrix(xe)) colnames(xe) else attributes(xe)$name
  alt <- !is.null(xa)

  egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], attrnames = xe, stringsAsFactors=FALSE), c(egoIDcol, attrnames))
  alters <- NVL2(xa, setNames(data.frame(egoIDcol = alters[[egoIDcol]], attrnames = xa, stringsAsFactors=FALSE), c(egoIDcol, attrnames)),
                     setNames(data.frame(egoIDcol = alters[[egoIDcol]], stringsAsFactors=FALSE), egoIDcol))
  
  ties.all <-merge(egos[c(egoIDcol,attrnames)],alters[c(egoIDcol,if(alt) attrnames)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties.all) <- c(egoIDcol, paste0(attrnames, ".ego"), if(alt) paste0(attrnames, ".alter"))
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties.all[[egoIDcol]])] 

  h.all <- NULL
  for(attrname in attrnames) {
    ties <- ties.all[c(egoIDcol, paste0(attrname, c(".ego", if(alt) ".alter")))]
    names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  
    ties <- data.frame(egoID=c(ties[[egoIDcol]],if(alt) ties[[egoIDcol]],isolates),x=c(ties$.e,if(alt) ties$.a,rep(0,length(isolates))),stringsAsFactors=FALSE)
  
    h <- cbind(sapply(tapply(ties$x,list(egoID=ties$egoID),FUN=sum),identity)) / if(alt) 2 else 1
    
    colname <- "nodecov"
    if(is.matrix(xe)) colname <- paste(colname, attributes(xe)$name, sep = ".")
    colname <- paste(colname, attrname, sep = ".")
    colnames(h) <- colname
  
    h.all <- cbind(h.all, h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE])
  }

  h.all
}


#' @export
EgoStat.nodefactor <- function(egodata, attr, base=1, levels=LEVELS_BASE1){
  if(!missing(base)) message("In term `nodefactor' in package `ergm.ego': Argument \"base\" has been superseded by \"levels\" and it is recommended to use the latter.  Note that its interpretation may be different.")

  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  xe <- ergm.ego_get_vattr(attr, egos)
  xa <- ERRVL(try(ergm.ego_get_vattr(attr, alters), silent=TRUE), NULL)

  attrname <- attributes(xe)$name
  alt <- !is.null(xa)
  
  levs <- ergm.ego_attr_levels(levels, c(xe, xa), egodata, sort(unique(c(xe, xa))))

  xe <- match(xe, levs, 0)
  if(alt) xa <- match(xa, levs, 0)

  egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], attrname = xe, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  alters <- NVL2(xa, setNames(data.frame(egoIDcol = alters[[egoIDcol]], attrname = xa, stringsAsFactors=FALSE), c(egoIDcol, attrname)),
                     setNames(data.frame(egoIDcol = alters[[egoIDcol]], stringsAsFactors=FALSE), egoIDcol))
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,if(alt) attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",if(alt) ".a")
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  ties <- data.frame(egoID=c(ties[[egoIDcol]],if(alt) ties[[egoIDcol]],isolates),x=c(ties$.e,if(alt) ties$.a,rep(0,length(isolates))),stringsAsFactors=FALSE)

  h <- t(sapply(tapply(ties$x, list(egoID=ties$egoID), FUN=tabulate, nbins=length(levs)),identity)) / if(alt) 2 else 1
  
  # handles case of one level
  if(length(levs) == 1) h <- t(h)
  
  colnames(h) <- paste("nodefactor",attrname,levs,sep=".")  

  if(length(base)==0 || base==0 || !missing(levels)) h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
  else h[match(egodata$egos[[egoIDcol]],rownames(h)),-base,drop=FALSE]
}

#' @export
EgoStat.nodematch <- function(egodata, attr, diff=FALSE, keep=NULL, levels=NULL){
  if(!missing(keep)) message("In term `nodematch' in package `ergm.ego': Argument \"keep\" has been superseded by \"levels\" and it is recommended to use the latter.  Note that its interpretation may be different.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  xe <- ergm.ego_get_vattr(attr, egos)
  xa <- ergm.ego_get_vattr(attr, alters)

  attrname <- attributes(xe)$name
  
  levs <- ergm.ego_attr_levels(levels, c(xe, xa), egodata, sort(unique(c(xe, xa))))

  xe <- match(xe, levs, 0)
  xa <- match(xa, levs, 0)

  egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], attrname = xe, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  alters <- setNames(data.frame(egoIDcol = alters[[egoIDcol]], attrname = xa, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  names(ties) <- c(egoIDcol,".e",".a")
  ties$match <- ifelse(ties$.e==ties$.a, as.integer(ties$.e), 0)
  if(!is.null(keep) && missing(levels)) ties$match[!(ties$match%in%keep)] <- 0
  if(!diff) ties$match[ties$match!=0] <- 1
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 

  ties <- data.frame(egoID=c(ties[[egoIDcol]],isolates), match=c(ties$match,rep(0,length(isolates))),stringsAsFactors=FALSE)

  h <- t(rbind(sapply(tapply(ties$match, list(egoID=ties$egoID), FUN=tabulate, nbins=if(diff) length(levs) else 1),identity)))

  colnames(h) <- if(diff) paste("nodematch",attrname,levs,sep=".") else paste("nodematch",attrname,sep=".")
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),if(!is.null(keep) && diff) keep else TRUE,drop=FALSE]/2
}


#' @export
EgoStat.nodemix <- function(egodata, attr, base=NULL, levels=NULL, levels2=NULL){
  if(!missing(base)) message("In term `nodemix' in package `ergm.ego': Argument \"base\" has been superseded by \"levels2\" and it is recommended to use the latter.  Note that its interpretation may be different.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  xeval <- ergm.ego_get_vattr(attr, egos)
  xaval <- ergm.ego_get_vattr(attr, alters)

  attrname <- attributes(xeval)$name
  
  levs <- ergm.ego_attr_levels(levels, c(xeval, xaval), egodata, sort(unique(c(xeval, xaval))))

  xe <- match(xeval, levs, 0)
  xa <- match(xaval, levs, 0)

  egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], attrname = xe, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  alters <- setNames(data.frame(egoIDcol = alters[[egoIDcol]], attrname = xa, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  
  ties<-merge(egos[c(egoIDcol,attrname)],alters[c(egoIDcol,attrname)],by=egoIDcol,suffixes=c(".ego",".alter"))
  
  # only consider ties where both ego and alter have attr value included in levs
  # others do not contribute to the statistics we are generating
  ties<-ties[ties[,2] > 0 & ties[,3] > 0,]
  
  names(ties) <- c(egoIDcol,".e",".a")
  
  isolates <- egos[[egoIDcol]][!(egos[[egoIDcol]]%in%ties[[egoIDcol]])] 
  
  nr <- length(levs)
  nc <- length(levs)

  levels2.list <- transpose(expand.grid(row = levs, col = levs, stringsAsFactors=FALSE))
  indices2.grid <- expand.grid(row = 1:nr, col = 1:nc)
  uun <- as.vector(outer(levs,levs,paste,sep="."))
    
  rowleqcol <- indices2.grid$row <= indices2.grid$col
  levels2.list <- levels2.list[rowleqcol]
  indices2.grid <- indices2.grid[rowleqcol,]
  uun <- uun[rowleqcol]
   
  levels2.sel <- ergm.ego_attr_levels(levels2, list(row = c(xeval, xaval), col = c(xaval, xeval)), egodata, levels2.list)

  if(!is.null(base) && !identical(base,0) && missing(levels2)) levels2.sel <- levels2.sel[-base]
    
  rows2keep <- match(levels2.sel,levels2.list, NA)
  rows2keep <- rows2keep[!is.na(rows2keep)]
  
  u <- indices2.grid[rows2keep,]
  namevec <- uun[rows2keep]
  
  mat <- t(apply(cbind(ties[,2:3]),1,sort))

  h <- table(ties[,1],paste(levs[mat[,1]],levs[mat[,2]],sep="."))
  
  h <- t(apply(h,1,function(x)x[namevec,drop=FALSE]))
  
  # handles case of single retained category
  if(length(namevec) == 1) h <- t(h)
  
  if(length(isolates)){
    isolates.mat <- matrix(0,nrow=length(isolates),ncol=length(namevec))
    rownames(isolates.mat) <- isolates	
    h <- rbind(h,isolates.mat)
  }
  
  h <- h[order(as.numeric(rownames(h))),,drop=FALSE]
  h[is.na(h)] <- 0
  colnames(h) <- paste("mix",attrname,namevec,sep=".")
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]/2
}

#' @export
EgoStat.absdiff <- function(egodata, attr, pow=1){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  xe <- ergm.ego_get_vattr(attr, egos, accept = "numeric")
  xa <- ergm.ego_get_vattr(attr, alters, accept = "numeric")
  
  attrname <- attributes(xe)$name

  egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], attrname = xe, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  alters <- setNames(data.frame(egoIDcol = alters[[egoIDcol]], attrname = xa, stringsAsFactors=FALSE), c(egoIDcol, attrname))
  
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
EgoStat.degree <- function(egodata, d, by=NULL, homophily=FALSE, levels=NULL){
  ## if(any(d==0)) warning("degree(0) (isolate) count statistic depends strongly on the specified population network size.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    xa <- ERRVL(try(ergm.ego_get_vattr(by, alters), silent=TRUE), NULL)
    
    by <- attributes(xe)$name
  
    egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], by = xe, stringsAsFactors=FALSE), c(egoIDcol, by))
    alters <- NVL2(xa, setNames(data.frame(egoIDcol = alters[[egoIDcol]], by = xa, stringsAsFactors=FALSE), c(egoIDcol, by)),
                       setNames(data.frame(egoIDcol = alters[[egoIDcol]], stringsAsFactors=FALSE), egoIDcol))

    levs <- ergm.ego_attr_levels(levels, c(xe, xa), egodata, sort(unique(c(xe, xa))))
  }

  alt <- !is.null(by) && !is.null(alters[[by]])
  if(homophily && !alt) stop("Attribute ", sQuote(by), " must be observed on alters if homophily=TRUE.")
    
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
EgoStat.degrange <- function(egodata, from=NULL, to=Inf, by=NULL, homophily=FALSE, levels=NULL){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
  
  to <- ifelse(to==Inf, .Machine$integer.max, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    xa <- ERRVL(try(ergm.ego_get_vattr(by, alters), silent=TRUE), NULL)
    
    by <- attributes(xe)$name
  
    egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], by = xe, stringsAsFactors=FALSE), c(egoIDcol, by))
    alters <- NVL2(xa, setNames(data.frame(egoIDcol = alters[[egoIDcol]], by = xa, stringsAsFactors=FALSE), c(egoIDcol, by)),
                       setNames(data.frame(egoIDcol = alters[[egoIDcol]], stringsAsFactors=FALSE), egoIDcol))

    levs <- ergm.ego_attr_levels(levels, c(xe, xa), egodata, sort(unique(c(xe, xa))))
  }

  alt <- !is.null(by) && !is.null(alters[[by]])
  if(homophily && !alt) stop("Attribute ", sQuote(by), " must be observed on alters if homophily=TRUE.")
  
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
EgoStat.concurrent <- function(egodata, by=NULL, levels=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    xa <- ERRVL(try(ergm.ego_get_vattr(by, alters), silent=TRUE), NULL)
    
    by <- attributes(xe)$name
  
    egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], by = xe, stringsAsFactors=FALSE), c(egoIDcol, by))
    alters <- NVL2(xa, setNames(data.frame(egoIDcol = alters[[egoIDcol]], by = xa, stringsAsFactors=FALSE), c(egoIDcol, by)),
                       setNames(data.frame(egoIDcol = alters[[egoIDcol]], stringsAsFactors=FALSE), egoIDcol))

    levs <- ergm.ego_attr_levels(levels, c(xe, xa), egodata, sort(unique(c(xe, xa))))
  }

  alt <- !is.null(by) && !is.null(alters[[by]])
   
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
EgoStat.concurrentties <- function(egodata, by=NULL, levels=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    xa <- ERRVL(try(ergm.ego_get_vattr(by, alters), silent=TRUE), NULL)
    
    by <- attributes(xe)$name
  
    egos <- setNames(data.frame(egoIDcol = egos[[egoIDcol]], by = xe, stringsAsFactors=FALSE), c(egoIDcol, by))
    alters <- NVL2(xa, setNames(data.frame(egoIDcol = alters[[egoIDcol]], by = xa, stringsAsFactors=FALSE), c(egoIDcol, by)),
                       setNames(data.frame(egoIDcol = alters[[egoIDcol]], stringsAsFactors=FALSE), egoIDcol))

    levs <- ergm.ego_attr_levels(levels, c(xe, xa), egodata, sort(unique(c(xe, xa))))
  }

  alt <- !is.null(by) && !is.null(alters[[by]])
  
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
EgoStat.degree1.5 <- function(egodata){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol
    
  ties<-merge(egos[egoIDcol],alters[egoIDcol],by=egoIDcol,suffixes=c(".ego",".alter"))

  alterct <- as.data.frame(table(ties[[egoIDcol]]),stringsAsFactors=FALSE)
  colnames(alterct)<-c(egoIDcol,".degree")

  egos <- merge(egos[c(egoIDcol)],alterct,by=egoIDcol,all=TRUE)
  egos$.degree[is.na(egos$.degree)]<-0

  h <- cbind(egos$.degree^(3/2))
  colnames(h) <- "degree1.5"
  rownames(h) <- egos[[egoIDcol]]
  
  h[match(egodata$egos[[egoIDcol]],rownames(h)),,drop=FALSE]
}

#' @export
EgoStat.mm <- function(egodata, attrs, levels=NULL, levels2=NULL){
  egos <- egodata$egos
  alters <- egodata$alters
  egoIDcol <- egodata$egoIDcol

  # Some preprocessing steps are the same, so run together:
  #' @import purrr
  #' @importFrom utils relist
  #' @importFrom methods is
  spec <-
    list(attrs = attrs, levels = levels) %>%
    map_if(~!is(., "formula"), ~call("~", .)) %>% # Embed into RHS of formula.
    map_if(~length(.)==2, ~call("~", .[[2]], .[[2]])) %>% # Convert ~X to X~X.
    map(as.list) %>% map(~.[-1]) %>% # Convert to list(X,X).
    map(set_names, c("row", "col")) %>% # Name elements rowspec and colspec.
    transpose() %>%
    unlist(recursive=FALSE) %>% # Convert into a flat list.
    map_if(~is.name(.)&&.==".", ~NULL) %>% # If it's just a dot, convert to NULL.
    map_if(~is.call(.)||(is.name(.)&&.!="."), ~as.formula(call("~", .))) %>% # If it's a call or a symbol, embed in formula.
    relist(skeleton=list(row=c(attrs=NA, levels=NA), col=c(attrs=NA, levels=NA))) %>% # Reconstruct list.
    transpose()

  if(is(attrs, "formula"))
    spec[["attrs"]] <- lapply(spec[["attrs"]], function(x){if(is(x,"formula")) environment(x) <- environment(attrs); x})
  if(is(levels, "formula"))
    spec[["levels"]] <- lapply(spec[["levels"]], function(x){if(is(x,"formula")) environment(x) <- environment(levels); x})
  spec <- transpose(spec)
  
  # Extract attribute values.
  attrval <-
    spec %>%
    imap(function(spec, whose){
      if(is.null(spec$attrs)){
        list(valcodes =
               rep(0L,
                   nrow(alters)*2
                   ),
             name = ".",
             levels = NA,
             levelcodes = 0,
             id = rep(merge(data.frame(i=egos[[egoIDcol]]),
                            data.frame(i=alters[[egoIDcol]]))$i,
                      2)
             )
      }else{
        xe <- ERRVL(ec <- try(ergm.ego_get_vattr(spec$attrs, egos), silent=TRUE), NULL)
        xa <- ERRVL(try(ergm.ego_get_vattr(spec$attrs, alters), silent=TRUE), NULL)
        name <- attr(NVL(xe,xa), "name")
        if(is.null(xe)&&is.null(xa)) stop(attr(ec, "condition"), call.=FALSE) # I.e., they were both errors. => propagate error message.
        xe <- NVL2(xe,
                   data.frame(i=egos[[egoIDcol]], xe=xe, stringsAsFactors=FALSE),
                   data.frame(i=egos[[egoIDcol]], stringsAsFactors=FALSE))
        xa <- NVL2(xa,
                   data.frame(i=alters[[egoIDcol]], xa=xa, stringsAsFactors=FALSE),
                   data.frame(i=alters[[egoIDcol]], stringsAsFactors=FALSE))
        xae <- merge(xe,xa)
        x <- switch(whose,
                    row = c(NVL(xae$xe,xae$xa),NVL(xae$xa,xae$xe)),
                    col = c(NVL(xae$xa,xae$xe),NVL(xae$xe,xae$xa)))
        list(name=name, id=rep(xae$i,length.out=length(x)), val=x, levels=spec$levels, unique=sort(unique(x)))
      }
    })

  # Undirected unipartite networks with identical attribute
  # specification produce square, symmetric mixing matrices. All
  # others do not.
  symm <- identical(spec$row$attrs, spec$col$attrs)
  # Are we evaluating the margin?
  marg <- length(attrval$row$unique)==0 || length(attrval$col$unique)==0
  
  # Filter the final level set and encode the attribute values.
  attrval <- attrval %>%
    map_if(~is.null(.$levelcodes), function(v){
      v$levels <- ergm.ego_attr_levels(v$levels, v$val, egodata, levels=v$unique)
      v$levelcodes <- seq_along(v$levels)
      v$valcodes <- match(v$val, v$levels, nomatch=0)
      v
    })

  # Construct all pairwise level combinations (table cells) and their numeric codes.
  levels2codes <- expand.grid(row=attrval$row$levelcodes, col=attrval$col$levelcodes) %>% transpose()
  levels2vals <- expand.grid(row=attrval$row$levels, col=attrval$col$levels, stringsAsFactors=FALSE) %>% transpose()

  # Drop redundant table cells if symmetrising.
  if(symm){
    levels2keep <- levels2codes %>% map_lgl(with, row <= col)
    levels2codes <- levels2codes[levels2keep]
    levels2vals <- levels2vals[levels2keep]
  }

  # Run the table cell list through the cell filter.
  levels2sel <- ergm.ego_attr_levels(levels2, list(row=attrval$row$val, col=attrval$col$val), egodata, levels=levels2vals)
  levels2codes <- levels2codes[match(levels2sel,levels2vals, NA)]
  levels2vals <- levels2sel; rm(levels2sel)

  # Construct the level names
  levels2names <-
    levels2vals %>%
    transpose() %>%
    map(unlist) %>%
    with(paste0(
      "[",
      if(length(attrval$row$levels)>1)
        paste0(attrval$row$name, "=", .$row)
      else ".",
      ",",
      if(length(attrval$col$levels)>1)
        paste0(attrval$col$name, "=", .$col)
      else ".",
      "]"))
  
  coef.names <- paste0("mm",levels2names)

  h <- attrval %>%
    map("valcodes") %>%
    transpose() %>%
    match(levels2codes) %>%
    map(tabulate, length(levels2codes)) %>%
    do.call(rbind,.) %>%
    aggregate(by=list(i=attrval$row$id), FUN=sum)

  i <- h$i
  h <- as.matrix(h)[,-1,drop=FALSE]
  colnames(h) <- coef.names

  if(symm){
    selff <- 1+map_lgl(levels2codes, all_identical)
    h <- sweep(h, 2, selff, `/`)
  }
 
  h <- h[match(egos[[egoIDcol]], i),,drop=FALSE]/2
  h[is.na(h)] <- 0
  h
}
