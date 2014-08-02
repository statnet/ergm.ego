as.egodata <- function(object, ..., egoIDcol="egoID"){
  UseMethod("as.egodata")
}

as.egodata.data.frame <- function(object, alters, egoWt = 1, ..., egoIDcol="egoID"){
  egoWt <- rep(egoWt, length.out=nrow(object))
  out <- list(egos=object, alters=.prune.alters(object, alters, egoIDcol), egoWt = egoWt, egoIDcol=egoIDcol)  
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
  if(any(duplicated(egoIDs))){
    warning("Non-unique ego IDs; using 1..N.")
    egoIDs <- seq_along(egoIDs)
  }
  
  egos<-list()
  egos[[egoIDcol]]<-egoIDs
  
  for(a in list.vertex.attributes(object))
    if(!(a %in% special.cols)) egos[[a]]<-get.vertex.attribute(object,attrname=a)

  el<-as.edgelist(object)
  el<-rbind(el,el[,2:1])
  alterS<-unlist(tapply(el[,2],INDEX=el[,1],FUN=c,simplify=FALSE))
  alter.eID<-egoIDs[unlist(tapply(el[,1],INDEX=el[,1],FUN=c,simplify=FALSE))]
  
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

  if(any(todel)) alters[!todel,,drop=FALSE]
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
  y0 %n% "ego.inds" <- egoinds
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

# Note: The following functions use parts of na.omit.data.frame() and
# subset.data.frame() from the R's stats package under the terms of
# the GNU GPL v3.

`[.egodata` <- function(x, i, j, ..., dup.action=c("make.unique", "fail", "number")){
  subset(x, i, j, ..., dup.action=dup.action)
}

subset.egodata <- function(x, subset, select, ..., dup.action=c("make.unique", "fail", "number")){
  if(missing(subset)) subset <- TRUE
  
  if (missing(select)) 
    egovars <- altervars <- TRUE
  else {
    if((is.numeric(select) || is.logical(select) || is.integer(select))
       && !isTRUE(all.equal(names(x$egos),names(x$alters)))){
      stop("Logical or numeric column subsets cannot be used if column orderings for egos and alters are not identical. Use column names instead.")
    }
    
    nl <- as.list(seq_along(x$egos))
    names(nl) <- names(x$egos)
    vars <- eval(substitute(select), nl, parent.frame())
    eIDi <- which(names(x$egos)==x$egoIDcol)
    egovars <- union(eIDi,vars)

    nl <- as.list(seq_along(x$alters))
    names(nl) <- names(x$alters)
    vars <- eval(substitute(select), nl, parent.frame())
    eIDi <- which(names(x$alters)==x$egoIDcol)
    altervars <- union(eIDi,vars)
  }

  egoIDs <- switch(mode(subset),
                   numeric =,
                   logical =,
                   integer = x$egos[[x$egoIDcol]][subset],
                   character = x$egos[[x$egoIDcol]][match(subset, x$egos[[x$egoIDcol]])])

  dup.action <- match.arg(dup.action)
  unique.egoIDs <- switch(dup.action,
                          fail=stop("Selected subset calls for duplicating egos."),
                          numeric=seq_along(egoIDs),
                          make.unique=make.unique(as.character(egoIDs)))
  
  egos <- cbind(x$egos[match(egoIDs,x$egos[[x$egoIDcol]]),egovars,drop=FALSE], .unique.egoIDs = unique.egoIDs, stringsAsFactors=FALSE)
  alters <- merge(egos[c(x$egoIDcol,".unique.egoIDs")], x$alters[altervars], by=x$egoIDcol)
  egoWt <- x$egoWt[match(egoIDs,x$egos[[x$egoIDcol]])]

  # If we have duplicated egoIDs, we have to handle it per dup.action.
  
  if(any(duplicated(egoIDs))){
    egos[[x$egoIDcol]] <- egos[[".unique.egoIDs"]]
    alters[[x$egoIDcol]] <- alters[[".unique.egoIDs"]]
  }
  egos[[".unique.egoIDs"]] <- NULL
  alters[[".unique.egoIDs"]] <- NULL

  out <- list(egos=egos, alters=alters, egoWt = egoWt, egoIDcol=x$egoIDcol)  
  class(out) <- "egodata"
  out
}

# TODO: A more efficient implementation of this.
na.omit.egodata <- function(object, relevant=TRUE, ...){
  # Create a subdataset containing only the relevant variables.
  obj <- subset(object,select=relevant)
  
  n <- length(obj$egos)
  omit <- FALSE
  vars <- seq_len(n)
  for (j in vars) {
    x <- obj$egos[[j]]
    if (!is.atomic(x)) 
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L) 
      omit <- omit | x
    else for (ii in 1L:d[2L]) omit <- omit | x[, ii]
  }
  ego.omit <- obj$egos[[obj$egoIDcol]][omit]

  n <- length(obj$alters)
  omit <- FALSE
  vars <- seq_len(n)
  for (j in vars) {
    x <- obj$alters[[j]]
    if (!is.atomic(x)) 
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L) 
      omit <- omit | x
    else for (ii in 1L:d[2L]) omit <- omit | x[, ii]
  }
  
  alter.omit <- obj$alters[[obj$egoIDcol]][omit]

  subset(object, -match(union(ego.omit,alter.omit),object$egos[[object$egoIDcol]]))  
}

dimnames.egodata <- function(x){
  list(x$egos[[x$egoIDcol]], names(x$egos), names(x$alters))
}

dim.egodata <- function(x){
  c(nrow(x$egos), ncol(x$egos), ncol(x$alters))
}

# Not really a generic function, but perhaps should be.

sample <- function(x, size, replace=FALSE, prob=NULL, ...) UseMethod("sample")
sample.default <- function(x, ...) base::sample(x, ...)

sample.egodata <- function(x, size, replace=FALSE, prob=NULL){
  if(missing(size)) size <- nrow(x)
  
  is <- sample.int(nrow(x), size, replace, prob)

  out <- subset(x, is)

  if(is.null(prob)) prob <- rep(1, nrow(x))
  
  out$egoWt <- x$egoWt[is]/prob[is]
  out
}
