#  File R/EgoStat.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2021 Statnet Commons
################################################################################
# An EgoStat.* function takes an egor object and returns a matrix of
# h(e[i]) values, with egos in rows and elements of h(e[i]) in
# columns.

#' @title Helper functions of EgoStats. Some may be documented and exported in the future.
#' 
#' @name EgoStat-internal
#' @keywords internal
#'
#' @param egor an [egor()] object.
#'
NULL

#' @describeIn EgoStat-internal
#'
#' Construct and throw an error message about ego or alter attributes being missing.
.attrErr <- function(term, attrname, req = c("one","both")){
  req <- match.arg(req)
  if(req=="one") stop("In EgoStat ",sQuote(term)," attribute ", sQuote(attrname), " is observed for neither egos nor alters.", call.=FALSE)
  else stop("In EgoStat ",sQuote(term)," attribute ", sQuote(attrname), " must be observed for both egos and alters.", call.=FALSE)
}

#' @describeIn EgoStat-internal
#'
#' Convert a factor to its "ordinary" vector representation.
#' @param x,bin a vector.
.unfactor <- function(x){
  if(is.factor(x)){
    levels(x)[as.integer(x)]
  }else x
}

#' @describeIn EgoStat-internal
#'
#' Degree sequence of the egos.
.degreeseq <- function(egor){
  tabulate(.alterEgos(egor), nbins=nrow(egor$ego))
}

#' @describeIn EgoStat-internal
#'
#' Ego index (i.e., row) associated with each alter.
.alterEgos <- function(egor){
  match(egor$alter$.egoID, as_tibble(egor$ego)$.egoID)
}

#' @describeIn EgoStat-internal
#'
#' As [sum()], but "extrapolating" the `NA`s.
.exsum <- function(x){
  if(length(x)) sum(x, na.rm=TRUE)/mean(!is.na(x)) else 0
}

#' @describeIn EgoStat-internal
#'
#' As [tabulate()], but "extrapolating" the `NA`s.
.extabulate <- function(bin, nbins = max(1, bin, na.rm = TRUE)){
  if(length(bin)) tabulate(bin, nbins)/mean(!is.na(bin)) else numeric(nbins)
}

#' @describeIn EgoStat-internal
#'
#' As [base::match()], but `NA`s are passed through.
.matchNA <- function(x, table, nomatch = NA_integer_, incomparables = NULL){
  as.integer(ifelse(is.na(x), NA, match(x, table, nomatch, incomparables)))
}

#' @describeIn EgoStat-internal
#'
#' As [paste()], but `NA`s are passed through rather than typeset.
.pasteNA <- function(..., sep = " ", collapse = NULL){
  ifelse(apply(do.call(cbind,lapply(list(...), is.na)),1,any), NA, paste(..., sep=sep, collapse=collapse))
}

#' @describeIn EgoStat-internal
#'
#' As [sapply()], but if the function returns a vector, return those
#' vectors in rows, and if it returns a scalar, return a column
#' vector.
.sapply_col <- function(...,simplify=TRUE){
  if(!simplify) stop("This helper function does not make sense with simplify=FALSE.")
  o <- sapply(..., simplify=simplify)
  if(is.null(dim(o))) cbind(o) else t(o)
}

#' @describeIn EgoStat-internal
#'
#' As [mapply()], but if the function returns a vector, return those
#' vectors in rows, and if it returns a scalar, return a column
#' vector.
.mapply_col <- function(..., SIMPLIFY=TRUE){
  if(!SIMPLIFY) stop("This helper function does not make sense with SIMPLIFY=FALSE.")
  o <- mapply(..., SIMPLIFY=SIMPLIFY)
  if(is.null(dim(o))) cbind(o) else t(o)
}

#' @describeIn EgoStat-internal
#'
#' Split an arbitrary vector or matrix by egos.
split_egos_by_ego <- function(x, egor){
  f <- factor(seq_len(nrow(egor$ego)))
  if(is.matrix(x) || is.array(x)) split(x, f, margin=1) else split(x, f)
}

#' @describeIn EgoStat-internal
#'
#' Split an arbitrary vector or matrix by egos.
split_alters_by_ego <- function(x, egor){
  f <- factor(egor$alter$.egoID, levels = as_tibble(egor$ego)$.egoID)
  if(is.matrix(x) || is.array(x)) split(x, f, margin=1) else split(x, f)
}

#' @describeIn EgoStat-internal
#'
#' Split an arbitrary vector or matrix by egos.
split_aaties_by_ego <- function(x, egor){
  f <- factor(egor$aatie$.egoID, levels = as_tibble(egor$ego)$.egoID)
  if(is.matrix(x) || is.array(x)) split(x, f, margin=1) else split(x, f)
}

#' @describeIn EgoStat-internal
#'
#' Pass through the input if it does not contain mising values; otherwise stop with given error. Arguments after the first are passed through to [stop()].
#' 
#' @param ... Additional arguments to subroutines.
#' 
.checkNA <- function(x, ...) if(any(is.na(x))) stop(...) else x

#' Generate the ego contribution matrix from an [egor()] object.
#'
#' @param egor an [egor()] object containing the egocentric dataset.
#' @param h a function taking a row of a [tibble()] and outputting a
#'   numeric vector of dyadwise contributions.
#' @param cn a vector of column names for the output.
#'
#' @return A numeric matrix with `nrow(egor$ego)` rows.
#' @noRd
.eval.h <- function(egor, h, cn, order=1){
  h <- apply(strip_ego_design(egor), 1, h)
  if(is.matrix(h)) h <- t(h) # apply() builds a matrix with egos in columns
  else h <- cbind(h)
  colnames(h) <- cn
  attr(h, "order") <- order
  h
}

#' \code{\link[ergm]{ergm}} Terms Implemented for
#' \code{\link{egor}}
#' 
#' This page describes the \code{\link[ergm]{ergm}} terms (and hence network
#' statistics) for which inference based on egocentrically sampled data is
#' implemented in \code{ergm.ego} package. Other packages may add their own
#' terms. These functions should not be called by the end-user.
#' 
#' The current recommendation for any package implementing additional
#' egocentric calculator terms is to create a help file with a name or alias
#' \code{ergm.ego-terms}, so that \code{help("ergm.ego-terms")} will
#' list egocentric ERGM terms available from all loaded packages.
#' 
#' 
#' @name ergm.ego-terms
#' @aliases ergm.ego-terms terms-ergm.ego ergm.ego.terms terms.ergm.ego
#' ergm-terms ergm.terms terms-ergm terms.ergm EgoStat
#' @docType methods
#' @section Currently implemented egocentric statistics: For each of these,
#' please see their respective package's \code{ergm-terms} help for meaning and
#' parameters. The simplest way to do this is usually via \code{? TERM}.
#' 
#' \describe{ \item{Special-purpose terms:}{ \describe{
#' \item{netsize.adj(edges=+1, mutual=0, transitiveties=0,
#' cyclicalties=0)}{A special-purpose term equivalent to a linear
#' combination of \code{\link[ergm]{edges}},
#' \code{\link[ergm]{mutual}}, \code{\link[ergm]{transitiveties}}, and
#' \code{\link[ergm]{cyclicalties}}, to house the network-size
#' adjustment offset. This term is added to the model automatically
#' and should not be used in the model formula directly.  } } }
#' 
#' \item{ergm:}{
#' * `offset`
#' * `edges`
#' * `nodecov`
#' * `nodefactor`
#' * `nodematch`
#' * `nodemix`
#' * `absdiff`
#' * `degree`
#' * `degrange`
#' * `concurrent`
#' * `concurrentties`
#' * `degree1.5`
#' * `transitiveties`
#' * `cyclicalties`
#' * `esp`
#' * `gwesp`
#' * `gwdegree`
#' * `mm`
#' * `meandeg`*
#' }
#' 
#' \item{tergm:}{
#' * `mean.age`*
#' }
#' }
#'
#' Starred terms are *nonscaling*, in that while they can be
#' evaluated, some inferential results and standard error calculation
#' methods may not be applicable.
#'
#' @seealso \code{\link[ergm]{ergm-terms}}
#' @keywords models
NULL

# copied from ergm
nodecov_names <- function(nodecov, prefix=NULL){
  cn <- if(is.matrix(nodecov)){
          cn <- colnames(nodecov)
          if(is.null(cn) || all(cn==seq_along(cn))) paste(attr(nodecov, "name"), seq_along(cn), sep=".")
          else cn
        }else attr(nodecov, "name")
  NVL3(prefix, paste0(prefix,".",cn), cn)
}
LEVELS_BASE1 <- NULL

EgoStat.offset <- function(egor, trm){
  trm <- substitute(trm)
  if(is.call(trm)){
    init.call <- list(as.name(paste("EgoStat.", trm[[1]],sep="")),egor=egor)
    init.call <- c(init.call,as.list(trm[-1]))
  }else{
    init.call <- list(as.name(paste("EgoStat.", trm,sep="")),egor=egor)
  }
  h <- eval(as.call(init.call))
  
  # needed to match ergm's naming convention for offsets  
  colnames(h) <- paste0("offset(", colnames(h), ")")
  h
}

EgoStat.edges <- function(egor){
  structure(matrix(.degreeseq(egor)/2, dimnames=list(NULL, "edges")), order=1)
}

EgoStat.nodecov <- function(egor, attr){
  egos <- as_tibble(egor$ego)
  alters <- egor$alter

  xe <- ergm.ego_get_vattr(attr, egos, accept = "numeric", multiple = "matrix")
  xa <- ERRVL(try(ergm.ego_get_vattr(attr, alters, accept = "numeric", multiple = "matrix"), silent=TRUE), NULL)

  attrnames <- if(is.matrix(xe)) colnames(xe) else attributes(xe)$name
  alt <- !is.null(xa)

  if(alt){
    xal <- split_alters_by_ego(as.matrix(xa), egor) %>% map(as.matrix)
    xe <- as.matrix(xe)*sapply(xal,nrow)
    xal <- .sapply_col(xal, apply, 2, .exsum)
  }else{
    xe <- as.matrix(xe)*.degreeseq(egor)
    xal <- 0
  }

  structure((xe + xal)/if(alt) 2 else 1, dimnames = list(NULL, paste("nodecov",attrnames,sep=".")), order=1)
}

EgoStat.nodefactor <- function(egor, attr, base=1, levels=LEVELS_BASE1){
  if(!missing(base)) message("In term `nodefactor' in package `ergm.ego': Argument \"base\" has been superseded by \"levels\" and it is recommended to use the latter.  Note that its interpretation may be different.")

  egos <- as_tibble(egor$ego)
  alters <- egor$alter

  xe <- ergm.ego_get_vattr(attr, egos)
  xa <- ERRVL(try(ergm.ego_get_vattr(attr, alters), silent=TRUE), NULL)

  attrname <- attributes(xe)$name
  alt <- !is.null(xa)
  
  levs <- ergm.ego_attr_levels(levels, c(xe, xa), egor, sort(unique(c(xe, xa))))
  if(!is.null(base) && !identical(base,0) && missing(levels)) levs <- levs[-base]

  if(alt){
    xe <- match(xe, levs, 0)
    xa <- .matchNA(xa, levs, 0)
    xal <- split_alters_by_ego(xa, egor)

    xe <- .sapply_col(xe, .extabulate, length(levs))*sapply(xal,length)
    xal <- .sapply_col(xal, .extabulate, length(levs))
  }else{
    xe <- .sapply_col(match(xe, levs, 0), .extabulate, length(levs))*.degreeseq(egor)
    xal <- 0
  }

  structure((xe + xal)/if(alt) 2 else 1, dimnames = list(NULL, paste("nodefactor",attrname,levs,sep=".")), order=1)
}

EgoStat.nodematch <- function(egor, attr, diff=FALSE, keep=NULL, levels=NULL){
  if(!missing(keep)) message("In term `nodematch' in package `ergm.ego': Argument \"keep\" has been superseded by \"levels\" and it is recommended to use the latter.  Note that its interpretation may be different.")
  
  egos <- as_tibble(egor$ego)
  alters <- egor$alter

  xe <- ergm.ego_get_vattr(attr, egos)
  xa <- ergm.ego_get_vattr(attr, alters)

  attrname <- attributes(xe)$name
  
  levs <- ergm.ego_attr_levels(levels, c(xe, xa), egor, sort(unique(c(xe, xa))))
  if(!is.null(keep) && missing(levels)) levs <- levs[keep]

  xe <- match(xe, levs, 0)
  xa <- .matchNA(xa, levs, 0)

  xal <- split_alters_by_ego(xa, egor)

  nlevs <- length(levs)
  combine <- if(diff) identity else sum
  h <- .mapply_col(function(e,a) combine(tabulate(e, nbins=nlevs)*.exsum(e==a)/2), xe, xal, SIMPLIFY=TRUE)
  colnames(h) <- if(diff) paste("nodematch",attrname,levs,sep=".") else paste("nodematch",attrname,sep=".")
  attr(h, "order") <- 1
  h
}

EgoStat.nodemix <- function(egor, attr, base=NULL, levels=NULL, levels2=-1){
  if(!missing(base)) message("In term `nodemix' in package `ergm.ego': Argument \"base\" has been superseded by \"levels2\" and it is recommended to use the latter.  Note that its interpretation may be different.")
  
  egos <- as_tibble(egor$ego)
  alters <- egor$alter
  
  xeval <- ergm.ego_get_vattr(attr, egos)
  xaval <- ergm.ego_get_vattr(attr, alters)

  attrname <- attributes(xeval)$name
  
  levs <- ergm.ego_attr_levels(levels, c(xeval, xaval), egor, sort(unique(c(xeval, xaval))))

  xe <- match(xeval, levs, 0)
  xa <- .matchNA(xaval, levs, 0)
 
  nr <- length(levs)
  nc <- length(levs)

  levels2.list <- transpose(expand.grid(row = levs, col = levs, stringsAsFactors=FALSE))
  indices2.grid <- expand.grid(row = 1:nr, col = 1:nc)
  uun <- as.vector(outer(levs,levs,paste,sep="."))
    
  rowleqcol <- indices2.grid$row <= indices2.grid$col
  levels2.list <- levels2.list[rowleqcol]
  indices2.grid <- indices2.grid[rowleqcol,]
  uun <- uun[rowleqcol]

  levels2.sel <- if(!is.null(base) && !identical(base,0) && missing(levels2)) levels2.list[-base]
                 else ergm.ego_attr_levels(levels2, list(row = c(xeval, xaval), col = c(xaval, xeval)), egor, levels2.list)

  rows2keep <- match(levels2.sel,levels2.list, NA)
  rows2keep <- rows2keep[!is.na(rows2keep)]
  
  u <- indices2.grid[rows2keep,,drop=FALSE]
  namevec <- uun[rows2keep]

  xal <- split_alters_by_ego(xa, egor)

  xeal <- mapply(function(e,a) unlist(apply(cbind(row=pmin(e,a),col=pmax(e,a)),1,list),recursive=FALSE), xe, xal, SIMPLIFY=FALSE)
  # xeal is now a list (over egos) of lists (over alters) of sorted pairs of edge category codes.
  # Missing alter attribute causes both to be NA.
  xeal <- lapply(xeal, lapply, function(xea) if(identical(xea,c(row=NA_integer_,col=NA_integer_))) NA else xea)
  # And now they've been replaced by NA.

  u <- lapply(split(u, seq_len(nrow(u))), unlist)
  h <- .sapply_col(xeal, function(xea) .extabulate(.matchNA(xea,u,0), nbins=length(u))/2)
  colnames(h) <- paste("mix",attrname,namevec,sep=".")
  attr(h, "order") <- 1
  h
}

EgoStat.absdiff <- function(egor, attr, pow=1){
  egos <- as_tibble(egor$ego)
  alters <- egor$alter
  
  xe <- ergm.ego_get_vattr(attr, egos, accept = "numeric")
  xa <- ergm.ego_get_vattr(attr, alters, accept = "numeric")
  
  attrname <- attributes(xe)$name

  xal <- split_alters_by_ego(xa, egor)

  h <- .mapply_col(function(e,a) .exsum(abs(e-a)^pow)/2, xe, xal, SIMPLIFY=TRUE)
  colnames(h) <- if(pow==1) paste("absdiff",attrname,sep=".") else paste("absdiff",pow,".",attrname,sep="")
  attr(h, "order") <- 1
  h
}

EgoStat.degree <- function(egor, d, by=NULL, homophily=FALSE, levels=NULL){
  ## if(any(d==0)) warning("degree(0) (isolate) count statistic depends strongly on the specified population network size.")

  egos <- as_tibble(egor$ego)
  alters <- egor$alter

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    xa <- ERRVL(try(ergm.ego_get_vattr(by, alters), silent=TRUE), NULL)
    
    by <- attributes(xe)$name
  
    levs <- ergm.ego_attr_levels(levels, c(xe, xa), egor, sort(unique(c(xe, xa))))

    xe <- match(xe, levs, 0)
    xa <- NVL3(xa, match(., levs, 0))
  }

  alt <- !is.null(by) && !is.null(alters[[by]])
  if(homophily && !alt) stop("Attribute ", sQuote(by), " must be observed on alters if homophily=TRUE.")

  if(!is.null(by) && homophily){
    xal <- split_alters_by_ego(xa, egor)
    xal <- mapply(function(e,a) a[e==a], xe, xal, SIMPLIFY=FALSE)
    alterct <- sapply(xal, length)
  }else{
    alterct <- .degreeseq(egor)
  }

  if(!is.null(by) && !homophily){
    bys <- rep(seq_along(levs),each=length(d))
    bynames <- rep(levs,each=length(d))
    degs <- rep(d,length(levs))
    cn <- paste0("deg",degs,".",by,bynames)
    h <- .mapply_col(function(e,c) as.numeric(c==degs & e==bys), xe, alterct)
  }else if(homophily){
    cn <-  paste0("deg",d,".homophily.",by)
    h <- .sapply_col(alterct, function(c) as.numeric(c==d))
  }else{
    cn <-  paste0("degree",d)
    h <- .sapply_col(alterct, function(c) as.numeric(c==d))
  }

  colnames(h) <- cn
  attr(h, "order") <- 1
  h
}

EgoStat.degrange <- function(egor, from=NULL, to=Inf, by=NULL, homophily=FALSE, levels=NULL){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  
  egos <- as_tibble(egor$ego)
  alters <- egor$alter
  
  to <- ifelse(to==Inf, .Machine$integer.max, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    xa <- ERRVL(try(ergm.ego_get_vattr(by, alters), silent=TRUE), NULL)
    
    by <- attributes(xe)$name
  
    levs <- ergm.ego_attr_levels(levels, c(xe, xa), egor, sort(unique(c(xe, xa))))

    xe <- match(xe, levs, 0)
    xa <- NVL3(xa, match(., levs, 0))
  }

  alt <- !is.null(by) && !is.null(xa)
  if(homophily && !alt) stop("Attribute ", sQuote(by), " must be observed on alters if homophily=TRUE.")

  if(!is.null(by) && homophily){
    xal <- split_alters_by_ego(xa, egor)
    xal <- mapply(function(e,a) a[e==a], xe, xal, SIMPLIFY=FALSE)
    alterct <- sapply(xal, length)
  }else{
    alterct <- .degreeseq(egor)
  }

  if(!is.null(by) && !homophily){
    bys <- rep(seq_along(levs),each=length(from))
    bynames <- rep(levs,each=length(from))
    froms <- rep(from,length(levs))
    tos <- rep(to,length(levs))

    cn <- ifelse(tos>=.Machine$integer.max,
                 paste0("deg", from, "+.",          by, bynames),
                 paste0("deg", from, "to", to, ".", by, bynames))
    h <- .mapply_col(function(e,c) as.numeric(c>=froms & c<tos & e==bys), xe, alterct)
  }else if(homophily){
    cn <- ifelse(to>=.Machine$integer.max,
                 paste0("deg", from,  "+",     ".homophily.", by),
                 paste0("deg", from, "to", to, ".homophily.", by))
    h <- .sapply_col(alterct, function(c) as.numeric(c>=from & c<to))
  }else{
    cn <- ifelse(to>=.Machine$integer.max,
                 paste0("deg", from,  "+"),
                 paste0("deg", from, "to", to))
    h <- .sapply_col(alterct, function(c) as.numeric(c>=from & c<to))
  }

  colnames(h) <- cn
  attr(h, "order") <- 1
  h
}

EgoStat.concurrent <- function(egor, by=NULL, levels=NULL){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  
  egos <- as_tibble(egor$ego)
  alters <- egor$alter

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    by <- attributes(xe)$name  
    levs <- ergm.ego_attr_levels(levels, xe, egor, sort(unique(c(xe))))
    xe <- match(xe, levs, 0)
  }

  alterct <- .degreeseq(egor)

  if(!is.null(by)){
    cn <- paste0("concurrent.",by,levs)
    h <- .mapply_col(function(e,c) as.numeric(c>=2 & e==seq_along(levs)), xe, alterct)
  }else{
    cn <-  "concurrent"
    h <- .sapply_col(alterct, function(c) as.numeric(c>=2))
  }

  colnames(h) <- cn
  attr(h, "order") <- 1
  h
}

EgoStat.concurrentties <- function(egor, by=NULL, levels=NULL){
  
  egos <- as_tibble(egor$ego)
  alters <- egor$alter

  if(!is.null(by)) {
    xe <- ergm.ego_get_vattr(by, egos)
    by <- attributes(xe)$name  
    levs <- ergm.ego_attr_levels(levels, xe, egor, sort(unique(c(xe))))
    xe <- match(xe, levs, 0)
  }

  alterct <- .degreeseq(egor)

  if(!is.null(by)){
    cn <- paste0("concurrentties.",by,levs)
    h <- .mapply_col(function(e,c) max(c-1,0)*as.numeric(e==seq_along(levs)), xe, alterct)
  }else{
    cn <-  "concurrentties"
    h <- .sapply_col(alterct, function(c) max(c-1,0))
  }

  colnames(h) <- cn
  attr(h, "order") <- 1
  h
}

EgoStat.degree1.5 <- function(egor){
  egor <- as_nested_egor(egor)
  h <- function(e) nrow(e$.alts)^(3/2)
  .eval.h(egor, h, "degree1.5", order=1)
}

EgoStat.transitiveties <- function(egor, attr=NULL, diff=FALSE, levels=TRUE){
  egos <- as_tibble(egor$ego)
  alters <- egor$alter
  aaties <- egor$aatie

  combine <- if(diff) identity else sum

  aal <- split_aaties_by_ego(as.matrix(aaties[,c(".srcID",".tgtID"),drop=FALSE]), egor)

  if(!is.null(attr)){
    xe <- ergm.ego_get_vattr(attr, egos)
    xa <- ergm.ego_get_vattr(attr, alters)
    attrname <- attributes(xe)$name  
    levs <- ergm.ego_attr_levels(levels, xe, egor, sort(unique(c(xe, xa))))
    xe <- match(xe, levs, 0)
    xa <- match(xa, levs, 0)
    xa <- tibble(.altID=alters$.altID, x=xa)
    xal <- NVL3(xa, split_alters_by_ego(., egor))

    if(is.null(xe) || is.null(xa)) .attrErr("transitiveties and cyclicalties", attrname, "both")
    aal <- mapply(function(e, a, aa) aa[e==a$x[match(aa[,1],a$.altID)] & e==a$x[match(aa[,2],a$.altID)], , drop=FALSE], xe, xal, aal, SIMPLIFY=FALSE)

    nlevs <- length(levs)
    # Implement Krivitsky and Morris (2017, p. 490) This works
    # because we want to count how many alters have at least one
    # alter-alter tie, thus forming a transitive tie.
    h <- .mapply_col(function(e, aa)
      combine(length(unique(union(aa[,1], aa[,2])))*tabulate(e,nlevs))/2, xe, aal)
    colnames(h) <- paste("transitiveties",attrname,sep=".")
  }else{
    h <- .sapply_col(aal, function(aa) length(unique(union(aa[,1], aa[,2])))/2)
    colnames(h) <- "transitiveties"
  }
  attr(h, "order") <- 3
  h
}

EgoStat.cyclicalties <- EgoStat.transitiveties

EgoStat.esp <- function(egor, d){
  egor <- as_nested_egor(egor)
  h <- function(e){
    aaties <- unique(
      cbind(pmin(e$.aaties$.srcID,e$.aaties$.tgtID),
            pmax(e$.aaties$.srcID,e$.aaties$.tgtID))
    )
    sp <- sapply(e$.alts$.altID, # for each alter
           function(a) length(unique(c(aaties[aaties[,1]==a,2],aaties[aaties[,2]==a,1]))) # Number of shared partners
           )
    sapply(d, function(k) sum(sp==k))/2
  }
    
  .eval.h(egor, h,
          paste0("esp",d),
          3
          )
}

EgoStat.gwesp <- function(egor, decay=NULL, fixed=FALSE, cutoff=30, alpha=NULL){
  maxesp <- cutoff # Hopefully, network.size > cutoff

  esp <- EgoStat.esp(egor, 1:maxesp)

  if(fixed==FALSE){
    colnames(esp) <- paste0("esp#", 1:maxesp)
    return(esp)
  }

  eta <- exp(decay)*(1-(1-exp(-decay))^(1:maxesp))
  hv <- esp%*%eta
  colnames(hv) <- paste0("gwesp.fixed.",decay)
  attr(hv, "order") <- 3
  hv
}

EgoStat.gwdegree <- function(egor, decay=NULL, fixed=FALSE, cutoff=30){
  maxdeg <- cutoff # Hopefully, network.size > cutoff

  deg <- EgoStat.degree(egor, 1:maxdeg)

  if(fixed==FALSE){
    colnames(deg) <- paste0("gwdegree#", 1:maxdeg)
    return(deg)
  }
  
  eta <- exp(decay)*(1-(1-exp(-decay))^(1:maxdeg))
  hv <- deg%*%eta
  colnames(hv) <- paste0("gwdeg.fixed.",decay)
  attr(hv, "order") <- 1
  hv
}

EgoStat.mm <- function(egor, attrs, levels=NULL, levels2=-1){
  egor <- as_nested_egor(egor)

  aei <- rep(seq_len(nrow(egor)), map_int(egor$.alts, nrow))
  
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
                   sum(sapply(egor$.alts, nrow))*2
                   ),
             name = ".",
             levels = 0,
             levelcodes = 0,
             id = rep(merge(data.frame(i=seq_len(nrow(egor))),
                            data.frame(i=aei))$i,
                      2)
             )
      }else{
        xe <- ERRVL(ec <- try(ergm.ego_get_vattr(spec$attrs, egor), silent=TRUE), NULL)
        xa <- ERRVL(try(ergm.ego_get_vattr(spec$attrs, bind_rows(egor$.alts)), silent=TRUE), NULL)
        name <- attr(NVL(xe,xa), "name")
        if(is.null(xe)&&is.null(xa)) stop(attr(ec, "condition"), call.=FALSE) # I.e., they were both errors. => propagate error message.
        name <- NVL(attr(xe, "name"),attr(xa, "name"))
        xe <- NVL2(xe,
                   data.frame(i=seq_len(nrow(egor)), xe=xe, stringsAsFactors=FALSE),
                   data.frame(i=seq_len(nrow(egor)), stringsAsFactors=FALSE))
        xa <- NVL2(xa,
                   data.frame(i=aei, xa=xa, stringsAsFactors=FALSE),
                   data.frame(i=aei, stringsAsFactors=FALSE))
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
      v$levels <- ergm.ego_attr_levels(v$levels, v$val, egor, levels=v$unique)
      v$levelcodes <- seq_along(v$levels)
      v$valcodes <- .matchNA(v$val, v$levels, nomatch=0)
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
  levels2sel <- ergm.ego_attr_levels(levels2, list(row=attrval$row$val, col=attrval$col$val), egor, levels=levels2vals)
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

  h1 <- attrval %>%
    map("valcodes") %>%
    transpose() %>%
    map(~if(is.na(.$row) || is.na(.$col)) NA else .) %>%
    .matchNA(levels2codes, 0L) %>%
    tapply(., rep(aei, length.out=length(.)), identity, simplify=FALSE) %>%
    map(.extabulate, length(levels2codes)) %>%
    do.call(rbind,.)

  h <- matrix(0, nrow(egor), ncol=ncol(h1))
  h[sort(unique(aei)),] <- h1/2
  colnames(h) <- coef.names

  if(symm){
    selff <- 1+map_lgl(levels2codes, all_identical)
    h <- sweep(h, 2, selff, `/`)
  }
 
  attr(h, "order") <- 1
  h
}

EgoStat.meandeg <- function(egor){
  out <- summary(egor~edges, individual=FALSE, scaleto=2)
  names(out) <- "meandeg"
  attr(out, "order") <- 0
  out
}
