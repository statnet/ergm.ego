#  File R/EgoStat.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################
# An EgoStat.* function takes an egor object and returns a matrix of
# h(e[i]) values, with egos in rows and elements of h(e[i]) in
# columns.


.allAlters <- function(egor){
  do.call(rbind, egor$.alters)
}

.attrErr <- function(term, attrname, req = c("one","both")){
  req <- match.arg(req)
  if(req=="one") stop("In EgoStat ",sQuote(term)," attribute ", sQuote(attrname), " is observed for neither egos nor alters.", call.=FALSE)
  else stop("In EgoStat ",sQuote(term)," attribute ", sQuote(attrname), " must be observed for both egos and alters.", call.=FALSE)
}

#' Generate the ego contribution matrix from an [egor()] object.
#'
#' @param egor an [egor()] object containing the egocentric dataset.
#' @param h a function taking a row of a [tibble()] and outputting a
#'   numeric vector of dyadwise contributions.
#' @param cn a vector of column names for the output.
#'
#' @return A numeric matrix with `nrow(egor)` rows.
#' @noRd
.eval.h <- function(egor, h, cn){
  h <- apply(egor, 1, h)
  if(is.matrix(h)) h <- t(h) # apply() builds a matrix with egos in columns
  else h <- cbind(h)
  colnames(h) <- cn
  h
}

#' \code{\link[ergm]{ergm}} Terms Implemented for
#' \code{\link{egor}}
#' 
#' This page describes the \code{\link[ergm]{ergm}} terms (and hence network
#' statistics) for which inference based on egocentrically sampled data is
#' implemented in \code{ergm.ego} package. Other packages may add their own
#' terms.
#' 
#' The current recommendation for any package implementing additional
#' egocentric calculator terms is to create a help file with a name or alias
#' \code{ergm.ego-terms}, so that \code{help("ergm.ego-terms")} will
#' list egocentric ERGM terms available from all loaded packages.
#' 
#' 
#' @name ergm.ego-terms
#' @aliases ergm.ego-terms terms-ergm.ego ergm.ego.terms terms.ergm.ego
#' ergm-terms ergm.terms terms-ergm terms.ergm EgoStat EgoStat.edges
#' EgoStat.nodecov EgoStat.nodefactor EgoStat.nodematch EgoStat.nodemix
#' EgoStat.absdiff EgoStat.degree EgoStat.degrange EgoStat.concurrent
#' EgoStat.concurrentties EgoStat.degreepopularity EgoStat.cyclicalties EgoStat.transitiveties EgoStat.mean.age netsize.adj
#' InitErgmTerm.netsize.adj
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
EgoStat.edges <- function(egor){
  h <- function(e) nrow(e$.alters)/2
  .eval.h(egor, h, "edges")
}

#' @export
EgoStat.nodecov <- function(egor, attrname){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alters[[1]]))

  if(nattr==0) .attrErr("nodecov", attrname, "one")
  
  h <- function(e) (sum(e[[attrname]])*nrow(e$.alters) + sum(e$.alters[[attrname]]))/nattr
  .eval.h(egor, h, paste("nodecov",attrname,sep="."))
}

#' @export
EgoStat.nodefactor <- function(egor, attrname, base=1){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alters[[1]]))
  if(nattr==0) .attrErr("nodefactor", attrname, "one")

  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "base" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  if(length(base)!=0 && !identical(as.integer(base),as.integer(0))) l <- l[-base]
  nl <- length(l)

  h <- function(e)
  (tabulate(match(e[[attrname]],l,0), nbins=nl)*nrow(e$.alters)
    + tabulate(match(e$.alters[[attrname]],l,0), nbins=nl))/nattr
  
  .eval.h(egor, h, paste("nodefactor",attrname,l,sep="."))
}

#' @export
EgoStat.nodematch <- function(egor, attrname, diff=FALSE, keep=NULL){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alters[[1]]))
  if(nattr==0) .attrErr("nodematch", attrname, "both")
  
  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "non-keep" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  l <- l[NVL(keep,TRUE)]
  nl <- length(l)

  combine <- if(diff) identity else sum
  h <- function(e)
    combine(tabulate(match(e[[attrname]],l,0), nbins=nl)*sum(e[[attrname]]==e$.alters[[attrname]])/2)

  .eval.h(egor, h,
          if(diff) paste("nodematch",attrname,l,sep=".")
          else paste("nodematch",attrname,sep="."))
}


#' @export
EgoStat.nodemix <- function(egor, attrname, base=NULL){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alters[[1]]))
  if(nattr==0) .attrErr("nodemix", attrname, "both")

  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "base" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  l <- outer(l,l,paste,sep=".")
  l <- l[upper.tri(l,diag=TRUE)]
  if (length(base) && !identical(as.integer(base),as.integer(0))) l <- l[-base]
  nl <- length(l)
  h <- function(e)
    tabulate(match(paste(pmin(e[[attrname]],e$.alters[[attrname]]),
                         pmax(e[[attrname]],e$.alters[[attrname]]),
                         sep="."),l,0), nbins=nl)/2
  .eval.h(egor, h,
          paste("mix",attrname,l,sep="."))
}

#' @export
EgoStat.absdiff <- function(egor, attrname, pow=1){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alters[[1]]))
  if(nattr==0) .attrErr("absdiff", attrname, "both")
  
  h <- function(e)
    sum(abs(e[[attrname]]-e$.alters[[attrname]])^pow)/2

  .eval.h(egor, h,
          if(pow==1) paste("absdiff",attrname,sep=".")
          else paste("absdiff",pow,".",attrname,sep=""))
}

#' @export
EgoStat.degree <- function(egor, d, by=NULL, homophily=FALSE){
  ## if(any(d==0)) warning("degree(0) (isolate) count statistic depends strongly on the specified population network size.")

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  alt <- !is.null(by) && !is.null(.allAlters(egor)[[by]])
  if(homophily && !alt) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on both egos and alters if homophily=TRUE.", call.=FALSE)
  
  if(!is.null(by)){
    l <- sort(unique(c(egor[[by]],.allAlters(egor)[[by]])))
    nl <- length(l)
  }

  if(!is.null(by) && !homophily){
    bys <- rep(l,each=length(d))
    degs <- rep(d,nl)
    cn <- paste0("deg",degs,".",by,bys)
    h <- function(e) as.numeric(nrow(e$.alters)==degs & e[[by]]==bys)
  }else if(homophily){
    cn <-  paste0("deg",d,".homophily.",by)
    h <- function(e) as.numeric(sum(e[[by]]==e$.alters[[by]])==d)
  }else{
    cn <-  paste0("degree",d)
    h <- function(e) as.numeric(nrow(e$.alters)==d)
  }

  .eval.h(egor, h, cn)
}

#' @export
EgoStat.degrange <- function(egor, from=NULL, to=Inf, by=NULL, homophily=FALSE){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  
  to <- ifelse(to==Inf, .Machine$integer.max, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  alt <- !is.null(by) && !is.null(.allAlters(egor)[[by]])
  if(homophily && !alt) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on both egos and alters if homophily=TRUE.", call.=FALSE)
  
  if(!is.null(by)){
    l <- sort(unique(c(egor[[by]],.allAlters(egor)[[by]])))
    nl <- length(l)
  }

  if(!is.null(by) && !homophily){
    bys <- rep(l,each=length(from))
    froms <- rep(from,nl)
    tos <- rep(to,nl)

    cn <- ifelse(tos>=.Machine$integer.max,
                 paste0("deg", from, "+.",          by, bys),
                 paste0("deg", from, "to", to, ".", by, bys))
    h <- function(e) as.numeric(nrow(e$.alters)>=froms & nrow(e$.alters)<tos & e[[by]]==bys)
  }else if(homophily){
    cn <- ifelse(to>=.Machine$integer.max,
                 paste0("deg", from,  "+",     ".homophily.", by),
                 paste0("deg", from, "to", to, ".homophily.", by))
    h <- function(e) as.numeric(sum(e[[by]]==e$.alters[[by]])>=from & sum(e[[by]]==e$.alters[[by]])<to)
  }else{
    cn <- ifelse(to>=.Machine$integer.max,
                 paste0("deg", from,  "+"),
                 paste0("deg", from, "to", to))
    h <- function(e) as.numeric(nrow(e$.alters)>=from & nrow(e$.alters)<to)
  }

  .eval.h(egor, h, cn)
}

#' @export
EgoStat.concurrent <- function(egor, by=NULL){

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("concurrent")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  if(!is.null(by)){
    l <- sort(unique(c(egor[[by]],.allAlters(egor)[[by]])))
    nl <- length(l)
  }  

  if(!is.null(by)){
    bys <- l
    cn <- paste0("concurrent.", by, bys)
    h <- function(e) as.numeric(nrow(e$.alters)>=2 & e[[by]]==bys)
  }else{
    cn <-  "concurrent"
    h <- function(e) nrow(e$.alters)>=2
  }

  .eval.h(egor, h, cn)
}

#' @export
EgoStat.concurrentties <- function(egor, by=NULL){
  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("concurrent")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  if(!is.null(by)){
    l <- sort(unique(c(egor[[by]],.allAlters(egor)[[by]])))
    nl <- length(l)
  }  

  if(!is.null(by)){
    bys <- l
    cn <- paste0("concurrentties.", by, bys)
    h <- function(e) max(nrow(e$.alters)-1,0)*as.numeric(e[[by]]==bys)
  }else{
    cn <-  "concurrentties"
    h <- function(e) max(nrow(e$.alters)-1,0)
  }

  .eval.h(egor, h, cn)
}

#' @export
EgoStat.degreepopularity <- function(egor){

  h <- function(e) nrow(e$.alters)^(3/2)

  .eval.h(egor, h, "degreepopularity")
}

#' @export
EgoStat.transitiveties <- function(egor, attrname=NULL){
  if(!is.null(attrname)){
    nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alters[[1]]))
    if(nattr!=2) .attrErr("transitiveties and cyclicalties", attrname, "both")
    # Note: alterID API is subject to change.
    egor$.matchAttr <- egor[[attrname]]
    egor <- subset(egor,
                   .matchAttr==.alters[[attrname]][.aaties$.srcIx] &
                   .matchAttr==.alters[[attrname]][.aaties$.dstIx],
                   aspect="ties")
  }
  h <- function(e)
    # Implement Krivitsky and Morris (2017, p. 490) This works
    # because we want to count how many alters have at least one
    # alter-alter tie, thus forming a transitive tie.
    length(unique(union(e$.aaties$Source, e$.aaties$Target)))/2

  .eval.h(egor, h,
          if(is.null(attrname)) paste("transitiveties",sep=".")
          else paste("transitiveties",attrname,sep="."))
}

#' @export
EgoStat.cyclicalties <- EgoStat.transitiveties
