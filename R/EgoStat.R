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
  do.call(rbind, egor$.alts)
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
#' \item{ergm:}{ \itemize{ \item \code{edges} \item \code{nodecov}
#' \item \code{nodefactor} \item \code{nodematch} \item \code{nodemix} \item
#' \code{absdiff} \item \code{degree} \item \code{degrange} \item
#' \code{concurrent} \item \code{concurrentties} \item \code{degreepopularity}
#' \item \code{transitiveties} \item \code{cyclicalties} \item \code{esp} \item \code{gwesp}
#' } }
#' 
#' \item{tergm:}{ \itemize{ \item \code{mean.age} } } }
#'
#' @param egor,attrname,base,diff,keep,pow,d,by,homophily,from,to,decay,fixed,cutoff,alpha,emptyval,nw,arglist,... arguments to terms. See \code{\link[ergm]{ergm-terms}}.
#' @seealso \code{\link[ergm]{ergm-terms}}
#' @keywords models
NULL

#' @export
#' @rdname ergm.ego-terms
EgoStat.edges <- function(egor){
  h <- function(e) nrow(e$.alts)/2
  .eval.h(egor, h, "edges")
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.nodecov <- function(egor, attrname){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))

  if(nattr==0) .attrErr("nodecov", attrname, "one")
  
  h <- function(e) (sum(e[[attrname]])*nrow(e$.alts) + sum(e$.alts[[attrname]]))/nattr
  .eval.h(egor, h, paste("nodecov",attrname,sep="."))
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.nodefactor <- function(egor, attrname, base=1){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
  if(nattr==0) .attrErr("nodefactor", attrname, "one")

  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "base" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  if(length(base)!=0 && !identical(as.integer(base),as.integer(0))) l <- l[-base]
  nl <- length(l)

  h <- function(e)
  (tabulate(match(e[[attrname]],l,0), nbins=nl)*nrow(e$.alts)
    + tabulate(match(e$.alts[[attrname]],l,0), nbins=nl))/nattr
  
  .eval.h(egor, h, paste("nodefactor",attrname,l,sep="."))
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.nodematch <- function(egor, attrname, diff=FALSE, keep=NULL){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
  if(nattr==0) .attrErr("nodematch", attrname, "both")
  
  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "non-keep" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  l <- l[NVL(keep,TRUE)]
  nl <- length(l)

  combine <- if(diff) identity else sum
  h <- function(e)
    combine(tabulate(match(e[[attrname]],l,0), nbins=nl)*sum(e[[attrname]]==e$.alts[[attrname]])/2)

  .eval.h(egor, h,
          if(diff) paste("nodematch",attrname,l,sep=".")
          else paste("nodematch",attrname,sep="."))
}


#' @export
#' @rdname ergm.ego-terms
EgoStat.nodemix <- function(egor, attrname, base=NULL){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
  if(nattr==0) .attrErr("nodemix", attrname, "both")

  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "base" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  l <- outer(l,l,paste,sep=".")
  l <- l[upper.tri(l,diag=TRUE)]
  if (length(base) && !identical(as.integer(base),as.integer(0))) l <- l[-base]
  nl <- length(l)
  h <- function(e)
    tabulate(match(paste(pmin(e[[attrname]],e$.alts[[attrname]]),
                         pmax(e[[attrname]],e$.alts[[attrname]]),
                         sep="."),l,0), nbins=nl)/2
  .eval.h(egor, h,
          paste("mix",attrname,l,sep="."))
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.absdiff <- function(egor, attrname, pow=1){
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
  if(nattr==0) .attrErr("absdiff", attrname, "both")
  
  h <- function(e)
    sum(abs(e[[attrname]]-e$.alts[[attrname]])^pow)/2

  .eval.h(egor, h,
          if(pow==1) paste("absdiff",attrname,sep=".")
          else paste("absdiff",pow,".",attrname,sep=""))
}

#' @export
#' @rdname ergm.ego-terms
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
    h <- function(e) as.numeric(nrow(e$.alts)==degs & e[[by]]==bys)
  }else if(homophily){
    cn <-  paste0("deg",d,".homophily.",by)
    h <- function(e) as.numeric(sum(e[[by]]==e$.alts[[by]])==d)
  }else{
    cn <-  paste0("degree",d)
    h <- function(e) as.numeric(nrow(e$.alts)==d)
  }

  .eval.h(egor, h, cn)
}

#' @export
#' @rdname ergm.ego-terms
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
    h <- function(e) as.numeric(nrow(e$.alts)>=froms & nrow(e$.alts)<tos & e[[by]]==bys)
  }else if(homophily){
    cn <- ifelse(to>=.Machine$integer.max,
                 paste0("deg", from,  "+",     ".homophily.", by),
                 paste0("deg", from, "to", to, ".homophily.", by))
    h <- function(e) as.numeric(sum(e[[by]]==e$.alts[[by]])>=from & sum(e[[by]]==e$.alts[[by]])<to)
  }else{
    cn <- ifelse(to>=.Machine$integer.max,
                 paste0("deg", from,  "+"),
                 paste0("deg", from, "to", to))
    h <- function(e) as.numeric(nrow(e$.alts)>=from & nrow(e$.alts)<to)
  }

  .eval.h(egor, h, cn)
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.concurrent <- function(egor, by=NULL){

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("concurrent")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  if(!is.null(by)){
    l <- sort(unique(c(egor[[by]],.allAlters(egor)[[by]])))
    nl <- length(l)
  }  

  if(!is.null(by)){
    bys <- l
    cn <- paste0("concurrent.", by, bys)
    h <- function(e) as.numeric(nrow(e$.alts)>=2 & e[[by]]==bys)
  }else{
    cn <-  "concurrent"
    h <- function(e) nrow(e$.alts)>=2
  }

  .eval.h(egor, h, cn)
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.concurrentties <- function(egor, by=NULL){
  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("concurrent")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  if(!is.null(by)){
    l <- sort(unique(c(egor[[by]],.allAlters(egor)[[by]])))
    nl <- length(l)
  }  

  if(!is.null(by)){
    bys <- l
    cn <- paste0("concurrentties.", by, bys)
    h <- function(e) max(nrow(e$.alts)-1,0)*as.numeric(e[[by]]==bys)
  }else{
    cn <-  "concurrentties"
    h <- function(e) max(nrow(e$.alts)-1,0)
  }

  .eval.h(egor, h, cn)
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.degreepopularity <- function(egor){

  h <- function(e) nrow(e$.alts)^(3/2)

  .eval.h(egor, h, "degreepopularity")
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.transitiveties <- function(egor, attrname=NULL){
  if(!is.null(attrname)){
    nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
    if(nattr!=2) .attrErr("transitiveties and cyclicalties", attrname, "both")
    # Note: alterID API is subject to change.
    egor$.matchAttr <- egor[[attrname]]
    egor <- subset(egor,
                   .matchAttr==.alts[[attrname]][.aaties$.srcRow] &
                   .matchAttr==.alts[[attrname]][.aaties$.tgtRow],
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
#' @rdname ergm.ego-terms
EgoStat.cyclicalties <- EgoStat.transitiveties

#' @export
#' @rdname ergm.ego-terms
EgoStat.esp <- function(egor, d){
  h <- function(e){
    aaties <- unique(
      cbind(pmin(e$.aaties$Source,e$.aaties$Target),
            pmax(e$.aaties$Source,e$.aaties$Target))
    )
    sp <- sapply(e$.alts$alterID, # for each alter
           function(a) length(unique(c(aaties[aaties[,1]==a,2],aaties[aaties[,2]==a,1]))) # Number of shared partners
           )
    sapply(d, function(k) sum(sp==k))/2
  }
  
  .eval.h(egor, h,
          paste0("esp",d)
          )
}

#' @export
#' @rdname ergm.ego-terms
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
  hv
}

