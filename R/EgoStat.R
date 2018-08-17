#  File R/EgoStat.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################
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
#' Return a [`tibble`] containing all alters.
.allAlters <- function(egor){
  do.call(rbind, egor$.alts)
}

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
#' As [sum()], but "extrapolating" the `NA`s.
.exsum <- function(x){
  if(length(x)) sum(x, na.rm=TRUE)/mean(!is.na(x)) else 0
}

#' @describeIn EgoStat-internal
#'
#' As [tabulate()], but "extrapolating" the `NA`s.
.extabulate <- function(bin, nbins = max(1, bin, na.rm = TRUE)){
  if(length(bin)) tabulate(bin, nbins)/mean(!is.na(bin)) else rep(0, nbins)
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
#' Pass through the input if it does not contain mising values; otherwise stop with given error. Arguments after the first are passed through to [stop()].
#' 
#' @param ... Additional arguments to subroutines.
#' 
.checkNA <- function(x, ...) if(any(is.na(x))) stop(...) else x

.preproc_factor <- function(egor, attrname){
  if(!is.null(attrname)){
    if(length(attrname)>1){
      # If there are multiple attributes, concatenate their names with a
      # dot and concatenate their values with a dot.
      attrnamename <- paste(attrname, collapse=".")
      NAlist <- apply(is.na(as.tibble(egor)[,attrname]), 1, any)
      egor[[attrnamename]] <- ifelse(NAlist, NA, .unfactor(do.call(paste,c(as.list(as.tibble(egor)[,attrname]),list(sep=".")))))
      egor$.alts <- lapply(egor$.alts, function(a){
        NAlist <- apply(is.na(a[,attrname]), 1, any)
        a[[attrnamename]] <- ifelse(NAlist, NA, .unfactor(do.call(paste,c(as.list(a[,attrname]),list(sep=".")))))
        a
      })
      attrname <- attrnamename
    }else{
      egor[[attrname]] <- .unfactor(egor[[attrname]])
      egor$.alts <- lapply(egor$.alts, function(a){
        a[[attrname]] <- .unfactor(a[[attrname]])
        a
      })
    }
  }
  list(egor=egor, attrname=attrname)
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
.eval.h <- function(egor, h, cn, order=1){
  h <- apply(egor, 1, h)
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
#' }
#' 
#' \item{tergm:}{
#' * `mean.age`
#' }
#' }
#'
#' @param egor The [`egor`] object whose statistics are to be evaluated.
#' @param attrname,base,diff,keep,pow,d,by,homophily,from,to,decay,fixed,cutoff,alpha,emptyval,nw,arglist,levels,levels2,attrs,... arguments to terms. See \code{\link[ergm]{ergm-terms}}.
#' @import zeallot
#' @seealso \code{\link[ergm]{ergm-terms}}
#' @keywords models
NULL

#' @export
#' @param trm [`ergm`] terms to be offset.
#' @rdname ergm.ego-terms
EgoStat.offset <- function(egor, trm){
  trm <- substitute(trm)
  if(is.call(trm)){
    init.call <- list(as.name(paste("EgoStat.", trm[[1]],sep="")),egor=egor)
    init.call <- c(init.call,as.list(trm[-1]))
  }else{
    init.call <- list(as.name(paste("EgoStat.", trm,sep="")),egor=egor)
  }
  eval(as.call(init.call))
}

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
  
  h <- function(e) (sum(e[[attrname]])*nrow(e$.alts) + .exsum(e$.alts[[attrname]]))/nattr
  .eval.h(egor, h, paste("nodecov",attrname,sep="."))
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.nodefactor <- function(egor, attrname, base=1, levels=NULL){
  c(egor, attrname) %<-% .preproc_factor(egor, attrname)
  
  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
  if(nattr==0) .attrErr("nodefactor", attrname, "one")
  
  l <- NVL(levels, sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]]))))
  # Note that all "base" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  if(length(base)!=0 && !identical(as.integer(base),as.integer(0))) l <- l[-base]
  nl <- length(l)

  h <- function(e)
  (tabulate(match(e[[attrname]],l,0), nbins=nl)*nrow(e$.alts)
    + .extabulate(.matchNA(e$.alts[[attrname]],l,0), nbins=nl))/nattr
  
  .eval.h(egor, h, paste("nodefactor",attrname,l,sep="."))
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.nodematch <- function(egor, attrname, diff=FALSE, keep=NULL){
  c(egor, attrname) %<-% .preproc_factor(egor, attrname)

  nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
  if(nattr==0) .attrErr("nodematch", attrname, "both")
  
  l <- sort(unique(c(egor[[attrname]],.allAlters(egor)[[attrname]])))
  # Note that all "non-keep" levels will be matched to 0 and therefore
  # excluded from the tabulation below.
  l <- l[NVL(keep,TRUE)]
  nl <- length(l)

  combine <- if(diff) identity else sum
  h <- function(e)
    combine(tabulate(match(e[[attrname]],l,0), nbins=nl)*.exsum(e[[attrname]]==e$.alts[[attrname]])/2)

  .eval.h(egor, h,
          if(diff) paste("nodematch",attrname,l,sep=".")
          else paste("nodematch",attrname,sep="."))
}


#' @export
#' @rdname ergm.ego-terms
EgoStat.nodemix <- function(egor, attrname, base=NULL){
  c(egor, attrname) %<-% .preproc_factor(egor, attrname)

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
    .extabulate(.matchNA(.pasteNA(pmin(e[[attrname]],e$.alts[[attrname]]),
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
    .exsum(abs(e[[attrname]]-e$.alts[[attrname]])^pow)/2

  .eval.h(egor, h,
          if(pow==1) paste("absdiff",attrname,sep=".")
          else paste("absdiff",pow,".",attrname,sep=""))
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.degree <- function(egor, d, by=NULL, homophily=FALSE, levels=NULL){
  ## if(any(d==0)) warning("degree(0) (isolate) count statistic depends strongly on the specified population network size.")
  c(egor, by) %<-% .preproc_factor(egor, by)

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  alt <- !is.null(by) && !is.null(.allAlters(egor)[[by]])
  if(homophily && !alt) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on both egos and alters if homophily=TRUE.", call.=FALSE)
  
  if(!is.null(by)){
    l <- NVL(levels, sort(unique(c(egor[[by]],.allAlters(egor)[[by]]))))
    nl <- length(l)
  }

  if(!is.null(by) && !homophily){
    bys <- rep(l,each=length(d))
    degs <- rep(d,nl)
    cn <- paste0("deg",degs,".",by,bys)
    h <- function(e) as.numeric(nrow(e$.alts)==degs & e[[by]]==bys)
  }else if(homophily){
    cn <-  paste0("deg",d,".homophily.",by)
    h <- function(e) as.numeric(sum(.checkNA(e[[by]]==e$.alts[[by]], "Degree distribution within group by attribute cannot be estimated when alter attributes are missing at this time."))==d)
  }else{
    cn <-  paste0("degree",d)
    h <- function(e) as.numeric(nrow(e$.alts)==d)
  }

  .eval.h(egor, h, cn)
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.degrange <- function(egor, from=NULL, to=Inf, by=NULL, homophily=FALSE, levels=NULL){
  ## if(any(from==0)) warning("degrange(0,...) (isolate) count depends strongly on the specified population network size.")
  c(egor, by) %<-% .preproc_factor(egor, by)
  
  to <- ifelse(to==Inf, .Machine$integer.max, to)

  if(length(to)==1 && length(from)>1) to <- rep(to, length(from))
  else if(length(from)==1 && length(to)>1) from <- rep(from, length(to))
  else if(length(from)!=length(to)) stop("The arguments of term degrange must have arguments either of the same length, or one of them must have length 1.")
  else if(any(from>=to)) stop("Term degrange must have from<to.")

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  alt <- !is.null(by) && !is.null(.allAlters(egor)[[by]])
  if(homophily && !alt) stop("For term ",sQuote("degree")," attribute ", sQuote(by), " must be observed on both egos and alters if homophily=TRUE.", call.=FALSE)
  
  if(!is.null(by)){
    l <- NVL(levels, sort(unique(c(egor[[by]],.allAlters(egor)[[by]]))))
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
    h <- function(e) as.numeric(sum(.checkNA(e[[by]]==e$.alts[[by]], "Degree distribution within group by attribute cannot be estimated when alter attributes are missing at this time."))>=from & sum(.checkNA(e[[by]]==e$.alts[[by]], "Degree distribution within group by attribute cannot be estimated when alter attributes are missing at this time."))<to)
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
EgoStat.concurrent <- function(egor, by=NULL, levels=NULL){
  c(egor, by) %<-% .preproc_factor(egor, by)

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("concurrent")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  if(!is.null(by)){
    l <- NVL(levels, sort(unique(c(egor[[by]],.allAlters(egor)[[by]]))))
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
EgoStat.concurrentties <- function(egor, by=NULL, levels=NULL){
  c(egor, by) %<-% .preproc_factor(egor, by)

  if(!is.null(by) && !by %in% names(egor)) stop("For term ",sQuote("concurrent")," attribute ", sQuote(by), " must be observed on egos.", call.=FALSE)
  
  if(!is.null(by)){
    l <- NVL(levels, sort(unique(c(egor[[by]],.allAlters(egor)[[by]]))))
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
EgoStat.degree1.5 <- function(egor){
  h <- function(e) nrow(e$.alts)^(3/2)
  .eval.h(egor, h, "degree1.5")
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.transitiveties <- function(egor, attrname=NULL){
  c(egor, attrname) %<-% .preproc_factor(egor, attrname)

  if(!is.null(attrname)){
    nattr <- (attrname %in% names(egor)) + (attrname %in% names(egor$.alts[[1]]))
    if(nattr!=2) .attrErr("transitiveties and cyclicalties", attrname, "both")
    egor <- subset(egor,
                   function(r, attrname)
                     .checkNA(r[[attrname]]==r$.alts[[attrname]][r$.aaties$.srcRow] &
                     r[[attrname]]==r$.alts[[attrname]][r$.aaties$.tgtRow], "Transitive ties count with grouping by attribute cannot be estimated when alter attributes are missing at this time."),
                   attrname=attrname,
                   unit="aatie")
  }
  h <- function(e)
    # Implement Krivitsky and Morris (2017, p. 490) This works
    # because we want to count how many alters have at least one
    # alter-alter tie, thus forming a transitive tie.
    length(unique(union(e$.aaties$.srcID, e$.aaties$.tgtID)))/2

  .eval.h(egor, h,
          if(is.null(attrname)) paste("transitiveties",sep=".")
          else paste("transitiveties",attrname,sep="."),
          3)
}

#' @export
#' @rdname ergm.ego-terms
EgoStat.cyclicalties <- EgoStat.transitiveties

#' @export
#' @rdname ergm.ego-terms
EgoStat.esp <- function(egor, d){
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
  attr(hv, "order") <- 3
  hv
}

#' @export
#' @rdname ergm.ego-terms
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

#' @export
#' @rdname ergm.ego-terms
EgoStat.mm <- function(egor, attrs, levels=NULL, levels2=NULL){

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
        xa <- ERRVL(try(ergm.ego_get_vattr(spec$attrs, .allAlters(egor)), silent=TRUE), NULL)
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
