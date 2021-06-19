#  File R/summary.statistics.egor.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2021 Statnet Commons
################################################################################
#' Calculation of ERGM-style summary statistics for \code{\link{egor}}
#' objects.
#' 
#' Used to calculate the specified network statistics inferred from a
#' \code{\link{egor}} object.
#' 
#' 
#' @aliases summary_formula.egodata summary summary_formula
#' @param object An \code{\link[ergm]{ergm}}-style formula with a
#' \code{\link{egor}} object as the LHS.
#' 
#' For a list of currently implemented egocentric terms for the RHS, see
#' \code{\link{ergm.ego-terms}}.
#' @param \dots Not used at this time.
#' @param basis An optional \code{\link{egor}} object relative to which the
#' statistics should be calculated.
#' @param individual If \code{FALSE} (the default), calculate the estimated
#' per-capita statistics, weighted according to the ego weights, then scale
#' them up to a network of size \code{scaleto}.
#' 
#' If \code{TRUE}, calculate each ego's individual contribution to the
#' specified network statistics.
#' @param scaleto Size of a hypothetical network to which to scale the
#' statistics. Defaults to the number of egos in the dataset.
#' @return If \code{individual==FALSE}, an `ergm.ego_svystat` object, which is a subclass of [svystat][survey::svymean]---effectively a named vector of statistics. If
#' \code{individual==TRUE}, a matrix with a row for each ego, giving that ego's
#' contribution to the network statistic.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link[ergm]{summary_formula}},
#' \code{\link[ergm]{summary_formula.ergm}}
#' @references
#'
#' * Pavel N. Krivitsky and Martina Morris (2017). "Inference for social network models from egocentrically sampled data, with application to understanding persistent racial disparities in HIV prevalence in the US." *Annals of Applied Statistics*, 11(1): 427–455. \doi{10.1214/16-AOAS1010}
#'
#' * Pavel N. Krivitsky, Mark S. Handcock, and Martina Morris (2011). "Adjusting for
#' Network Size and Composition Effects in Exponential-Family Random Graph
#' Models." \emph{Statistical Methodology}, 8(4): 319–339. \doi{10.1016/j.stamet.2011.01.005}
#'
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#' (nw.summ <- summary(faux.mesa.high~edges+degree(0:3)+nodematch("Race")+
#'                     nodematch("Sex")+absdiff("Grade")+nodemix("Grade")))
#' 
#' (ego.summ <- summary(fmh.ego~edges+degree(0:3)+nodematch("Race")+nodematch("Sex")+
#'                      absdiff("Grade")+nodemix("Grade"),
#'                      scaleto=network.size(faux.mesa.high)))
#' 
#' stopifnot(isTRUE(all.equal(as.vector(nw.summ),as.vector(ego.summ))))
#'
#' @importFrom survey svymean
#' @export
summary_formula.egor <- function(object,..., basis=NULL, individual=FALSE, scaleto=NULL){
  egor <-
    if(!is.null(basis)) basis
    else eval_lhs.formula(object)

  stats <- lapply(list_rhs.formula(object), function(trm){
    if(is.call(trm)){
      egostat <- locate_prefixed_function(trm[[1]], "EgoStat", "Egocentric statistic", env=environment(object))
      init.call <- list(egostat, egor=egor)
      init.call <- c(init.call,as.list(trm[-1]))
    }else{
      egostat <- locate_prefixed_function(trm, "EgoStat", "Egocentric statistic", env=environment(object))
      init.call <- list(egostat, egor=egor)
    }
    stat <- eval(as.call(init.call), environment(object))
    if(attr(stat, "order")==0 && individual)
      stop("Nonscaling statistic detected. Individual contributions are meaningless.")
    stat
  })

  orders <- sapply(stats, attr, "order")

  scaling <- rep.int(orders != 0, ifelse(orders == 0, lengths(stats), sapply(stats, ncol)))
  
  stats <-
    if(!individual){
      if(any(scaling)){
        s <- do.call(cbind, stats[orders!=0])
        s <- NVL3(ego_design(egor),
                              svymean(s, ., ...),
                              structure(colMeans(s),
                                        var=cov(s)/nrow(s),
                                        statistic="mean", class="svystat")
                       )
      }else s <- NULL

      if(any(!scaling)) s <- combine_stats(c(list(s), stats[orders==0]))

      # Permutation that maps scaling and nonscaling stats back to their original positions:
      p <- integer(length(scaling))
      p[scaling] <- seq_len(sum(scaling))
      p[!scaling] <- sum(scaling) + seq_len(sum(!scaling))

      s <- structure(s[p], var = attr(s, "var")[p,p],
                     statistic="scaled mean", class="svystat")

      attr(s,"order") <- unique(orders)
      attr(s,"scaling") <- scaling
      class(s) <- c("ergm.ego_svystat", "svystat")

      scaleto <- if(is.null(scaleto)) nrow(egor$ego) else scaleto
      s * scaleto
    }else{
      structure(do.call(cbind, stats), order = unique(orders))
    }
}

combine_stats <- function(l){
  l <- lapply(compact(l), function(stat){
    if(inherits(stat, "svystat")) stat
    else structure(stat,
                   var=matrix(NA, length(stat), length(stat)),
                   statistic="mean", class="svystat")
  })

  lens <- lengths(l)
  starts <- cumsum(c(0, lens))
  v <- matrix(NA, sum(lens), sum(lens))
  for(i in seq_along(l)){
    v[starts[i]+seq_len(lens[i]), starts[i]+seq_len(lens[i])] <- attr(l[[i]], "var")
  }

  stats <- do.call(c, l)
  rownames(v) <- colnames(v) <- names(stats)

  structure(stats, var = v,
            statistic = "mean", class="svystat")
}

#' A scalar multiplication method for `svystat`
#'
#' Multiply the values of survey statistics by a specified vector elementwise, adjusting the variance.
#'
#' @param x an object of class `[svystat][survey::svymean]`.
#' @param y a numeric vector equal in length to `x`; shorter vectors will be recycled.
#'
#' @return a `[svystat][survey::svymean]` object with the updated statistics and variance-covariance matrix.
#'
#' @examples
#' library(survey)
#' data(api)
#' # From example(svymean):
#' dclus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
#'
#' (m1 <- svymean(~api99, dclus1))
#' (v1 <- vcov(m1))
#'
#' # Scale the suvery stat object by a factor of two:
#' (m2 <- m1 * 2)
#' (v2 <- vcov(m2))
#'
#' \dontshow{
#' stopifnot(isTRUE(all.equal(as.vector(m2), as.vector(m1)*2, check.attributes=FALSE)))
#' stopifnot(isTRUE(all.equal(v2, v1*4)))
#' }
#' @export
`*.svystat` <- function(x, y){
  if(!is.numeric(y)) stop("At this time, only multiplication of ",sQuote("svystat")," by a numeric vector is supported.")
  if(length(x) %% length(y)) warning("length of ", sQuote("x"), " is not an integer multiple of length of ", sQuote("y"))
  y <- rep_len(y, length(x))
  o <- unclass(x);
  if(!startsWith(attr(o, "statistic"), "scaled")) attr(o, "statistic") <- paste("scaled", attr(o, "statistic"))
  o <- o * y
  attr(o, "var") <- t(attr(x, "var") * y) * y
  class(o) <- class(x)
  o
}


#' @describeIn summary_formula.egor A multiplication method that takes into account which statistics are scalable.
#'
#' @param x,y see [`*.svystat`].
#'
#' @examples
#'
#' (ego.summ2 <- summary(fmh.ego ~ edges + meandeg + degree(0:2)))
#' vcov(ego.summ2)
#'
#' ego.summ2 * 2 # edges and degrees scales, meandeg doesn't
#' vcov(ego.summ2 * 2)
#'
#' @export
`*.ergm.ego_svystat` <- function(x, y){
  if(!is.numeric(y)) stop("At this time, only multiplication of ",sQuote("svystat")," by a numeric vector is supported.")
  if(length(x) %% length(y)) warning("length of ", sQuote("x"), " is not an integer multiple of length of ", sQuote("y"))
  if(length(y) == length(x) && any(y != 1 & !attr(x, "scaling"))) warning("attempting to scale a nonscalable ego statistic: scale will be ignored")
  y <- rep_len(y, length(x))
  y[!attr(x, "scaling")] <- 1
  NextMethod()
}
