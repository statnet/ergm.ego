#  File R/summary.statistics.egor.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################
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
#' @return If \code{individual==FALSE}, a [svystat][survey::svymean] object---effectively a named vector of statistics. If
#' \code{individual==TRUE}, a matrix with a row for each ego, giving that ego's
#' contribution to the network statistic.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link[ergm]{summary_formula}},
#' \code{\link[ergm]{summary_formula.ergm}}
#' @references Pavel N. Krivitsky and Martina Morris. Inference for Social
#' Network Models from Egocentrically-Sampled Data, with Application to
#' Understanding Persistent Racial Disparities in HIV Prevalence in the US.
#' Thechnical Report. National Institute for Applied Statistics Research
#' Australia, University of Wollongong, 2015(05-15).
#' \doi{10.1214/16-AOAS1010}
#' 
#' Pavel N. Krivitsky, Mark S. Handcock, and Martina Morris. Adjusting for
#' Network Size and Composition Effects in Exponential-Family Random Graph
#' Models. \emph{Statistical Methodology}, 2011, 8(4), 319-339.
#' \doi{10.1016/j.stamet.2011.01.005}
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
    else get(as.character(object[[2]]), envir=environment(object))

  scaling.stats <- NULL
  scaling.pos <- c(0)
  nonscaling.stats <- c()
  nonscaling.pos <- c(0)
  orders <- c()
  
  for(trm in list_rhs.formula(object)){
    if(is.call(trm)){
      init.call <- list(as.name(paste("EgoStat.", trm[[1]],sep="")),egor=egor)
      init.call <- c(init.call,as.list(trm[-1]))
    }else{
      init.call <- list(as.name(paste("EgoStat.", trm,sep="")),egor=egor)
    }
    stat <- eval(as.call(init.call), environment(object))
    if(attr(stat, "order")==0){
      if(individual) stop("Nonscaling statistic detected. Individual contributions are meaningless.")
      nonscaling.stats <- c(nonscaling.stats, stat)
      nonscaling.pos <- c(nonscaling.pos, max(scaling.pos,nonscaling.pos) + seq_len(length(stat)))
    }else{
      orders <- c(orders, attr(stat, "order"))
      scaling.stats<-cbind(scaling.stats,stat)
      scaling.pos <- c(scaling.pos, max(scaling.pos,nonscaling.pos) + seq_len(ncol(stat)))
    }
  }
  
  stats <-
    if(!individual){
      if(length(scaling.stats)){
        scaleto <- if(is.null(scaleto)) nrow(egor) else scaleto
        scaling.stats <- svymean(scaling.stats, ego_design(egor), ...)
        scaling.stats <- scaling.stats*scaleto
      }
      
      stats <- numeric(max(scaling.pos,nonscaling.pos))
      scaling.pos <- scaling.pos[scaling.pos>0]
      nonscaling.pos <- nonscaling.pos[nonscaling.pos>0]
      
      stats[scaling.pos] <- scaling.stats
      stats[nonscaling.pos] <- nonscaling.stats
      
      names(stats)[scaling.pos] <- names(scaling.stats)
      names(stats)[nonscaling.pos] <- names(nonscaling.stats)

      attr(stats, "var") <- matrix(NA, max(scaling.pos,nonscaling.pos), max(scaling.pos,nonscaling.pos))
      attr(stats, "var")[scaling.pos,scaling.pos] <- attr(scaling.stats, "var")
      dimnames(attr(stats, "var")) <- rep(list(names(stats)), 2)

      attr(stats, "statistic") <- paste("scaled mean")
      
      class(stats) <- class(scaling.stats)
      
      stats
    }else{
      scaling.stats
    }

  attr(stats,"order") <- unique(orders)
  stats
}

#' A scalar multiplication method for `svystat`
#'
#' Multiply the values of survey statistics by the specified number, adjusting the variance.
#'
#' @param x an object of class `[svystat][survey::svymean]`.
#' @param y a scalar (numeric vector of length 1).
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
  if(!is.numeric(y) || length(y)!=1) stop("At this time, only scalar multiplication of ",sQuote("svystat")," objects is supported.")
  o <- unclass(x);
  attr(o, "statistic") <- paste("scaled", attr(o, "statistic"))
  o <- o * y
  attr(o, "var") <- attr(x, "var") * y^2
  class(o) <- class(x)
  o
}
