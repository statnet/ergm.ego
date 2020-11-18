#  File R/summary.statistics.egodata.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################


#' Calculation of ERGM-style summary statistics for \code{\link{egodata}}
#' objects.
#' 
#' Used to calculate the specified network statistics inferred from a
#' \code{\link{egodata}} object.
#' 
#' 
#' @aliases summary_formula.egodata summary summary_formula
#' @param object An \code{\link[ergm]{ergm}}-style formula with a
#' \code{\link{egodata}} object as the LHS.
#' 
#' For a list of currently implemented egocentric terms for the RHS, see
#' \code{\link{ergm.ego-terms}}.
#' @param \dots Not used at this time.
#' @param basis An optional \code{\link{egodata}} object relative to which the
#' statistics should be calculated.
#' @param individual If \code{FALSE} (the default), calculate the estimated
#' per-capita statistics, weighted according to the ego weights, then scale
#' them up to a network of size \code{scaleto}.
#' 
#' If \code{TRUE}, calculate each ego's individual contribution to the
#' specified network statistics.
#' @param scaleto Size of a hypothetical network to which to scale the
#' statistics. Defaults to the number of egos in the dataset.
#' @return If \code{individual==FALSE}, a named vector of statistics. If
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
#' fmh.ego <- as.egodata(faux.mesa.high)
#' (nw.summ <- summary(faux.mesa.high~edges+degree(0:3)+nodematch("Race")+
#'                     nodematch("Sex")+absdiff("Grade")+nodemix("Grade")))
#' 
#' (ego.summ <- summary(fmh.ego~edges+degree(0:3)+nodematch("Race")+nodematch("Sex")+
#'                      absdiff("Grade")+nodemix("Grade"),
#'                      scaleto=network.size(faux.mesa.high)))
#' 
#' stopifnot(isTRUE(all.equal(nw.summ,ego.summ)))
#' 
#' @export
summary_formula.egodata <- function(object,..., basis=NULL, individual=FALSE, scaleto=NULL){
  egodata <-
    if(!is.null(basis)) basis
    else eval_lhs.formula(object)

  scaling.stats <- NULL
  scaling.pos <- c(0)
  nonscaling.stats <- c()
  nonscaling.pos <- c(0)
  
  
  for(trm in list_rhs.formula(object)){
    if(is.call(trm)){
      egostat <- locate_prefixed_function(trm[[1]], "EgoStat", "Egocentric statistic", env=environment(object))
      init.call <- list(egostat, egodata=egodata)
      init.call <- c(init.call,as.list(trm[-1]))
    }else{
      egostat <- locate_prefixed_function(trm, "EgoStat", "Egocentric statistic", env=environment(object))
      init.call <- list(egostat, egodata=egodata)
    }
    stat<-eval(as.call(init.call), environment(object))
    if(isTRUE(attr(stat, "nonscaling"))){
      if(individual) stop("Nonscaling statistic detected. Individual contributions are meaningless.")
      nonscaling.stats <- c(nonscaling.stats, stat)
      nonscaling.pos <- c(nonscaling.pos, max(scaling.pos,nonscaling.pos) + seq_len(length(stat)))
    }else{
      scaling.stats<-cbind(scaling.stats,stat)
      scaling.pos <- c(scaling.pos, max(scaling.pos,nonscaling.pos) + seq_len(ncol(stat)))
    }
  }
  
  if(!individual){
    if(length(scaling.stats)){
      scaleto <- if(is.null(scaleto)) nrow(egodata$egos) else scaleto
      scaling.stats <- colSums(scaling.stats*egodata$egoWt)/sum(rep(egodata$egoWt,length.out=nrow(scaling.stats)))
      scaling.stats <- scaling.stats*scaleto
    }
      
    stats <- numeric(max(scaling.pos,nonscaling.pos))
    scaling.pos <- scaling.pos[scaling.pos>0]
    nonscaling.pos <- nonscaling.pos[nonscaling.pos>0]

    stats[scaling.pos] <- scaling.stats
    stats[nonscaling.pos] <- nonscaling.stats
    
    names(stats)[scaling.pos] <- names(scaling.stats)
    names(stats)[nonscaling.pos] <- names(nonscaling.stats)
    
    stats
  }else{
    rownames(scaling.stats) <- egodata$egos[[egodata$egoIDcol]]
    scaling.stats
  }
}
