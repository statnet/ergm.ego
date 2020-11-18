#  File R/simulate.ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################


#' Simulate from a \code{\link{ergm.ego}} fit.
#' 
#' A wrapper around \code{\link[ergm]{simulate.formula}} to simulate networks
#' from an ERGM fit using \code{\link{ergm.ego}}.
#' 
#' 
#' @param object An \code{\link{ergm.ego}} fit.
#' @param nsim Number of realizations to simulate.
#' @param seed Random seed.
#' @param popsize Either network size to which to scale the model for
#'   simulation or a [`data.frame`] with at least those ego attributes
#'   required to estimate the model, to simulate over a specific set
#'   of actors.
#' @param control A \code{\link{control.simulate.ergm.ego}} control list.
#' 
#' @param output one of `"network"`, `"stats"`, `"edgelist"`,
#'   `"pending_update_network"`, or, for future compatibility,
#'   `"ergm_state"`. See help for [simulate.ergm()] for explanation.
#' 
#' @param constraints,\dots Additional arguments passed to \code{\link[ergm]{san}} and
#' \code{\link[ergm]{simulate.formula}}.
#' @param verbose Verbosity of output.
#' @return The ouput has the same format (with the same options) as
#' \code{\link[ergm]{simulate.formula}}. If \code{output="stats"} is passed, an
#' additional attribute, \code{"ppopsize"} is set, giving the actual size of
#' the network reconstructed, when the \code{pop.wt} control parameter is set
#' to \code{"round"} and \code{"popsize"} is not a multiple of the egocentric
#' sample size or the sampling weights.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link[ergm]{simulate.formula}},
#' \code{\link[ergm]{simulate.ergm}}
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
#' @keywords models
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egodata(faux.mesa.high)
#' egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'                +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
#'                popsize=network.size(faux.mesa.high))
#' colMeans(egosim <- simulate(egofit, popsize=300,nsim=50,
#'                        output="stats", control=control.simulate.ergm.ego(
#'                        simulate.control=control.simulate.formula(MCMC.burnin=2e6))))
#' colMeans(egosim)/attr(egosim,"ppopsize")*network.size(faux.mesa.high)
#' summary(faux.mesa.high~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'            +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"))
#'
#' @importFrom stats simulate
#' @export
simulate.ergm.ego <- function(object, nsim = 1, seed = NULL, constraints=object$constraints, popsize=if(object$popsize==1) object$ppopsize else object$popsize, control=control.simulate.ergm.ego(), output=c("network","stats","edgelist","pending_update_network", "ergm_state"), ..., verbose=FALSE){
  statnet.common::check.control.class("simulate.ergm.ego", "simulate.ergm.ego")
  output <- match.arg(output)
  
  egodata <- object$egodata
  if(is.data.frame(popsize)){ # If pseudopoluation composition is given in popsize, use that.
    pegos <- popsize
    pegos[[".pegoID"]] <- 1:nrow(pegos)
    pdata <- egodata(pegos, data.frame(.pegoID=c()), egoIDcol = ".pegoID")
    popsize <- nrow(popsize)
    popnw <- as.network.egodata(pdata, popsize)
  }else{
    popnw <-
      if(popsize == object$ppopsize) object$newnetwork
      else as.network(egodata, popsize, scaling=control$ppop.wt)
  }

  ppopsize <- if(network.size(popnw)!=popsize){
      message("Note: Constructed network has size ", network.size(popnw), " different from requested ", popsize,". Simulated statistics may need to be rescaled.")
    network.size(popnw)
  }else popsize

  san.stats <-
    if(length(object$target.stats)>nparam(object, offset=FALSE)) object$target.stats[!object$etamap$offsettheta]
    else object$target.stats
  if(popsize != object$ppopsize) popnw <- san(object$formula, target.stats = san.stats/object$ppopsize*ppopsize,verbose=verbose, constraints=constraints, basis=popnw, control=control$SAN.control, ..., output=if(utils::packageVersion("ergm")>="4") "ergm_state" else "pending_update_network")
  ergm.formula <- nonsimp_update.formula(object$formula,.~offset(netsize.adj)+.)

  out <- simulate(ergm.formula, nsim=nsim, seed=seed, verbose=verbose, coef=c(netsize.adj=-log(ppopsize/object$popsize),object$coef[-1]), constraints=constraints, control=control$simulate.control, basis=popnw, ..., output=output)
  if(is.matrix(out)){
    out <- out[,-1,drop=FALSE]
    attr(out, "ppopsize") <- ppopsize
  }
  out
}
