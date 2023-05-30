#  File R/simulate.ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2023 Statnet Commons
################################################################################


#' Simulate from a \code{\link{ergm.ego}} fit.
#' 
#' A wrapper around \code{\link[ergm]{simulate.formula}} to simulate networks
#' from an ERGM fit using \code{\link{ergm.ego}}.
#' 
#' 
#' @param object An \code{\link{ergm.ego}} fit.
#' @param nsim Number of realizations to simulate.
#' @template seed
#' @param popsize,basis A network size to which to scale the model for
#'   simulation; a [`data.frame`] with at least those ego attributes
#'   used to estimate the model to simulate over a specific set of
#'   actors; or a [`network`] object to use as is. `basis` is provided
#'   for consistency with [ergm()], [ergm.ego()], [simulate.ergm()],
#'   and others. If both are specified, `popsize` overrules.
#' @param control A \code{\link{control.simulate.ergm.ego}} control list.
#' @param output one of `"network"`, `"stats"`, `"edgelist"`,
#'   `"pending_update_network"`, or, for future compatibility,
#'   `"ergm_state"`. See help for [simulate.ergm()] for explanation.
#' 
#' @param constraints,\dots Additional arguments passed to \code{\link[ergm]{san}} and
#' \code{\link[ergm]{simulate.formula}}.
#' @template verbose
#' @return The ouput has the same format (with the same options) as
#' \code{\link[ergm]{simulate.formula}}. If \code{output="stats"} is passed, an
#' additional attribute, \code{"ppopsize"} is set, giving the actual size of
#' the network reconstructed, when the \code{pop.wt} control parameter is set
#' to \code{"round"} and \code{"popsize"} is not a multiple of the egocentric
#' sample size or the sampling weights.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link[ergm]{simulate.formula}},
#' \code{\link[ergm]{simulate.ergm}}
#' @references
#'
#' * Pavel N. Krivitsky and Martina Morris (2017). "Inference for social network models from egocentrically sampled data, with application to understanding persistent racial disparities in HIV prevalence in the US." *Annals of Applied Statistics*, 11(1): 427–455. \doi{10.1214/16-AOAS1010}
#'
#' * Pavel N. Krivitsky, Martina Morris, and Michał Bojanowski (2019). "Inference for Exponential-Family Random Graph Models from Egocentrically-Sampled Data with Alter–Alter Relations." NIASRA Working Paper 08-19. \url{https://www.uow.edu.au/niasra/publications/}
#' 
#' * Pavel N. Krivitsky, Mark S. Handcock, and Martina Morris (2011). "Adjusting for
#' Network Size and Composition Effects in Exponential-Family Random Graph
#' Models." \emph{Statistical Methodology}, 8(4): 319–339. \doi{10.1016/j.stamet.2011.01.005}
#'
#' @keywords models
#' @examples
#' 
#' data(faux.mesa.high)
#' data(fmhfit)
#' colMeans(egosim <- simulate(fmhfit, popsize=300,nsim=50,
#'                        output="stats", control=control.simulate.ergm.ego(
#'                        simulate=control.simulate.formula(MCMC.burnin=2e6))))
#' colMeans(egosim)/attr(egosim,"ppopsize")*network.size(faux.mesa.high)
#' summary(faux.mesa.high~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'            +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"))
#'
#' @importFrom stats simulate
#' @export
simulate.ergm.ego <- function(object, nsim = 1, seed = NULL, constraints=object$constraints, popsize=if(object$popsize==1 || object$popsize==0 || is(object$popsize, "AsIs")) object$ppopsize else object$popsize, control=control.simulate.ergm.ego(), output=c("network","stats","edgelist","pending_update_network", "ergm_state"), ..., basis=NULL, verbose=FALSE){
  statnet.common::check.control.class("simulate.ergm.ego", "simulate.ergm.ego")
  output <- match.arg(output)

  nsa <- !is.null(object$netsize.adj)

  if(!is.null(basis)) popsize <- basis
  if(is.network(popsize)){
    popnw <- popsize
    popsize <- network.size(popsize)
  }else if(is.data.frame(popsize)){ # If pseudopoluation composition is given in popsize, use that.
    popnw <- template_network(popsize, nrow(popsize))
    popsize <- nrow(popsize)
  }else{
    popnw <-
      if(popsize == object$ppopsize) object$newnetwork
      else template_network(object$egor, popsize, scaling=control$ppop.wt)
  }

  if(!nsa && network.size(popnw)!=object$ppopsize) warning("Simulation network size ", network.size(popnw), " different from fitted ", object$ppopsize,", with network size adjustment disabled. Simulation may not be meaningful.", immediate.=TRUE)

  ppopsize <- if(network.size(popnw)!=popsize){
    message("Note: Constructed network has size ", network.size(popnw), " different from requested ", popsize,". Simulated statistics may need to be rescaled.")
    network.size(popnw)
  }else popsize

  # TODO: Make it work with ergm_state output.
  if(popsize != object$ppopsize) popnw <- san(object$formula, target.stats = object$target.stats/object$ppopsize*ppopsize,verbose=verbose, constraints=constraints, basis=popnw, control=control$SAN, ..., output="network")

  ergm.formula <- if(nsa) nonsimp_update.formula(object$formula,object$netsize.adj) else object$formula

  out <- simulate(ergm.formula, nsim=nsim, seed=seed, verbose=verbose,
                  coef=c(netsize.adj=if(nsa) -log(ppopsize/object$popsize),
                         coef(object)[if(nsa) -1 else TRUE]),
                  constraints=constraints, control=control$simulate, basis=popnw, ..., output=output)
  if(is.matrix(out)){
    if(nsa)
      out <- structure(out[,-1, drop=FALSE],
                       monitored = attr(out, "monitored")[-1],
                       formula = attr(out,"formula"), .Basis = attr(out,".Basis"), monitor = attr(out,"monitor"),
                       constraints = attr(out,"constraints"), reference = attr(out,"reference"))

    attr(out, "ppopsize") <- ppopsize
  }
  out
}
