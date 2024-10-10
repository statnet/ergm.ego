#  File R/gof.ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################
# This file contains the following 8 functions for assessing goodness of fit
#         <gof>              <summary.gofobject>
#         <gof.default>      <plot.gofobject>
#         <gof.ergm>         <ergm.get.terms.formula>
#         <gof.formula>      <ergm.rhs.formula>
#=============================================================================



###############################################################################
# Each of the <gof.X> functions assesses the goodness of fit of X by comparison
# with 'control$nsim' ergm simulations of X
#
# --PARAMETERS--
#   object/formula: either an ergm object or a formula
#   ...           : additional parameters passed from within the program;
#                   these are ignored
#   init        : the parameters from which the simulations will be drawn;
#                   default=NULL;
#   control$nsim          : the number of simulated ergms, with which to compare X;
#                   default=100
#   burnin        : the number of proposals to disregard before any MCMC
#                   sampling is done; this is passed along to the simulation
#                   routines; default=10000
#   interval      : the number of proposals between sampled ergm statistics;
#                   this is passed along to the simulation rountines;
#                   default=1000
#   GOF           : a one-sided formula specifying which summary statistics
#                   should be used in the GOF comparison; choices include
#                       distance      espartners    dspartners
#                       odegree       idegree       degree
#                       triadcensus   model
#                   default=NULL; is internally mapped to 
#                   ~degree+espartners+distance if nw is undirected, and
#                   ~idegree+odegree+espartners+distance otherwise
#   constraints   : a one-sided formula of the constraint terms; options are
#                         bd        degrees        nodegrees
#                         edges     degreedist     idegreedist
#                         observed  odegreedist
#                   default="~ ."   
#   control       : a list of parameters for controlling GOF evaluation, as
#                   returned by <control.gof.X>; default=control.gof.X()
#                   (note that <control.gof.X> has different defaults 
#                    depending on the class of X)
#   seed          : an integer value at which to set the random generator;
#                   default=NULL
#   verbose       : whether to print information on the progress of the
#                   simulations; default=FALSE
#
# --RETURNED--
#   returnlist: a list with the following components for each term
#               G given in 'GOF'
#      summary.G: a matrix of summary statistics for the observed and
#                 simulated G's; if G takes on the values {G1, G2,...,Gq},
#                 the entries of 'summary.G' are
#         [i,1]-- the observed frequency of Gi
#         [i,2]-- the minimum value of Gi from the simulations
#         [i,3]-- the mean value of Gi from the simulations
#         [i,4]-- the maximum value of Gi from the simulations
#         [i,5]-- the p-value for the observed Gi estimated from the
#                 distribution of simulations
#      pobs.G   : a vector giving G's observed probability distribution
#      psim.G   : a matrix of G's simulated probability distributions; each
#                 row gives a distribution
#      bds.G    : the estimatd confidence interval, as the .025 and .975
#                 quantile values of the simulations
#      obs.G    : the vector of summary statistics for the observed X
#      sim.G    : the matrix of summary statistics for each simulated
#                 version of X
#
###############################################################################



#' Conduct Goodness-of-Fit Diagnostics on a Exponential Family Random Graph
#' Model fit to Egocentrically Sampled Data
#' 
#' [gof.ergm.ego()] implements the [gof()] method for
#' [`ergm.ego`] fit objects.
#' 
#' 
#' @param object An [`ergm.ego`] fit.
#' @param ... Additional arguments. Unused by [gof.ergm.ego()], passed to
#'   [ergm::plot.gof()] by [plot.gof.ergm.ego()]
#' @param GOF A string specifying the statistics whose goodness of fit is to be
#' evaluated. Currently, only \dQuote{degree}, \dQuote{espartners} and \dQuote{model} are
#' implemented; see [gof()] documentation for details.
#' @param control A list to control parameters, constructed using
#' [control.gof.formula()] or [control.gof.ergm()] (which
#' have different defaults).
#' @param verbose Provide verbose information on the progress of the
#' simulation.
#' @return An object of class [`gof.ergm.ego`], inheriting from [gof.ergm()].
#' @author Pavel N. Krivitsky
#' @seealso For examples, see [ergm.ego()].
#' @keywords models
#' @examples
#' \donttest{
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#' 
#' head(fmh.ego)
#' 
#' egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'                          +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
#'                           popsize=network.size(faux.mesa.high))
#' 
#' # Check whether the model "converged":
#' (modelgof <- gof(egofit, GOF="model"))
#' plot(modelgof)
#' 
#' # Check whether the model reconstructs the degree distribution:
#' (deggof <- gof(egofit, GOF="degree"))
#' plot(deggof)
#' }
#' @import ergm stats
#' @method gof ergm.ego
#'
#' @references
#'
#' * David R. Hunter, Steven M. Goodreau, and Mark S. Handcock (2008). "Goodness of Fit of Social Network Models." *Journal of the American Statistical Association*, 103:481: 248–258. \doi{10.1198/016214507000000446}
#'
#' @export
gof.ergm.ego <- function (object, ..., 
                          GOF=c("model","degree", "espartners"), 
                          control=control.gof.ergm(),
                          verbose=FALSE) {
  statnet.common::check.control.class(c("gof.ergm","gof.formula"), "gof.ergm.ego")
  
  #Set up the defaults, if called with GOF==NULL
  if(is.null(GOF)){
    GOF<- ~degree
  }

  egor <- object$egor
  formula <- object$formula
  
  control.transfer <- c("MCMC.burnin", "MCMC.prop", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

  # Rescale the interval by the ratio between the estimation sample size and the GOF sample size so that the total number of MCMC iterations would be about the same.
  NVL(control$MCMC.interval) <- max(ceiling(object$control$MCMC.interval*object$control$MCMC.samplesize/control$nsim),1)

  constraints <- object$constraints
  
  if(!is.null(control$seed)) {set.seed(as.integer(control$seed))}
  if (verbose) 
    cat("Starting GOF for the given ERGM formula.\n")

  GOF <- match.arg(GOF)
  
  # Calculate network statistics for the observed graph
  # Set up the output arrays of sim variables
  if(verbose)
    cat("Calculating observed network statistics.\n")

  n <- network.size(object$newnetwork)

  do.sim <- function(mon=NULL){
    sim <- simulate(object, nsim=control$nsim, seed=control$seed, popsize=object$ppopsize, control=control.simulate.ergm.ego(simulate=set.control.class("control.simulate.formula",control)),...,verbose=verbose, output="stats", monitor=mon)
    sim[, if(is.null(mon)) TRUE else attr(sim,"monitored"), drop=FALSE]/n
  }
  
  if ('model' %in% GOF) {
    obs.model <- summary(formula, scaleto=1, basis=egor)
    sim.model <- do.sim()
  }

  if ('degree' %in% GOF) {
    maxdeg <- max(.degreeseq(egor),3)*2
    f <- as.formula(if(maxdeg >= n-1) paste0("~degree(0:",n-1,")")
                    else paste0("~degree(0:",maxdeg-1,")+degrange(",maxdeg,")"))
    obs.deg <- summary(f, scaleto=1, basis=egor)
    sim.deg <- do.sim(f)
  }
  
  if("espartners" %in% GOF) {
    # Maximum esp = max(egodegree) + 1 - 2
    maxdeg <- max(.degreeseq(egor), 3)
    maxesp <- (maxdeg - 1)*2
    f <- as.formula(paste0("egor ~ esp(0:", maxesp, ")"))
    obs.esp <- summary(f, scaleto = 1, basis=egor)
    sim.esp <- do.sim(f)
  }

  if(verbose)
    cat("Starting simulations.\n")

  # calculate p-values
  
  returnlist <- list(network.size=n, GOF=as.formula(paste0("~",GOF)))
  
  if ('model' %in% GOF) {
    pval.model <- apply(sim.model <= obs.model[col(sim.model)],2,mean)
    pval.model.top <- apply(sim.model >= obs.model[col(sim.model)],2,mean)
    q.model <- (pval.model+1-pval.model.top)/2 # Handle ties by averaging.
    pval.model <- cbind(obs.model,apply(sim.model, 2,min), apply(sim.model, 2,mean),
                        apply(sim.model, 2,max), pmin(1,2*pmin(pval.model,pval.model.top)))
    dimnames(pval.model)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.model <- q.model
    psim.model <- apply(sim.model,2,rank)/nrow(sim.model)
    bds.model <- apply(psim.model,2,quantile,probs=c(0.025,0.975))

    returnlist$summary.model <- returnlist$pval.model <- pval.model
    returnlist$pobs.model <- pobs.model
    returnlist$psim.model <- psim.model
    returnlist$bds.model <- bds.model
    returnlist$obs.model <- obs.model
    returnlist$sim.model <- sim.model
  }

  if ('degree' %in% GOF) {
    pval.deg <- apply(sim.deg <= obs.deg[col(sim.deg)],2,mean)
    pval.deg.top <- apply(sim.deg >= obs.deg[col(sim.deg)],2,mean)
    pval.deg <- cbind(obs.deg,apply(sim.deg, 2,min), apply(sim.deg, 2,mean),
                      apply(sim.deg, 2,max), pmin(1,2*pmin(pval.deg,pval.deg.top)))
    dimnames(pval.deg)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.deg <- obs.deg/sum(obs.deg)
    psim.deg <- sweep(sim.deg,1,apply(sim.deg,1,sum),"/")
    psim.deg[is.na(psim.deg)] <- 1
    bds.deg <- apply(psim.deg,2,quantile,probs=c(0.025,0.975))
    
    returnlist$summary.deg <- returnlist$pval.deg <- pval.deg
    returnlist$pobs.deg <- pobs.deg
    returnlist$psim.deg <- psim.deg
    returnlist$bds.deg <- bds.deg
    returnlist$obs.deg <- obs.deg
    returnlist$sim.deg <- sim.deg
  }
  
  if("espartners" %in% GOF) {
    pval.esp <- apply(sim.esp <= obs.esp[col(sim.esp)], 2, mean)
    pval.esp.top <- apply(sim.esp >= obs.esp[col(sim.esp)], 2, mean)
    pval.esp <- cbind(
      obs.esp,
      apply(sim.esp, 2, min), 
      apply(sim.esp, 2, mean),
      apply(sim.esp, 2, max), 
      pmin(1, 2*pmin(pval.esp, pval.esp.top))
    )
    dimnames(pval.esp)[[2]] <- c("obs","min","mean","max","MC p-value")
    pobs.esp <- obs.esp / sum(obs.esp)
    psim.esp <- sweep(sim.esp, 1 ,apply(sim.esp, 1, sum), "/")
    psim.esp[is.na(psim.esp)] <- 1
    bds.esp <- apply(psim.esp, 2, quantile,probs=c(0.025,0.975))
    
    returnlist$summary.espart <- returnlist$pval.espart <- pval.esp
    returnlist$pobs.espart <- pobs.esp
    returnlist$psim.espart <- psim.esp
    returnlist$bds.espart <- bds.esp
    returnlist$obs.espart <- obs.esp
    returnlist$sim.espart <- sim.esp
    
  }

  class(returnlist) <- c("gof.ergm.ego", "gof.ergm", "gof")
  returnlist
}

#' @rdname gof.ergm.ego
#' @description An enhanced plotting method is also provided, giving uncertainty bars for the observed statistics as well.
#' @param x an object returned by [gof.ergm.ego()].
#' @param ego.conf.level confidence level for the observed statistic estimates as well.
#' @method plot gof.ergm.ego
#' @importFrom graphics lines
#' @export
plot.gof.ergm.ego <- function(x, ..., ego.conf.level=0.95){
  # Call the plotting method for gof objects in general.
  NextMethod()
  # Add confidence intervals based on some `obs` stats
  add_obs_ci <- function(obs, conf.level) {
    obs.CI <- confint(obs, level=conf.level) # svydesign() summary objects have the necessary components for a confint() method.
    xvals <- seq_len(nrow(obs.CI))
    lines(xvals, obs.CI[,1], lty=3)
    lines(xvals, obs.CI[,2], lty=3)
    obs.CI
  }
  # Now, overplot with our statistics if degrees.
  if(!is.null(x$obs.deg)) {
    invisible(add_obs_ci(x$obs.deg, conf.level = ego.conf.level))
  } 
  if(!is.null(x$obs.espart)) {
    invisible(add_obs_ci(x$obs.espart, conf.level = ego.conf.level))
  }
  invisible()
}
