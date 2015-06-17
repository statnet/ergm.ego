simulate.ergm.ego <- function(object, nsim = 1, seed = NULL, popsize=if(object$popsize==1) object$ppopsize else object$popsize, control=control.simulate.ergm.ego(), ..., verbose=FALSE){
  check.control.class()
  
  egodata <- object$egodata
  popnw <- as.network(egodata, popsize, scaling=control$ppop.wt)
  ppopsize <- if(network.size(popnw)!=popsize){
    message("Note: Constructed network has size ", network.size(popnw), " different from requested ", popsize,". Simulated statistics may need to be rescaled.")
    network.size(popnw)
  }else popsize

  
  ergm.formula <- ergm.update.formula(object$formula,popnw~.,from.new="popnw")

  popnw <- san(ergm.formula, target.stats = object$target.stats[-1]/object$ppopsize*ppopsize,verbose=verbose, control=control$SAN.control, ...)
  ergm.formula <- ergm.update.formula(object$formula,popnw~netsize.adj+.,from.new="popnw")

  out <- simulate(ergm.formula, nsim=nsim, seed=seed, verbose=verbose, coef=c(netsize.adj=-log(ppopsize/object$popsize),object$coef[-1]), control=control$simulate.control, ...)
  if(is.matrix(out)){
    out <- out[,-1,drop=FALSE]
    attr(out, "ppopsize") <- ppopsize
  }
  out
}
