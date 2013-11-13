simulate.ergm.ego <- function(object, nsim = 1, seed = NULL, popsize=object$popsize, control=control.simulate.ergm.ego(), ..., verbose=FALSE){
  egodata <- object$egodata
  popnw <- as.network(egodata, popsize)
  ergm.formula <- ergm.update.formula(object$formula,popnw~.,from.new="popnw")

  popnw <- san(ergm.formula, target.stats = object$ergm.fit$target.stats[-1]/object$popsize*popsize,verbose=verbose, control=control$SAN.control, ...)
  ergm.formula <- ergm.update.formula(object$formula,popnw~edges+.,from.new="popnw")

  out <- simulate(ergm.formula, nsim=nsim, seed=seed, verbose=verbose, coef=c(-log(popsize/object$popsize),object$coef), control=control$simulate.control, ...)
  if(is.matrix(out)) out <- out[,-1,drop=FALSE]
}
