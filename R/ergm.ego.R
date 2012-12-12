ergm.ego <- function(formula, popsize, offset.coef=NULL, ..., control=control.ergm()){
  control$force.main <- TRUE
  egodata <- get(as.character(formula[[2]]), envir=environment(formula))
  
  popnw <- as.network.egodata(x, popsize)

  # Get the sample h values.
  stats <- summary(remove.offset.formula(formula), individual=TRUE)

  # We can have summary() compute some of this, but that's likely to
  # be slow.
  w <- egodata$egoWt/sum(egodata$egoWt)
  m <- colSums(stats*w)
  # TODO: Include finite-population correction here:
  v <- crossprod(sweep(stats, 2, m, "-")*sqrt(w))/(1-sum(w^2))

  ergm.formula <- ergm.update.formula(formula,popnw~offset(edges)+.)

  ergm.fit <- ergm(ergm.formula, target.stats=m*popsize, offset.coef=c(-log(popsize),offset.coef),..., eval.loglik=FALSE,control=control)
  coef <- coef(ergm.fit)

  oi <- offset.info.formula(ergm.formula)

  DtDe <- -ergm.fit$hessian[!oi$theta,!oi$theta,drop=FALSE]

  vcov <- matrix(NA, length(coef), length(coef))
  
  
  vcov[!oi$theta,!oi$theta] <- solve(DtDe/popsize)%*%v%*%solve(DtDe/popsize)
  
  list(coef=coef, vcov=vcov)
}
