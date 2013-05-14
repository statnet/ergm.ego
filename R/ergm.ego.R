ergm.ego <- function(formula, popsize, ppopsize=popsize, offset.coef=NULL, na.action=na.fail, ..., control=control.ergm()){
  control$force.main <- TRUE
  egodata <- get(as.character(formula[[2]]), envir=environment(formula))
  
  popnw <- as.network(egodata, ppopsize)

  # Get the sample h values.
  stats <- summary(remove.offset.formula(formula), individual=TRUE)

  # h is just a matrix, so this will do the sensible thing.
  w <- egodata$egoWt
  tmp <- na.action(cbind(w,stats))
  w <- tmp[,1]
  stats <- tmp[,-1, drop=FALSE]

  # We can have summary() compute some of this, but that's likely to
  # be slow.
  w <- w/sum(w)
  m <- colSums(stats*w)
  # TODO: Include finite-population correction here:
  v <- crossprod(sweep(stats, 2, m, "-")*sqrt(w))/(1-sum(w^2))

  ergm.formula <- ergm.update.formula(formula,popnw~offset(edges)+.,from.new="popnw")

  ergm.fit <- ergm(ergm.formula, target.stats=m*ppopsize, offset.coef=c(-log(ppopsize/popsize),offset.coef),..., eval.loglik=FALSE,control=control)
  coef <- coef(ergm.fit)

  oi <- offset.info.formula(ergm.formula)

  DtDe <- -ergm.fit$hessian[!oi$theta,!oi$theta,drop=FALSE]

  vcov <- matrix(NA, length(coef), length(coef))
  
  vcov[!oi$theta,!oi$theta] <- solve(DtDe/ppopsize)%*%v%*%solve(DtDe/ppopsize)*sum(w^2)/sum(w)^2

  rownames(vcov) <- colnames(vcov) <- names(coef)
  
  list(coef=coef, vcov=vcov, ergm.fit=ergm.fit, DtDe=DtDe, v=v)
}
