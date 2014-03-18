ergm.ego <- function(formula, popsize, ppopsize=popsize, offset.coef=NULL, na.action=na.fail, stats.var = c("ratio", "naive"), ..., control=control.ergm()){
  control$force.main <- TRUE
  stats.var <- match.arg(stats.var)
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
  # TODO: Jacknife for bias-correction and variance estimation?
  w <- w/sum(w)
  m <- colSums(stats*w)
  # TODO: Include finite-population correction here:
  v <- switch(stats.var,
              naive = crossprod(sweep(stats, 2, m, "-")*sqrt(w))/(1-sum(w^2))*sum(w^2)/sum(w)^2,
              ratio = .ratio.var(stats, w)/length(w)
              )
  
  ergm.formula <- ergm.update.formula(formula,popnw~offset(edges)+.,from.new="popnw")

  ergm.fit <- ergm(ergm.formula, target.stats=m*ppopsize, offset.coef=c(-log(ppopsize/popsize),offset.coef),..., eval.loglik=FALSE,control=control)
  coef <- coef(ergm.fit)

  oi <- offset.info.formula(ergm.formula)

  DtDe <- -ergm.fit$hessian[!oi$theta,!oi$theta,drop=FALSE]

  vcov <- matrix(NA, length(coef), length(coef))
  
  vcov[!oi$theta,!oi$theta] <- solve(DtDe/ppopsize)%*%v%*%solve(DtDe/ppopsize)

  rownames(vcov) <- colnames(vcov) <- names(coef)
  
  out <- list(coef=coef[-1], netsize.offset=coef[1], vcov=vcov[-1,-1,drop=FALSE], ergm.fit=ergm.fit, DtDe=DtDe, v=v, formula=formula, egodata=egodata, ppopsize=ppopsize, popsize=popsize)
  class(out) <- "ergm.ego"
  out
}

.ratio.var <- function(X, w){
  n <- length(w)
  p <- ncol(X)
  wX <- cbind(w,sweep(X,1,w,"*"))
  m <- colSums(wX[,-1,drop=FALSE])/sum(w)
  S <- cov(wX)
  A <- 1/mean(w)*cbind(-m,diag(1,nrow=p))
  A%*%S%*%t(A)
}
