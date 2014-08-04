ergm.ego <- function(formula, popsize, ppopsize=popsize, offset.coef=NULL, na.action=na.fail, stats.est = c("asymptotic", "bootstrap", "jackknife", "naive"), stats.wt = c("data","ppop"), ppop.wt = "round", R=10000, ..., control=control.ergm(), do.fit=TRUE){
  stats.est <- match.arg(stats.est)
  stats.wt <- match.arg(stats.wt)
  egodata <- get(as.character(formula[[2]]), envir=environment(formula))

  message("Constructiong pseudopopulation network.")
  popnw <- as.network(egodata, ppopsize, scaling=ppop.wt)
  if(network.size(popnw)!=ppopsize){
    message("Note: Constructed network has size ", network.size(popnw), ", different from requested ", ppopsize,". Estimation should not be meaningfully affected.")
    ppopsize <- network.size(popnw)
  }
  

  w <- switch(stats.wt,
              data=egodata$egoWt,
              ppop=tabulate(popnw %n% "ego.inds", nbins=nrow(egodata))
              )
  
  # Get the sample h values.
  stats <- summary(remove.offset.formula(formula), individual=TRUE)

  # h is just a matrix, so this will do the sensible thing.
  tmp <- na.action(cbind(w,stats))
  w <- tmp[,1]
  stats <- tmp[,-1, drop=FALSE]
  n <- length(w)
  
  wmean <- function(w,s){
    w <- w/sum(w)
    colSums(s*w)
  }
  m <- wmean(w,stats)
  
  if(stats.est=="bootstrap"){
    m.b <- t(replicate(R,{
      i <- sample.int(length(w),replace=TRUE)
      wmean(w[i],stats[i,,drop=FALSE])
    }))
    m <- m - (colMeans(m.b)-m)
    
  }else if(stats.est=="jackknife"){
    m.j <- t(sapply(seq_len(n), function(i){
      wmean(w[-i],stats[-i,,drop=FALSE])
    }))
    m <- n*m - (n-1)*colMeans(m.j)
  }    
  
  # TODO: Include finite-population correction here:
  v <- switch(stats.est,
              bootstrap = cov(m.b),
              jackknife = (n-1)/n*crossprod(sweep(m.j,2,colMeans(m.j))),
              naive = {w <- w/sum(w); crossprod(sweep(stats, 2, m, "-")*sqrt(w))/(1-sum(w^2))*sum(w^2)},
              asymptotic = .asymptotic.var(stats, w)/length(w)
              )

  out <- list(v=v, m=m, formula=formula, egodata=egodata, ppopsize=ppopsize, popsize=popsize)
  
  if(do.fit){
    ergm.formula <- ergm.update.formula(formula,popnw~offset(edges)+.,from.new="popnw")

    ergm.fit <- ergm(ergm.formula, target.stats=m*ppopsize, offset.coef=c(-log(ppopsize/popsize),offset.coef),..., eval.loglik=FALSE,control=control)
    coef <- coef(ergm.fit)

    oi <- offset.info.formula(ergm.formula)

    DtDe <- -ergm.fit$hessian[!oi$theta,!oi$theta,drop=FALSE]

    vcov <- matrix(NA, length(coef), length(coef))
  
    vcov[!oi$theta,!oi$theta] <- solve(DtDe/ppopsize)%*%v%*%solve(DtDe/ppopsize)

    rownames(vcov) <- colnames(vcov) <- names(coef)
  
    out <- c(out, list(coef=coef[-1], netsize.offset=coef[1], vcov=vcov, ergm.fit=ergm.fit, DtDe=DtDe))
  }
  class(out) <- "ergm.ego"
  out
}

.asymptotic.var <- function(X, w){
  n <- length(w)
  p <- ncol(X)
  wX <- cbind(w,sweep(X,1,w,"*"))
  m <- colSums(wX[,-1,drop=FALSE])/sum(w)
  S <- cov(wX)
  A <- 1/mean(w)*cbind(-m,diag(1,nrow=p))
  A%*%S%*%t(A)
}

