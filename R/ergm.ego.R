#  File R/ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2016 Statnet Commons
#######################################################################


#' Inference for Exponential-Family Random Graph Models based on Egocentrically
#' Sampled Data
#' 
#' A wrapper around the \code{\link[ergm]{ergm}} to fit an ERGM to an
#' \code{\link{egor}}.
#' 
#' 
#' @param formula An \code{\link{formula}} object, of the form \code{e
#'   ~ <model terms>}, where \code{e} is a \code{\link{egor}}
#'   object. See \code{\link[ergm]{ergm}} for details and examples.
#' 
#' For a list of currently implemented egocentric terms for the RHS, see
#' \code{\link{ergm.ego-terms}}.
#' @param constraints A one-sided formula \code{\link{formula}} giving
#'   the sample space constraints. See \code{\link[ergm]{ergm}} for
#'   details and examples.
#' @param popsize The size \eqn{|N|} of the finite population network from
#' which the egocentric sample was taken; only affects the shift in the
#' coefficients of the terms modeling the overall propensity to have ties.
#' Setting it to 1 (the default) essentially uses the \eqn{-\log |N'|} offset
#' on the edges term.
#' @param offset.coef A vector of coefficients for the offset terms.
#' @param na.action How to handle missing actor attributes in egos or alters,
#' when the terms need them.
#' @param \dots Additional arguments passed to \code{\link[ergm]{ergm}}.
#' @param control A \code{\link{control.ergm.ego}} control list.
#' @param do.fit Whether to actually call \code{\link[ergm]{ergm}}
#' @return An object of class \code{ergm.ego} inheriting from
#' \code{\link[ergm]{ergm}}, with the following additional or overridden
#' elements: \item{"v"}{Variance-covariance matrix of the estimate of the
#' sufficient statistics} \item{"m"}{Estimate of the sufficient
#' statistics} \item{"egor"}{The egor object passed}
#' \item{"popsize"}{Population network size and pseudopopulation size
#' used, respectively}\item{, }{Population network size and pseudopopulation
#' size used, respectively}\item{"ppopsize"}{Population network size and
#' pseudopopulation size used, respectively} \item{"coef"}{The
#' coefficients, along with the network size adjustment \code{netsize.adj}
#' coefficient.} \item{"covar"}{Pseudo-MLE estimate of the
#' variance-covariance matrix of the parameter estimates under repeated
#' egocentric sampling} \item{"ergm.covar"}{ The variance-covariance matrix of
#' parameter estimates under the ERGM superpopulation process (without
#' incorporating sampling).  }
#' \item{"DtDe"}{Estimated Jacobian of the expectation of the sufficient
#' statistics with respect to the model parameters}
#' @author Pavel N. Krivitsky
#' @references
#' 
#' Pavel N. Krivitsky and Martina Morris. Inference for Social Network Models
#' from Egocentrically-Sampled Data, with Application to Understanding
#' Persistent Racial Disparities in HIV Prevalence in the US. Thechnical
#' Report. National Institute for Applied Statistics Research Australia,
#' University of Wollongong, 2015(05-15).
#' \url{http://niasra.uow.edu.au/publications/UOW190187.html}
#' @keywords models
#' @examples
#' 
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#' 
#' head(fmh.ego)
#' 
#' egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'                          +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade")+transitiveties, 
#'                           popsize=network.size(faux.mesa.high))
#' 
#' # Run convergence diagnostics
#' mcmc.diagnostics(egofit)
#' 
#' # Estimates and standard errors
#' summary(egofit)
#' 
#' # Note that we recover the ergm() parameters
#' \dontrun{
#' coef(ergm(faux.mesa.high~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'                          +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade")+transitiveties,
#'           eval.loglik=FALSE))
#' }
#' 
#' rbind(c(0, -2.9868, -0.154, 0.1162, -0.1869, 0.1431, -1.1494, -0.8786,
#'            -2.2698, -0.6101, 0.7227, -0.1537, 0.5557, -1.0783, 1.5875),
#'       coef(egofit))
#'
#' \dontshow{
#'
#' stopifnot(isTRUE(all.equal(c(0, -2.9868, -0.154, 0.1162, -0.1869,
#' 0.1431, -1.1494, -0.8786,-2.2698, -0.6101, 0.7227, -0.1537, 0.5557,
#' -1.0783, 1.5875), coef(egofit), tolerance=.05,
#' check.attributes=FALSE)))
#'
#' }
#' @import ergm stats
#' @importFrom utils modifyList
#' @export
ergm.ego <- function(formula, popsize=1, offset.coef=NULL, constraints=~.,..., control=control.ergm.ego(), na.action=na.fail, do.fit=TRUE){
  statnet.common::check.control.class()
  
  stats.est <- control$stats.est
  stats.wt <- control$stats.wt
  egor <- get(as.character(formula[[2]]), envir=environment(formula))

  sampsize <- nrow(egor)
  ppopsize <-
    if(is.data.frame(control$ppopsize)) nrow(control$ppopsize)
    else if(is.numeric(control$ppopsize)) control$ppopsize
    else switch(control$ppopsize,
                auto = if(missing(popsize) || popsize==1) sampsize*control$ppopsize.mul else popsize*control$ppopsize.mul,  
                samp = sampsize*control$ppopsize.mul,
                pop = popsize*control$ppopsize.mul)
  
  if(ppopsize < sampsize && !is.data.frame(control$ppopsize)) warning("Using a smaller pseudopopulation size than sample size usually does not make sense.")
  else if(ppopsize == sampsize && var(weights(egor))>sqrt(.Machine$double.eps))
    warning("Using pseudopoulation size equal to sample size under weighted sampling: results may be highly biased. Recommend increasing popsize.mul control parameter.")
  
  message("Constructing pseudopopulation network.")
  popnw <-
    if(is.data.frame(control$ppopsize)){ # If pseudopoluation composition is given in popsize, use that.
      pegos <- control$ppopsize
      pdata <- egor(alters.df=rep(list(data.frame()), nrow(pegos)), egos.df=pegos)
      as.network(pdata, ppopsize)
    }else{
      as.network(egor, ppopsize, scaling=control$ppop.wt)
    }

  if(network.size(popnw)!=ppopsize){
    message("Note: Constructed network has size ", network.size(popnw), ", different from requested ", ppopsize,". Estimation should not be meaningfully affected.")
    ppopsize <- network.size(popnw)
  }

  w <- switch(stats.wt,
              data=weights(egor),
              ppop=tabulate(popnw %v% "ego.ind", nbins=nrow(egor))
              )
  
  # Get the sample h values.
  stats <- try(summary(remove.offset.formula(formula), individual=TRUE))
  ord <- attr(stats, "order")  
  adj.update <- call("~",as.name("."),call("+", call("offset", call("netsize.adj", +(1%in%ord), -(2%in%ord), -1/3*(3%in%ord))), as.name(".")))

  if(!inherits(stats,"try-error")){
    # h is just a matrix, so this will do the sensible thing.
    tmp <- na.action(cbind(w,stats))
    w <- tmp[,1]
    stats <- tmp[,-1, drop=FALSE] * ppopsize
    n <- length(w)
    
    wmean <- function(w,s){
      w <- w/sum(w)
      colSums(s*w)
    }
    m <- wmean(w,stats)
    
    if(stats.est=="bootstrap"){
      m.b <- t(rbind(replicate(control$boot.R,{
                           i <- sample.int(length(w),replace=TRUE)
                           wmean(w[i],stats[i,,drop=FALSE])
                         })))
      m <- m - (colMeans(m.b)-m)
      
    }else if(stats.est=="jackknife"){
      m.j <- t(rbind(sapply(seq_len(n), function(i){
                        wmean(w[-i],stats[-i,,drop=FALSE])
                      })))
      m <- n*m - (n-1)*colMeans(m.j)
    }
    
    # TODO: Include finite-population correction here:
    v <- switch(stats.est,
                bootstrap = cov(m.b),
                jackknife = (n-1)/n*crossprod(sweep(m.j,2,colMeans(m.j))),
                naive = {w <- w/sum(w); crossprod(sweep(stats, 2, m, "-")*sqrt(w))/(1-sum(w^2))*sum(w^2)},
                asymptotic = .asymptotic.var(stats, w)/length(w)
                )
  }else{
    if(stats.est %in% c("naive","asymptotic"))
      stop("Non-scaling statistic detected: use bootstrap or jackknife variance estimator.")
    if(do.fit && popsize!=ppopsize)
      warning("Non-scaling statistic detected when trying to fit a model: network-size invariant parametrization probably does not exist so pseudopopulation size should equal the population size.")

    n <- nrow(egor)
    m <- summary(remove.offset.formula(formula), basis=egor, individual=FALSE, scaleto=ppopsize)
    ord <- attr(m, "order")  
    adj.update <- call("~",as.name("."),call("+", call("offset", call("netsize.adj", +(1%in%ord), -(2%in%ord), -1/3*(3%in%ord))), as.name(".")))
      
    if(stats.est=="bootstrap"){
      m.b <- t(replicate(control$boot.R,{
                           i <- sample.int(length(w),replace=TRUE)
                           e <- egor[i,]
                           summary(remove.offset.formula(formula), basis=e, individual=FALSE, scaleto=ppopsize)
                         }))
      m <- m - (colMeans(m.b)-m)
      
    }else if(stats.est=="jackknife"){
      m.j <- t(sapply(seq_len(n), function(i){
                        e <- egor[-i,]
                        summary(remove.offset.formula(formula), basis=e, individual=FALSE, scaleto=ppopsize)
                      }))
      m <- n*m - (n-1)*colMeans(m.j)
    }
    
    # TODO: Include finite-population correction here:
    v <- switch(stats.est,
                bootstrap = cov(m.b),
                jackknife = (n-1)/n*crossprod(sweep(m.j,2,colMeans(m.j)))
                )
  }

  ergm.formula <- ergm.update.formula(formula,popnw~.,from.new="popnw")
  ergm.formula <- ergm.update.formula(ergm.formula,adj.update)

  ergm.offset.coef <- c(-log(ppopsize/popsize),offset.coef)

  # If nominations were limited, represent the cap a degree bound.
  constraints <- if(!control$ignore.max.alters && alter.design(egor)$max<Inf){
    newterm <- as.formula(substitute(~bd(maxout=.max), eval(list(.max=alter.design(egor,"max")))))
    if(constraints==~.)
      newterm
    else
      append.rhs.formula(constraints, term.list.formula(newterm))
  }else constraints
  
  out <- list(v=v, m=m, formula=formula, ergm.formula=ergm.formula, offset.coef=offset.coef, ergm.offset.coef=ergm.offset.coef, egor=egor, ppopsize=ppopsize, popsize=popsize, constraints=constraints, netsize.adj=adj.update)

  if(do.fit){

    ergm.fit <- ergm(ergm.formula, target.stats=m, offset.coef=ergm.offset.coef, constraints=constraints, ..., eval.loglik=FALSE,control=control$ergm.control)

    ## Workaround to keep mcmc.diagnostics from failing. Should be removed after fix is released.
    if(inherits(ergm.fit$sample,"mcmc.list")){
      for(thread in 1:coda::nchain(ergm.fit$sample))
        ergm.fit$sample[[thread]][,1] <- 0
    } else ergm.fit$sample[,1] <- 0
    ergm.fit$drop[1] <- 0

    coef <- coef(ergm.fit)

    oi <- offset.info.formula(ergm.formula)

    DtDe <- -ergm.fit$hessian[!oi$theta,!oi$theta,drop=FALSE]

    vcov <- matrix(NA, length(coef), length(coef))
  
    vcov[!oi$theta,!oi$theta] <- solve(DtDe)%*%v%*%solve(DtDe)
    
    rownames(vcov) <- colnames(vcov) <- names(coef)

    out <- c(out, list(covar=vcov, ergm.covar=ergm.fit$covar, DtDe=DtDe))
    out <- modifyList(ergm.fit, out)
  }
  class(out) <- c("ergm.ego","ergm")
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

