#  File R/ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################


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
#'
#' @param popsize The size \eqn{|N|} of the finite population network
#'   from which the egocentric sample was taken; only affects the
#'   shift in the coefficients of the terms modeling the overall
#'   propensity to have ties.  Setting it to 1 (the default)
#'   essentially uses the \eqn{-\log |N'|} offset on the edges
#'   term. Passing 0 disables network size adjustment and uses the
#'   egocentric sample size; passing [`I(N)`][I] uses the specified
#'   size `N` (though can be overridden by the `ppop`
#'   [control.ergm.ego()] option) and disables network size
#'   adjustment.
#'
#' @param offset.coef A vector of coefficients for the offset terms.
#' @param na.action How to handle missing actor attributes in egos or alters,
#' when the terms need them for  models that scale.
#' @param na.rm How to handle missing actor attributes in egos or alters,
#' when the terms need them for models that do not scale.
#' @param \dots Additional arguments passed to \code{\link[ergm]{ergm}}.
#' @param control A \code{\link{control.ergm.ego}} control list.
#' @param do.fit Whether to actually call \code{\link[ergm]{ergm}}
#' @return An object of class \code{ergm.ego} inheriting from
#' \code{\link[ergm]{ergm}}, with the following additional or overridden
#' elements:
#' \item{"v"}{Variance-covariance matrix of the estimate of the
#' sufficient statistics}
#' \item{"m"}{Estimate of the sufficient
#' statistics}
#' \item{"egor"}{The [`egor`] object passed}
#' \item{"popsize"}{Population network size used}
#' \item{"ppopsize"}{Pseudopopulation size used, see \code{\link{control.ergm.ego}}}
#' \item{"coef"}{The
#' coefficients, along with the network size adjustment \code{netsize.adj}
#' coefficient.}
#' \item{"covar"}{Pseudo-MLE estimate of the
#' variance-covariance matrix of the parameter estimates under repeated
#' egocentric sampling}
#' \item{"ergm.covar"}{ The variance-covariance matrix of
#' parameter estimates under the ERGM superpopulation process (without
#' incorporating sampling).  }
#' \item{"DtDe"}{Estimated Jacobian of the expectation of the sufficient
#' statistics with respect to the model parameters}
#' @author Pavel N. Krivitsky
#' @references
#' 
#' 
#' Pavel N. Krivitsky and Martina Morris (2017). "Inference for social network models from egocentrically sampled data, with application to understanding persistent racial disparities in HIV prevalence in the US." *Annals of Applied Statistics*, 11(1): 427–455. \doi{10.1214/16-AOAS1010}
#'
#' Pavel N. Krivitsky, Martina Morris, and Michał Bojanowski (2019). "Inference for Exponential-Family Random Graph Models from Egocentrically-Sampled Data with Alter–Alter Relations." NIASRA Working Paper 08-19. \url{https://www.uow.edu.au/niasra/publications/}
#'
#' Pavel N. Krivitsky, Michał Bojanowski, and Martina Morris (2020). "Impact of survey design on estimation of exponential-family random graph models from egocentrically-sampled data." *Social Networks*, to appear. \doi{10.1016/j.socnet.2020.10.001}
#' 
#' Pavel N. Krivitsky, Mark S. Handcock, and Martina Morris (2011). "Adjusting for
#' Network Size and Composition Effects in Exponential-Family Random Graph
#' Models." \emph{Statistical Methodology}, 8(4): 319–339. \doi{10.1016/j.stamet.2011.01.005}
#'
#' @keywords models
#' @seealso \code{\link[ergm]{ergm}()}
#' @examples
#' \donttest{
#' data(faux.mesa.high)
#' fmh.ego <- as.egor(faux.mesa.high)
#' 
#' head(fmh.ego)
#' 
#' egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
#'                          +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade")+gwesp(0,fix=TRUE), 
#'                           popsize=network.size(faux.mesa.high))
#' 
#' # Run convergence diagnostics
#' mcmc.diagnostics(egofit)
#' 
#' # Estimates and standard errors
#' summary(egofit)
#' }
#' 
#' @import ergm
#' @importFrom utils modifyList
#' @export
ergm.ego <- function(formula, popsize=1, offset.coef=NULL, constraints=~.,..., control=control.ergm.ego(), na.action=na.fail, na.rm=FALSE, do.fit=TRUE){
  statnet.common::check.control.class("ergm.ego","ergm.ego")

  ergm.ego_call <- match.call(ergm)

  stats.est <- control$stats.est
  stats.wt <- control$stats.wt
  egor <- eval_lhs.formula(formula)

  sampsize <- nrow(egor$ego)
  ppopsize <-
    if(is.network(control$ppopsize)) network.size(control$ppopsize)
    else if(is.data.frame(control$ppopsize)) nrow(control$ppopsize)
    else if(is.numeric(control$ppopsize)) control$ppopsize
    else switch(control$ppopsize,
                auto = if(is(popsize, "AsIs")) popsize
                else if(missing(popsize) || popsize %in% c(0,1)) sampsize*control$ppopsize.mul else popsize*control$ppopsize.mul,
                samp = sampsize*control$ppopsize.mul,
                pop = popsize*control$ppopsize.mul)

  nsa <- !is(popsize, "AsIs") && popsize > 0
  
  if(ppopsize < sampsize && !is.data.frame(control$ppopsize)) warning("Using a smaller pseudopopulation size than sample size usually does not make sense.")
  else if(ppopsize == sampsize && var(weights(egor))>sqrt(.Machine$double.eps))
    warning("Using pseudopoulation size equal to sample size under weighted sampling: results may be highly biased. Recommend increasing popsize.mul control parameter.")
  
  message("Constructing pseudopopulation network.")
  popnw <-
    if(is.network(control$ppopsize)){  # If pseudopopulation network is given in popsize, use that.
      control$ppopsize
    }else if(is.data.frame(control$ppopsize)){ # If pseudopoluation composition is given in popsize, use that.
      template_network(control$ppopsize, ppopsize)
    }else{
      template_network(egor, ppopsize, scaling=control$ppop.wt)
    }

  if(network.size(popnw)!=ppopsize){
    if(nsa) message("Note: Constructed network has size ", network.size(popnw), ", different from requested ", ppopsize,". Estimation should not be meaningfully affected.")
    else warning("Constructed network has size ", network.size(popnw), ", different from requested ", ppopsize,", with network size adjustment disabled; estimation may be affected.", immediate.=TRUE)
    ppopsize <- network.size(popnw)
  }

  w <- switch(stats.wt,
              data=weights(egor),
              ppop=tabulate(popnw %v% "ego.ind", nbins=nrow(egor))
              )

  deoffset <- function(f)
    filter_rhs.formula(f, function(x)
    (if(is.call(x)) x[[1]] else x)!="offset")
  
  # Try to get the sample h values.
  if(stats.est != "survey") stats <- try(summary(deoffset(formula), individual=TRUE))

  if(stats.est!="survey" && !inherits(stats,"try-error")){
    ord <- attr(stats, "order")  
    adj.update <- call("~",as.name("."),call("+", call("offset", call("netsize.adj", edges = +(1%in%ord), mutual = -(2%in%ord), transitiveties = -1/3*(3%in%ord))), as.name(".")))
    
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
    m <- summary(deoffset(formula), basis=egor, na.rm=na.rm, individual=FALSE, scaleto=ppopsize)
    ord <- attr(m, "order")  
    adj.update <- call("~",as.name("."),call("+", call("offset", call("netsize.adj", edges = +(1%in%ord), mutual = -(2%in%ord), transitiveties = -1/3*(3%in%ord))), as.name(".")))

    if(0 %in% ord && stats.est %in% c("survey", "naive","asymptotic"))
      stop("Non-scaling statistic detected: use bootstrap or jackknife variance estimator.")
    if(0 %in% ord && do.fit && popsize!=ppopsize)
      warning("Non-scaling statistic detected when trying to fit a model: network-size invariant parametrization probably does not exist so pseudopopulation size should equal the population size.")

    n <- nrow(egor)

    if(stats.est=="survey"){
      v <- vcov(m)
      m <- setNames(as.vector(m), names(m))
    }else if(stats.est=="bootstrap"){
      m.b <- t(replicate(control$boot.R,{
                           i <- sample.int(length(w),replace=TRUE)
                           e <- egor$ego[i,]
                           summary(deoffset(formula), basis=e, individual=FALSE, scaleto=ppopsize)
                         }))
      m <- m - (colMeans(m.b)-m)
      
    }else if(stats.est=="jackknife"){
      m.j <- t(sapply(seq_len(n), function(i){
                        e <- egor$ego[-i,]
                        summary(deoffset(formula), basis=e, individual=FALSE, scaleto=ppopsize)
                      }))
      m <- n*m - (n-1)*colMeans(m.j)
    }
    
    # TODO: Include finite-population correction for non-survey methods here:
    v <- switch(stats.est,
                survey = v,
                bootstrap = cov(m.b),
                jackknife = (n-1)/n*crossprod(sweep(m.j,2,colMeans(m.j)))
                )
  }
  
  ergm.formula <- nonsimp_update.formula(formula,popnw~.,from.new="popnw")

  ergm.names <- param_names(ergm_model(ergm.formula), canonical=TRUE, offset=FALSE)
  if(!setequal(ergm.names,names(m))){
    ergm.not.ts <- setdiff(ergm.names, names(m))
    ts.not.ergm <- setdiff(names(m), ergm.names)
    errstr <-
      if(length(ergm.not.ts)) paste("statistics", paste.and(sQuote(ergm.not.ts)), "required by the ERGM could not be estimated from data.")
      else paste("statistics", paste.and(sQuote(ts.not.ergm)), "estimated from data are extraneous to the ERGM.")
    stop("There appears to be a mismatch between estimated statistic and the sufficient statistic of the ERGM: ", errstr, " A common cause of this is that egos and alters do not have a consistent set of levels for one or more factors.")
  }
  
  if(nsa) ergm.formula <- nonsimp_update.formula(ergm.formula,adj.update)
  ergm.offset.coef <- c(if(nsa) -log(ppopsize/popsize),offset.coef)

  # If nominations were limited, represent the cap a degree bound.
  constraints <- if(!control$ignore.max.alters && alter_design(egor)$max<Inf){
    newterm <- as.formula(substitute(~bd(maxout=.max), eval(list(.max=alter_design(egor,"max")))))
    if(constraints==~.)
      newterm
    else
      append_rhs.formula(constraints, list_rhs.formula(newterm))
  }else constraints
  
  out <- list(v=v, m=m, formula=formula, ergm.formula=ergm.formula, offset.coef=offset.coef, ergm.offset.coef=ergm.offset.coef, egor=egor, ppopsize=ppopsize, popsize=popsize, constraints=constraints, netsize.adj=if(nsa) adj.update, call=ergm.ego_call)

  if(do.fit){
    ergm.fit <- ergm(ergm.formula, target.stats=m, offset.coef=ergm.offset.coef, constraints=constraints, ..., eval.loglik=FALSE,control=control$ergm)
    if(is.curved(ergm.fit)) warning("Theory of egocentric inference and particularly of variance calculation for curved ERGMs is not well understood; standard errors might not be reliable.")

    coef <- coef(ergm.fit)

    oi <- ergm.fit$etamap$offsettheta
    dropped <- ergm.fit$drop != 0
    
    DtDe <- -ergm.fit$hessian[!oi,!oi,drop=FALSE]

    novar <- diag(DtDe)<sqrt(.Machine$double.eps)

    if(any(novar)) warning("Unable to estimate standard errors for parameters ",paste.and(names(novar)[novar],oq="'",cq="'"),". Try increasing MCMC sample size and interval.")
    oi[!oi][novar] <- TRUE
    
    vcov <- matrix(NA, length(coef), length(coef))

    iDtDe <- solve(DtDe[!novar,!novar,drop=FALSE])

    # Augment the statistics covariance matrix with 0 rows and columns
    # for offsets, then pre- and post-multiply it by the derivative of
    # theta with to eta.
    vdropped <- dropped[!oi | dropped]
    vaug <- matrix(0, length(ergm.fit$etamap$offsetmap), length(ergm.fit$etamap$offsetmap))
    vaug[!ergm.fit$etamap$offsetmap, !ergm.fit$etamap$offsetmap] <- v[!vdropped,!vdropped]
    vaug <- ergm.etagradmult(coef, t(ergm.etagradmult(coef, vaug, ergm.fit$etamap)), ergm.fit$etamap)
    
    vcov[!oi,!oi] <- iDtDe%*%vaug[!oi,!oi,drop=FALSE]%*%iDtDe

    rownames(vcov) <- colnames(vcov) <- names(coef)

    out <- c(out, list(covar=vcov, ergm.covar=ergm.fit$covar, DtDe=DtDe, ergm.call=ergm.fit$call))
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

