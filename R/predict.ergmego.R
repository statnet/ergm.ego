#  File R/predict.ergmego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################
#' ERGM-based predicted tie probabilities for the pseudo-population network
#' 
#' @param object model fit as returned by [ergm.ego()]
#' @param ... other arguments passed to/from other methods
#' 
#' @return See [ergm::predict.ergm()]
#' 
#' @export

predict.ergm.ego <- function(object, ...) {
  # Extract network
  net <- object$network
  # Update formula with pseudo-population network
  frm <- statnet.common::nonsimp_update.formula(object$ergm.formula, net ~ .)
  assign("net", object$network, envir=environment(frm))
  # thetas
  th <- ergm.eta(coef(object), object$etamap)
  predict(frm, theta=th, ...)
}
