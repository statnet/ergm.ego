#  File R/predict.ergmego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2025 Statnet Commons
################################################################################
#' ERGM-based predicted tie probabilities for the pseudo-population network
#'
#' @param object model fit as returned by [ergm.ego()]
#' @param ... other arguments passed to/from other methods
#'
#' @return See [ergm::predict.ergm()]
#'
#' @importFrom stats predict
#' @export
predict.ergm.ego <- function(object, ...) {
  predict(object$ergm.formula, eta = ergm.eta(coef(object), object$etamap),
          basis = object$network, ...)
}
