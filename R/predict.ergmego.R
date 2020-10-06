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
  th <- ergm.eta(object$coef, object$etamap)
  predict(frm, theta=th, ...)
}
