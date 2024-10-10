#  File R/summary.ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################
## A thin wrapper around summary.ergm to get rid of a spurious error message.
#' @method summary ergm.ego
#' @export
summary.ergm.ego <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE){
  class(object) <- "ergm"
  summ <- NextMethod("summary")
  class(summ) <- c("summary.ergm.ego", "summary.ergm")
  summ
}

#' @method print summary.ergm.ego
#' @export
print.summary.ergm.ego <- function (x, ...){
  class(x) <- "summary.ergm"
  NextMethod("print", object=x, ..., print.deviances=FALSE)
}

#' @method logLik ergm.ego
#' @export
#' @noRd
logLik.ergm.ego <- function(object, ...){
  stop("Log-likelihood is not meaningful for moments-based inference used by ", sQuote("ergm.ego()"), ".")
}
