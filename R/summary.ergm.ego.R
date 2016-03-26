## A thin wrapper around summary.ergm to get rid of a spurious error message.
.summary.ergm.ego <- function (object, ..., 
                          digits = max(3, getOption("digits") - 3),
                          correlation=FALSE, covariance=FALSE,
                          total.variation=TRUE){
  class(object) <- "ergm"
  summ <- NextMethod("summary")
  class(summ) <- c("summary.ergm.ego", "summary.ergm")
  summ
}

.print.summary.ergm.ego <- function (x, ...){
  print.summary.ergm(x, ..., print.deviances=FALSE)
}

