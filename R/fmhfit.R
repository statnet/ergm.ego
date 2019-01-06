#' Fitted ergm.ego model object
#' 
#' This is an object with a fitted model to `faux.mesa.high` data using the code
#' shown below in the Examples section.
#'
#' @docType data
#' @name fmhfit
#'
#' @format An object of class `ergm.ego`.
#' 
#' @examples
#' \dontrun{
#' data(faux.mesa.high)
#' fmh.ego <- egor::as.egor(faux.mesa.high)
#' fmhfit <- ergm.ego(
#'   fmh.ego ~ edges + degree(0:3) + 
#'     nodefactor("Race") + nodematch("Race")
#'   + nodefactor("Sex")+nodematch("Sex")
#'   + absdiff("Grade") + transitiveties, 
#'   popsize = network.size(faux.mesa.high),
#'   control = control.ergm.ego(
#'     ergm.control = control.ergm(parallel=2)
#'   )
#' )
#' }
NULL
