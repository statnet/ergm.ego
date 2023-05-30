#  File R/fmhfit.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2023 Statnet Commons
################################################################################
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
#'   + absdiff("Grade") + gwesp(0, fix=TRUE), 
#'   popsize = network.size(faux.mesa.high),
#'   control = control.ergm.ego(
#'     ergm = control.ergm(parallel=2)
#'   )
#' )
#' }
NULL
