#  File R/control.simulate.ergm.ego.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################


#' Control parameters for \code{\link{simulate.ergm.ego}}.
#' 
#' Constructs and checks the list of control parameters for simulation by
#' \code{\link{simulate.ergm.ego}}.
#' 
#' 
#' @param ppop.wt Because each ego must be represented in the pseuodopopulation
#' network an integral number of times, if the sample is weighted (or the
#' target \eqn{|N'|} calculated from \code{ppopsize} and \code{ppopsize.mul} is
#' not a multiple of the sample size), it may not be possible, for a finite
#' \eqn{|N'|} to represent each ego exactly according to its relative weight,
#' and \code{ppop.wt} controls how the fractional egos are allocated:
#' \describe{ \item{"round"}{(default) Rather than treating \code{ppopsize} as
#' a hard setting, calculate \eqn{|N'| w_i / w_\cdot} for each ego \eqn{i} and
#' round it to the nearest integer. Then, the \eqn{|N'|} actually used will be
#' the sum of these rounded freqencies.}
#' \item{"sample"}{Resample in proportion to \eqn{w_i}.} }
#' 
#' @param SAN.control A list of control parameters for \code{\link[ergm]{san}}
#' constructed by \code{\link[ergm]{control.ergm}}, called to construct a
#' pseudopopulation network consistent with the data.
#' @param simulate.control A list of control parameters for
#' \code{\link[ergm]{simulate.formula}} constructed by
#' \code{\link[ergm]{control.simulate}}, called to simulate from the model fit.
#' @param \dots Not used at this time.
#' @return A list with arguments as components.
#' @author Pavel N. Krivitsky
#' @seealso control.simulate, control.san
#' @keywords models
#' @export
control.simulate.ergm.ego <- function(
  ppop.wt = c("round","sample"),
  SAN.control = control.san(),
  simulate.control = control.simulate(),
  ...){
  match.arg.pars <- c("ppop.wt")

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))
  if(length(list(...))) stop("Unrecognized control parameter: ",arg,".")

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class("control.simulate.ergm.ego")
}
