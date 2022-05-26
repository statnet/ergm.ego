#  File R/control.ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################


#' Control parameters for \code{\link{ergm.ego}}.
#' 
#' Constructs and checks the list of control parameters for estimation by
#' \code{\link{ergm.ego}}.
#' 
#' 
#' @param ppopsize,ppopsize.mul Parameters to determine the size
#'   \eqn{|N'|} of the pseudopopulation network. \code{ppopsize} can be
#'   \describe{
#' 
#' \item{"auto"}{If the \code{popsize} (\eqn{|N|}) argument is
#' specified and is different from 1, as if \code{"pop"}; otherwise,
#' as \code{"samp"}.}
#'
#' \item{"samp"}{set \eqn{|N'|} based on the sample size:
#' \eqn{|N'|=|S| \times \code{popsize.mul}}}
#' 
#' \item{"pop"}{set \eqn{|N'|} based on the population size:
#' \eqn{|N'|=|N| \times \code{popsize.mul}}}
#'
#' \item{a number}{set \eqn{|N'|} directly (\code{popsize.mul}
#' ignored)}
#'
#' \item{a [`network`] object}{use the specified network as the
#' pseudo-population network directly; use at your own risk}
#'
#' \item{a data frame}{use the specified data frame as the
#' pseudo-population; use at your own risk}}
#' 
#' The default is to use the same pseudopopulation size as the sample size,
#' but, particularly if there are sampling weights in the data, it should be
#' bigger.
#' 
#' Note that depending on \code{ppop.wt}, this may only be an approximate
#' target specification, with the actual constructed pseudopopulation network
#' being slightly bigger or smaller.
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
#' @param stats.wt Weight assigned to each ego's contribution to the ERGM's
#' sufficient statistic: \describe{ \item{"data"}{(default) Use weights
#' \eqn{|N'| w_i / w_\cdot} for each ego \eqn{i} as in the data.} \item{"ppop"}{Use weights ultimately used in the
#' pseudopopulation network.} }
#' @param stats.est,boot.R Method to be used to estimate the ERGM's sufficient
#' statistics and their variance: \describe{
#' \item{"survey"}{Variance estimator returned by [survey::svymean()], appropriate to the design of the dataset.}
#' \item{"asymptotic"}{Delta method, as derived by Krivitsky and
#' Morris (2017), assuming the ego weights are sampled alongside the
#' egos.}\item{ (default)}{Delta method, as derived by Krivitsky and Morris
#' (2017), assuming the ego weights are sampled alongside the egos.}
#' \item{"bootstrap"}{Nonparametric bootstrap with bias correction,
#' resampling egos, using \code{R} replications.}
#' \item{"jackknife"}{Jackknife with bias correction.}
#' \item{"naive"}{"Naive" estimator, assuming that weights are
#' fixed.} }
#' 
#' @param ignore.max.alters if `TRUE`, ignores any constraints on the
#'   number of nominations. Used to be `FALSE`, now `TRUE` in light of
#'   the findings of Krivitsky et. al (2020).
#' 
#' @param ergm Control parameters for the \code{\link[ergm]{ergm}()} call
#' to fit the model, constructed by \code{\link[ergm]{control.ergm}()}.
#' @param \dots Not used at this time.
#' @return A list with arguments as components.
#' @author Pavel N. Krivitsky
#' @seealso \code{\link[ergm]{control.ergm}()}
#' @references
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
#' @export
control.ergm.ego <- function(
  ppopsize = c("auto", "samp", "pop"),
  ppopsize.mul = 1,
  ppop.wt = c("round","sample"),                                      
  stats.wt = c("data","ppop"),
  stats.est = c("survey","asymptotic", "bootstrap", "jackknife", "naive"),
  boot.R = 10000,
  ignore.max.alters = TRUE,
  ergm = control.ergm(),
  ...){
  old.controls <- list(ergm.control="ergm")

  match.arg.pars <- c("stats.est", "ppop.wt", "stats.wt", if(is.character(ppopsize)) "ppopsize")

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))

  for(arg in names(list(...))){
    if(!is.null(old.controls[[arg]])){
      warning("Passing ",arg," to control.ergm.ego(...) is deprecated and may be removed in a future version. Specify it as control.ergm.ego(",old.controls[[arg]],"=...) instead.")
      control[old.controls[[arg]]]<-list(list(...)[[arg]])
    }else{
      stop("Unrecognized control parameter: ",arg,".")
    }
  }

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class("control.ergm.ego")
}
