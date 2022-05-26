#  File R/zzz.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################
#' @import statnet.common
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.ego", c("statnet"), TRUE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
}

.onLoad <- function(libname, pkgname){
  # . is used as a placeholder by stantet.common::NVL3().
  utils::globalVariables(".")

  eval(COLLATE_ALL_MY_CONTROLS_EXPR)
}


# TODO: Figure out some automatic way to keep this in sync with statnet.common.
#' @name snctrl
#'
#' @title Statnet Control
#'
#' @description A utility to facilitate argument completion of control lists, reexported from `statnet.common`.
#'
#' @section Currently recognised control parameters:
#' This list is updated as packages are loaded and unloaded.
#'
#' \Sexpr[results=rd,stage=render]{statnet.common::snctrl_names()}
#'
#' @seealso [statnet.common::snctrl()]
#' @docType import
NULL
#' @export
snctrl <- statnet.common::snctrl

eval(UPDATE_MY_SCTRL_EXPR)
