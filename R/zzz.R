#  File R/zzz.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################
#' @importFrom statnet.common statnetStartupMessage ERRVL NVL NVL3 all_identical despace eval_lhs.formula filter_rhs.formula list_rhs.formula nonsimp_update.formula paste.and set.control.class
.onAttach <- function(lib, pkg){
  sm <- statnetStartupMessage("ergm.ego", c("statnet"), TRUE)
  if(!is.null(sm)){
    packageStartupMessage(sm)
  }
  
}

.onLoad <- function(lib, pkg){
  # . is used in purrr and elsewhere.
  utils::globalVariables(".")
}
