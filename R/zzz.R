#  File R/zzz.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################
#' @import statnet.common
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





# To be used inside testthat tests, i.e. within test code of test_that(). Checks
# environment variable SKIP_ERGM_EGO_LONG_TEST. If it is set to 1 the test is
# skipped.
long_test <- function() {
  # check envar
  v <- Sys.getenv("SKIP_ERGM_EGO_LONG_TEST")
  do_skip <- isTRUE( as.numeric(v) == 1 )
  if(do_skip) {
    testthat::skip("Not running this lengthy test because SKIP_ERGM_EGO_LONG_TEST=1")
  }
}
