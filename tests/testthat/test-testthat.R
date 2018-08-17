#  File tests/testthat/test-testthat.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################
context("Dummy testing of testthat infrastructure")

test_that("there is ergm.ego()", {
  expect_true(exists("ergm.ego", where="package:ergm.ego"))
})