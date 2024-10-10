#  File tests/testthat/test-drop.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################

nw <- network.initialize(20,directed=FALSE)

# No homophilous ties; heterophilous density of 1/2.
nw %v% "a" <- rep(1:2, each=10)
nw[1:10,11:20] <- 0:1

test_that("dropped ergm terms", {
  out.coef <- coef(ergm.ego(as.egor(nw)~edges+nodematch("a")))
  out.coef <- c(sum(out.coef[1:2]),out.coef[3])
  expect_equal(coef(ergm(nw~edges+nodematch("a"))), out.coef, tolerance=0.1, ignore_attr=TRUE)
})
