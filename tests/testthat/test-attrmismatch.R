#  File tests/testthat/test-attrmismatch.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2025 Statnet Commons
################################################################################

egos <- data.frame(egoIDcol = 1:12, x = rep(1:3, 4))
alters <- data.frame(egoIDcol = sample(1:12, 24, TRUE), x = rep(1:3, 8))

e <- egor(egos=egos,
          alters=alters,
          ID.vars=list(ego="egoIDcol"))

test_that("no attribute mismatch error", {
  expect_error(ergm.ego(e ~ nodefactor("x")), NA)
})

test_that("attribute mismatch error", {
  e$alter$x[1] <- 4
  
  expect_warning(ergm.ego(e ~ nodefactor("x")), "In ego statistic 'nodefactor' from package 'ergm.ego', the set of levels of categorical attribute 'x' among alters (.*) is not a subset of its levels among the egos (.*); it is likely that levels and/or the pseudo population network will need to be specified explicitly.")
})
