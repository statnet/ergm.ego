#  File tests/testthat/test-table_ppop.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################
test_that("estimation and simulation work", {
  data(faux.mesa.high)
  fmh.ego <- as.egor(faux.mesa.high)
  ppop <- fmh.ego$ego[rep(seq_len(nrow(fmh.ego$ego)), each=2),]

  set.seed(0)
  (egofit_egopop <- ergm.ego(
     fmh.ego ~ edges + degree(0:3) + nodefactor("Race") + nodematch("Race")
     + nodefactor("Sex") + nodematch("Sex") + absdiff("Grade"),
     popsize = network.size(faux.mesa.high), estimate = "MPLE",
     control = control.ergm.ego(
       ppopsize=ppop
     )
   )) |> expect_error(NA) |> expect_warning(NA)

  set.seed(0)
  (egofit_scl <- ergm.ego(
     fmh.ego ~ edges + degree(0:3) + nodefactor("Race") + nodematch("Race")
     + nodefactor("Sex") + nodematch("Sex") + absdiff("Grade"),
     popsize = network.size(faux.mesa.high), estimate = "MPLE",
     control = control.ergm.ego(
       ppopsize=2*network.size(faux.mesa.high)
     )
   )) |> expect_error(NA) |> expect_warning(NA)

  expect_equal(coef(egofit_scl), coef(egofit_egopop))

  ppop <- ppop[sample.int(nrow(ppop), nrow(ppop)*1.5, replace=TRUE),]
  
  (egosim <- simulate(egofit_scl, popsize=ppop)) |> expect_silent()
})
