#  File tests/testthat/test-coef_recovery.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2021 Statnet Commons
################################################################################
test_that("complete ERGM and ergm.ego() give similar coef estimates",{
  ergm.ego:::long_test()

  data(faux.mesa.high)
  fmh.ego <- egor::as.egor(faux.mesa.high)
  
  fit <- ergm(
    faux.mesa.high ~ edges + degree(0:3) + nodefactor("Race") + nodematch("Race") + 
      nodefactor("Sex") + nodematch("Sex") + absdiff("Grade") + gwesp(0,fix=TRUE),
    eval.loglik=FALSE
  )
  
  egofit <- ergm.ego(
    fmh.ego ~ edges + degree(0:3) + nodefactor("Race") + nodematch("Race") + 
      nodefactor("Sex") + nodematch("Sex") + absdiff("Grade") + gwesp(0,fix=TRUE),
    popsize = network.size(faux.mesa.high),
    control = control.ergm.ego(
      ergm = control.ergm(
        parallel=1
      )
    )
  )
  
  # rbind(c(0, -2.9868, -0.154, 0.1162, -0.1869, 0.1431, -1.1494, -0.8786,
  #       -2.2698, -0.6101, 0.7227, -0.1537, 0.5557, -1.0783, 1.5875),
  #     coef(egofit))
  
  # rbind(
  #   coef(egofit),
  #   c(0, coef(fit))
  # )

  expect_equivalent(
    coef(egofit),
    c(0, coef(fit)),
    tolerance = 0.05
  )
})
