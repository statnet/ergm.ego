#  File tests/testthat/test-table_ppop.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################
test_that("estimation works", {
  
  skip("Not yet fixed")
  
  data(faux.mesa.high)
  fmh.ego <- as.egor(faux.mesa.high)
  ppop <- rbind(fmh.ego, fmh.ego)
  ppop$.alts <- NULL
  ppop$.aaties <- NULL
  ppop <- as.data.frame(ppop)
  
  set.seed(0)
  expect_silent(
    egofit <- ergm.ego(
      fmh.ego ~ edges + degree(0:3) + nodefactor("Race") + nodematch("Race")
      + nodefactor("Sex") + nodematch("Sex") + absdiff("Grade"), 
      popsize = network.size(faux.mesa.high),
      control = control.ergm.ego(
        ppopsize=ppop,
        ergm = control.ergm(
          MCMLE.maxit=2
        )
      )
    )
  )
})



test_that("coefs are equal", {
  skip("Not yet fixed")
  
  est <- coef(egofit)
  tgt <-   c(## -0.8407,
    2.3393, 1.4686, 0.6323, 0.5287, -1.3603, -1.0454,
    -2.4998, -0.7207, 0.833, -0.1823, 0.6357, -1.3513)[-(1:2)]
  
  # plot(est, tgt); abline(a=0, b=1, lty=2)
  
  expect_equal(ignore_attr=TRUE,
    unname(est)[-(1:2)],
    tgt,
    tolerance = 0.1
  )
})




test_that("simulating works", {
  
  skip("Not yet fixed")
  
  ppop <- fmh.ego[sample.int(nrow(fmh.ego), nrow(fmh.ego)*1.5, replace=TRUE),]
  ppop$.alts <- NULL
  ppop$.aaties <- NULL
  ppop <- as.data.frame(ppop)
  
  expect_silent(
    egosim <- simulate(egofit, popsize=ppop)
  )
})


