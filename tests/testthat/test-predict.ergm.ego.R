#  File tests/testthat/test-predict.ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2021 Statnet Commons
################################################################################

library(ergm.ego)

test_that("it just works for model without offsets", {
  data(faux.mesa.high, package="ergm")
  fmh.ego <- as.egor(faux.mesa.high)
  egofit <- ergm.ego(
    fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
    +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
    popsize=network.size(faux.mesa.high)
  )
  expect_silent(
    p <- predict(egofit)
  )
  expect_true(all(is.finite(p$p)))
})


test_that("it just works for model with offsets", {
  data("faux.mesa.high", package="ergm")
  fmhego <- as.egor(faux.mesa.high)
  fit <- ergm.ego(
    fmhego ~ edges 
    + nodefactor("Grade")
    + nodematch("Grade", diff=T)
    + offset(nodematch("Sex",
                       diff = TRUE,
                       levels = c(1, 2))),
    offset.coef = rep(-Inf, 2)
  )
  expect_silent(
    p <- predict(fit) # data frame
  )
  expect_true(all(is.finite(p$p)))
})
