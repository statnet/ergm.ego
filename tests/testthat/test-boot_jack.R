#  File tests/boot_jack.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################

context("Test if bootstrap and jackknife methods of suffstat estimation give same results")

data(faux.mesa.high)
fmh.ego <- as.egor(faux.mesa.high)


test_that("equality of results for a two-parameter model", {
  
  set.seed(0)
  
    egofit <- ergm.ego(fmh.ego ~ edges+nodematch("Sex"), 
                       popsize=network.size(faux.mesa.high))
    egofit.boot <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                            popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="bootstrap"))
    egofit.jack <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                            popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="jackknife"))

  expect_equivalent(
    coef(egofit),
    coef(egofit.boot),
    tolerance = 0.02
  )

  expect_equivalent(
    coef(egofit),
    coef(egofit.jack),
    tolerance = 0.02
  )

  expect_equivalent(
    vcov(egofit),
    vcov(egofit.boot),
    tolerance = 0.05
  )
  
  expect_equivalent(
    vcov(egofit),
    vcov(egofit.jack),
    tolerance = 0.05
  )
})





test_that("equality of results for a one-parameter model", {
  egofit <- ergm.ego(fmh.ego~edges, 
                     popsize=network.size(faux.mesa.high))
  
  egofit.boot <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="bootstrap"))
  
  egofit.jack <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="jackknife"))
  
  expect_equivalent(
    coef(egofit),
    coef(egofit.boot),
    tolerance = 0.02
  )
  
  expect_equivalent(
    coef(egofit),
    coef(egofit.jack),
    tolerance = 0.02
  )
  
  expect_equivalent(
    vcov(egofit),
    vcov(egofit.boot),
    tolerance = 0.05
  )
  
  expect_equivalent(
    vcov(egofit),
    vcov(egofit.jack),
    tolerance = 0.05
  )
})
