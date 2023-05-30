#  File tests/testthat/test-boot_jack.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2023 Statnet Commons
################################################################################

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

  expect_equal(ignore_attr=TRUE,
    coef(egofit),
    coef(egofit.boot),
    tolerance = 0.02
  )

  expect_equal(ignore_attr=TRUE,
    coef(egofit),
    coef(egofit.jack),
    tolerance = 0.02
  )

  expect_equal(ignore_attr=TRUE,
    vcov(egofit),
    vcov(egofit.boot),
    tolerance = 0.05
  )
  
  expect_equal(ignore_attr=TRUE,
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
  
  expect_equal(ignore_attr=TRUE,
    coef(egofit),
    coef(egofit.boot),
    tolerance = 0.02
  )
  
  expect_equal(ignore_attr=TRUE,
    coef(egofit),
    coef(egofit.jack),
    tolerance = 0.02
  )
  
  expect_equal(ignore_attr=TRUE,
    vcov(egofit),
    vcov(egofit.boot),
    tolerance = 0.05
  )
  
  expect_equal(ignore_attr=TRUE,
    vcov(egofit),
    vcov(egofit.jack),
    tolerance = 0.05
  )
})
