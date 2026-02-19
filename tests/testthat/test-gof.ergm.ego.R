#  File tests/testthat/test-gof.ergm.ego.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2026 Statnet Commons
################################################################################
data("fmhfit", package="ergm.ego")


test_that("GOF='model' works", {
  # Computing
  expect_warning(z <- gof(fmhfit, GOF = "model"), NA)

  # Plotting
  expect_warning({
    pdf(NULL)
    plot(z)
    dev.off()
  }, NA)
})


test_that("GOF='degree' works", {
  # Computing
  expect_warning(z <- gof(fmhfit, GOF = "degree"), NA)

  # Plotting
  expect_warning({
    pdf(NULL)
    plot(z)
    dev.off()
  }, NA)
})


test_that("GOF='espartners' works", {
  # Computing
  expect_warning(z <- gof(fmhfit, GOF = "espartners"), NA)
  
  # Plotting
  expect_warning({
    pdf(NULL)
    plot(z)
    dev.off()
  }, NA)
})



test_that("GOF='espartners' works if `esp` term is in the model", {
  data("faux.mesa.high", package="ergm")
  edata <- as.egor(faux.mesa.high)
  fit <- ergm.ego(
    edata ~ edges + esp(1), 
    control = control.ergm.ego(
      ergm=control.ergm(MCMLE.maxit=2)
    )
  )

  expect_warning(gofobj <- gof(fit, GOF="espartners"), NA)
})
