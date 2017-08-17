context("Dummy testing of testthat infrastructure")

test_that("there is ergm.ego()", {
  expect_true(exists("ergm.ego", where="package:ergm.ego"))
})