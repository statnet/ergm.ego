context("Testing gof.ergm.ego() with `fmhfit` model")

data("fmhfit", package="ergm.ego")


test_that("GOF='model' works", {
  # Computing
  expect_silent(
    z <- gof(fmhfit, GOF="model")
  )
  
  # Plotting
  expect_silent({
    pdf(NULL)
    plot(z)
    dev.off()
  })
})


test_that("GOF='degree' works", {
  # Computing
  expect_silent(
    z <- gof(fmhfit, GOF="degree")
  )
  
  # Plotting
  expect_silent({
    pdf(NULL)
    plot(z)
    dev.off()
  })
})



test_that("GOF='espartners' works", {
  # Computing
  expect_silent(
    z <- gof(fmhfit, GOF="espartners")
  )
  
  # Plotting
  expect_silent({
    pdf(NULL)
    plot(z)
    dev.off()
  })
})






context("Temporary testing gof.ergm.ego() using model fit to GSS data")

model_path <- "~/Projects/statnet/ergm.ego-papers/sn-work/estimation/main-10k-cd.rds"

test_that("computing and plotting gof(GOF='espartners') works", {
  testthat::skip_if_not(file.exists(model_path), "Can't find the model file, skipping.")
  model <- readRDS(model_path)
  expect_silent(
    gobject <- gof(model, GOF="espartners")
  )
  
  expect_silent({
    pdf(NULL)
    plot(gobject)
    dev.off()
  })
})
