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
