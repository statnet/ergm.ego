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



test_that("GOF='espartners' works if `esp` term is in the model", {
  data("faux.mesa.high", package="ergm")
  edata <- as.egor(faux.mesa.high)
  fit <- ergm.ego(
    edata ~ edges + esp(1), 
    control = control.ergm.ego(
      ergm=control.ergm(MCMLE.maxit=2)
    )
  )
  
  expect_silent(
    gofobj <- gof(fit, GOF="espartners")
  )
})
