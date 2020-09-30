context("Test predict.ergm.ego()")

test_that("it just works for model without offsets", {
  data(faux.mesa.high, package="ergm")
  fmh.ego <- as.egodata(faux.mesa.high)
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
  data("faux.magnolia.high", package="ergm")
  fmhego <- as.egodata(faux.magnolia.high)
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
