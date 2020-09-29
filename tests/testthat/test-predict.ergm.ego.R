context("Test predict.ergm.ego()")

test_that("it works", {
  data(faux.mesa.high, package="ergm")
  fmh.ego <- as.egodata(faux.mesa.high)
  egofit <- ergm.ego(
    fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
    +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
    popsize=network.size(faux.mesa.high)
  )
  expect_silent(
    predict(egofit)
  )
})
