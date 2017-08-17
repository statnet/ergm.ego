context("Testing degreedist()")

library(egor)
library(ergm.ego)


test_that("degreedist() works on egor::egor32 data", {
  data("egor32", package="egor")
  expect_silent(degreedist(egor32))
})

test_that("degreedist() works on egor::egor32 data with `by=sex`", {
  data("egor32", package="egor")
  expect_silent(degreedist(egor32, by="sex"))
})

test_that("degreedist() works on data based on faux.mesa.high", {
  data(faux.mesa.high, package="ergm")
  fmh.ego <- as.egor(faux.mesa.high)
  expect_silent(degreedist(fmh.ego))
})

test_that("degreedist() works on data based on faux.mesa.high with `by=Sex`", {
  data(faux.mesa.high, package="ergm")
  fmh.ego <- egor::as.egor(faux.mesa.high)
  expect_silent(degreedist(fmh.ego, by="Sex"))
})





if(FALSE) {
  debugonce(degreedist)
  
  data("egor32", package="egor")
  degreedist(egor32)
  degreedist(egor32, by="age")
  table(sapply(egor32$.alts, nrow))
  summary(egor32 ~ degree(0:19, by="sex"))
  
  
  
  data(faux.mesa.high)
  fmh.ego <- as.egor(faux.mesa.high)
  table(sapply(fmh.ego$.alts, nrow))
  summary(fmh.ego ~ degree(0:13, by="Sex"))
  s <- summary(
    fmh.ego ~ edges + degree(0:13, by="Sex"),
    scaleto=network.size(faux.mesa.high)
  )
  degreedist(fmh.ego)
  degreedist(fmh.ego, by="Sex")
  
  
  
  data(emon, package="network")
  w.ego <- as.egor(emon$Wichita)
  table(vapply(w.ego$.alts, nrow, numeric(1)))
  summary(w.ego ~ degree(0:18))
  summary(w.ego ~ degree(0:18, by="Location"))
  
  degreedist(w.ego, by="Location")
}
