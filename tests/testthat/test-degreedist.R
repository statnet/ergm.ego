context("Testing degreedist()")

library(egor)
library(ergm.ego)


test_that("degreedist() works on egor::egor32 data", {
  data("egor32", package="egor")
  expect_silent(degreedist(egor32))
})

test_that("degreedist() works on egor::egor32 data with `by=sex` (a factor)", {
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
  degreedist(egor32, by="sex")

  egor32a <- egor32
  egor32a$woman <- egor32a$sex == "w"
  degreedist(egor32a, by="woman")
  
  data(faux.mesa.high)
  fmh.ego <- as.egor(faux.mesa.high)
  degreedist(fmh.ego)
  degreedist(fmh.ego, by="Sex")
  
  
  
  data(emon, package="network")
  w.ego <- as.egor(emon$Wichita)
  degreedist(w.ego)
  degreedist(w.ego, by="Location")
}

## |----|--------|---|----------|-----------|
## | id | weight | x | # alters | alter a's |
## |----|--------|---|----------|-----------|
## | 1  | 2      | a | 1        | a         |
## | 2  | 1      | a | 2        | b, a      |
## | 3  | 1      | b | 1        | b         |
## | 4  | 2      | b | 2        | a, b      |
## |----|--------|---|----------|-----------|

library(tibble)
e <- egor(list(tibble(x="a"),
               tibble(x=c("b","a")),
               tibble(x="b"),
               tibble(x=c("a","b"))),
          tibble(x=letters[c(1,1,2,2)]),
          ego_design=list(~1,weights=c(2,1,1,2)))

test_that("weighted degreedist", {
  expect_equivalent(unclass(degreedist(e, plot=FALSE)), c(1/2,1/2))
})

test_that("weighted degreedist by attribute", {
  expect_equivalent(unclass(degreedist(e, plot=FALSE, by="x")), rbind(c(2/3,1/3),
                                                                      c(1/3,2/3)))
})
