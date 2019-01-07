context("Testing degreedist()")

library(egor)


test_that("degreedist() gives correct results on custom data", {
  
  library(tibble)
  
  ## |----|--------|---|----------|-----------|
  ## | id | weight | x | # alters | alter a's |
  ## |----|--------|---|----------|-----------|
  ## | 1  | 2      | a | 1        | a         |
  ## | 2  | 1      | a | 2        | b, a      |
  ## | 3  | 1      | b | 1        | b         |
  ## | 4  | 2      | b | 2        | a, b      |
  ## |----|--------|---|----------|-----------|
  
  e <- egor(
    list(
      tibble(x="a"),
      tibble(x=c("b","a")),
      tibble(x="b"),
      tibble(x=c("a","b"))
    ),
    tibble(y=letters[c(1,1,2,2)]),
    ID.vars = list(ego="x", alter="x"),
    ego_design=list(~1,weights=c(2,1,1,2))
  )

  expect_equivalent(
    unclass(degreedist(e, plot=FALSE)), 
    c(1/2,1/2)
  )

  expect_equivalent(
    unclass(degreedist(e, plot=FALSE, by="y")), 
    rbind(c(2/3,1/3), c(1/3,2/3))
  )
})







# Tests using egor:egor32 data --------------------------------------------


data("egor32", package="egor")


test_that("degreedist() works on egor::egor32 data", {
  pdf(NULL)
  expect_silent(degreedist(egor32))
  dev.off()
})

test_that("degreedist() works on egor::egor32 data with `by=sex` (a factor)", {
  pdf(NULL)
  expect_silent(degreedist(egor32, by="sex"))
  dev.off()
})



# Tests using ergm::faux.mesa.high data -----------------------------------

data(faux.mesa.high, package="ergm")
fmh.ego <- as.egor(faux.mesa.high)

test_that("degreedist() works on data based on faux.mesa.high", {
  pdf(NULL)
  expect_silent(degreedist(fmh.ego))
  dev.off()
})

test_that("degreedist() works on data based on faux.mesa.high with `by=Sex`", {
  pdf(NULL)
  expect_silent(degreedist(fmh.ego, by="Sex"))
  dev.off()
})
