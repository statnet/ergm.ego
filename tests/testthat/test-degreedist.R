#  File tests/testthat/test-degreedist.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2022 Statnet Commons
################################################################################

## |----|--------|---|----------|-----------|
## | id | weight | x | # alters | alter a's |
## |----|--------|---|----------|-----------|
## | 1  | 2      | a | 1        | a         |
## | 2  | 1      | a | 2        | b, a      |
## | 3  | 1      | b | 1        | b         |
## | 4  | 2      | b | 2        | a, b      |
## |----|--------|---|----------|-----------|

e <- egor(
  alters=tibble(x=c("a","b","a","b","a","b"),
                i=c(1L,2L,2L,3L,4L,4L)),
  egos=tibble(i=seq_len(4),x=letters[c(1,1,2,2)],
              w=c(2,1,1,2)),
  ID.vars = list(ego="i"),
  ego_design=list(weights="w")
)

test_that("degreedist() gives correct results on custom data", {
  expect_equal(
    c(unclass(degreedist(e, plot=FALSE))),
    c(1/2,1/2),
    ignore_attr=TRUE
  )

  expect_equal(
    unclass(degreedist(e, plot=FALSE, by="x")), 
    rbind(c(2/3,1/3), c(1/3,2/3)),
    ignore_attr=TRUE
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

test_that("weighted degreedist with weights disabled", {
  expect_equal(ignore_attr=TRUE,unclass(degreedist(e, plot=FALSE, weight=FALSE)), c(1/2,1/2))
})

test_that("weighted degreedist by attribute with weights disabled", {
  expect_equal(ignore_attr=TRUE,unclass(degreedist(e, plot=FALSE, by="x", weight=FALSE)), rbind(c(1/2,1/2),
                                                                               c(1/2,1/2)))
})
