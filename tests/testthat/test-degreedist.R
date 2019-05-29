#  File tests/testthat/test-degreedist.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2019 Statnet Commons
#######################################################################
context("test-degreedist.R")

## |----|--------|---|----------|-----------|
## | id | weight | x | # alters | alter a's |
## |----|--------|---|----------|-----------|
## | 1  | 2      | a | 1        | a         |
## | 2  | 1      | a | 2        | b, a      |
## | 3  | 1      | b | 1        | b         |
## | 4  | 2      | b | 2        | a, b      |
## |----|--------|---|----------|-----------|

e <- egodata(data.frame(egoID=1:4, x=letters[c(1,1,2,2)]),
             alters=data.frame(egoID=c(1L,2L,2L,3L,4L,4L), x=letters[c(1,2,1,2,1,2)]),
             egoWt=c(2,1,1,2),
             egoIDcol="egoID")
test_that("weighted degreedist", {
  expect_equivalent(unclass(degreedist(e, plot=FALSE)), c(1/2,1/2))
})

test_that("weighted degreedist by attribute", {
  expect_equivalent(unclass(degreedist(e, plot=FALSE, by="x")), rbind(c(2/3,1/3),
                                                                      c(1/3,2/3)))
})

test_that("weighted degreedist with weights disabled", {
  expect_equivalent(unclass(degreedist(e, plot=FALSE, weight=FALSE)), c(1/2,1/2))
})

test_that("weighted degreedist by attribute with weights disabled", {
  expect_equivalent(unclass(degreedist(e, plot=FALSE, by="x", weight=FALSE)), rbind(c(1/2,1/2),
                                                                                    c(1/2,1/2)))
})
