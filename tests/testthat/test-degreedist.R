context("test-degreedist.R")

## |----|--------|---|----------|-----------|
## | id | weight | x | # alters | alter a's |
## |----|--------|---|----------|-----------|
## | 1  | 2      | a | 1        | a         |
## | 2  | 1      | a | 2        | b, a      |
## | 3  | 1      | b | 1        | b         |
## | 4  | 2      | b | 2        | a, b      |
## |----|--------|---|----------|-----------|

e <- egodata(data.frame(egoID=1:4, a=letters[c(1,1,2,2)]),
             alters=data.frame(egoID=c(1L,2L,2L,3L,4L,4L), a=letters[c(1,2,1,2,1,2)]),
             egoWt=c(2,1,1,2),
             egoIDcol="egoID")
test_that("weighted degreedist", {
  expect_equivalent(degreedist(e, plot=FALSE), c(0,1/2,1/2))
})

test_that("weighted degreedist by attribute", {
  expect_equivalent(degreedist(e, plot=FALSE, by="a"), rbind(c(0,2/3,1/3),
                                                             c(0,1/3,2/3)))
})
