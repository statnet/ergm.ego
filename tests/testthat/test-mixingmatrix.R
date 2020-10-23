data(egor32, package="egor")

varnames <- c("sex", "age")


for( v in varnames ) {
  test_that(
    paste0("Sum of mixing matrix for ", v, " is equal to number of alters"), {
      expect_silent(mm <- mixingmatrix(egor32, v))
      expect_equal(sum(mm), nrow(egor32$alter))
    })
}


for( v in varnames ) {
  test_that(
    paste0(
      "Ordering of rows and columns for ", 
      sQuote(v), 
      " is identical to ordering of levels"
    ), {
      x <- dplyr::pull(egor32, v)
      mm <- mixingmatrix(egor32, v)
      expect_identical( rownames(mm), levels(x))
      expect_identical( colnames(mm), levels(x))
    })  
}

rm(varnames, v)


set.seed(666)
egodf <- dplyr::data_frame(
  id_ego = seq(1, 5),
  num = rnorm(max(id_ego)),
  int = as.integer(sample(1:3, max(id_ego), replace=TRUE)),
  ch = sample(letters[1:3], max(id_ego), replace=TRUE),
  fac = factor(sample(LETTERS[1:3], max(id_ego), replace=TRUE), levels=LETTERS[1:3])
)
degs <- rpois(nrow(egodf), 2)
alterdf <- dplyr::data_frame(
  id_ego = rep(egodf$id_ego, degs),
  id_alter = unlist(lapply(degs[degs>0], function(x) seq(1, x)))
) %>%
  dplyr::mutate(
    num = rnorm(n()),
    int = as.integer(sample(1:3, n(), replace=TRUE)),
    ch = sample(letters[1:4], n(), replace=TRUE),
    fac = factor(sample(LETTERS[1:4], n(), replace=TRUE), levels=LETTERS[4:1])
  )
edata <- egor::egor(alterdf, egodf, ID.vars=list(ego="id_ego", alter="id_alter"))

varnames <- c("int", "ch", "fac")
for ( v in varnames ) {
  test_that(
    paste0(
      "Sum of mixing matrix for ", 
      sQuote(v), 
      " is equal to the total number of alters"
    ), {
      expect_silent(mm <- mixingmatrix(edata, v))
      expect_equal(sum(mm), nrow(edata$alter))
    })  
}

# test_that("Rows of mm are properly ordered for factors", {
#   x <- dplyr::pull(edata, "fac")
#   expect_silent( mm <- mixingmatrix(edata, "fac"))
#   expect_identical( rownames(mm), levels(x))
# })

