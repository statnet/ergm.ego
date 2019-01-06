context("Testing mixingmatrix() with egor::egor32 data")

data(egor32, package="egor")

varnames <- c("sex", "age")


for( v in varnames ) {
  test_that(
    paste0("Sum of mixing matrix for ", v, " is equal to number of alters"), {
      expect_silent(mm <- mixingmatrix(egor32, v))
      expect_equal(sum(mm), sum(vapply(egor32$.alts, nrow, numeric(1))))
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
