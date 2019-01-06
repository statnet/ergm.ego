context("Test netsize.adj() on directed network")

data(sampson, package="ergm")

fmla <- samplike ~ netsize.adj(edges=ec, mutual=mc, transitiveties=tc, cyclicalties=cc) + 
  edges + mutual + transitiveties + cyclicalties

set.seed(0)
for(i in 1:10) {
  a <- rnorm(4)
  a[rbinom(4,1,.7)==0] <- 0 
  
  fmla <- as.formula(do.call(substitute, list(fmla, list(ec=a[1], mc=a[2], tc=a[3], cc=a[4]))))
  
  test_that(
    paste("it works for", format(frmla)),
    {
      expect_equivalent(
        0,
        as.vector(crossprod(c(-1,a),summary(fmla)))
      )
    }
  )
}











context("Test netsize.adj() on undirected network")


test_that("", {
  skip("for now")
  
  # Undirected
  data(faux.mesa.high)
  
  replicate(10, {
    a <- rnorm(2)
    a[rbinom(2,1,.5)==0] <- 0 
    
    fmla <- faux.mesa.high~netsize.adj(edges=ec, transitiveties=tc) + edges + transitiveties
    fmla <- as.formula(do.call(substitute, list(fmla, list(ec=a[1], tc=a[2]))))
    stopifnot(isTRUE(all.equal(0, as.vector(crossprod(c(-1,a),summary(fmla))))))
  })
  
})