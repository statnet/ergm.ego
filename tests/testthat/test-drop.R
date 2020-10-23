library(ergm.ego)
nw <- network.initialize(20,directed=FALSE)

# No homophilous ties; heterophilous density of 1/2.
nw %v% "a" <- rep(1:2, each=10)
nw[1:10,11:20] <- 0:1

test_that("dropped ergm terms", {
  out.coef <- coef(ergm.ego(as.egodata(nw)~edges+nodematch("a")))
  out.coef <- c(sum(out.coef[1:2]),out.coef[3])
  expect_equivalent(coef(ergm(nw~edges+nodematch("a"))), out.coef, tolerance=0.01)
})
