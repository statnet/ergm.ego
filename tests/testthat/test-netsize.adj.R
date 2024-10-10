#  File tests/testthat/test-netsize.adj.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################
data(sampson, package="ergm")

for(i in 1:10) {
  set.seed(i)
  a <- rnorm(4)
  a[rbinom(4,1,.7)==0] <- 0

  fmla <- as.formula(do.call("substitute", list(
                                             samplike ~ netsize.adj(edges=ec, mutual=mc, transitiveties=tc, cyclicalties=cc) +
                                               edges + mutual + gwesp(0,fix=TRUE) + cyclicalties,
                                             list(ec=a[1], mc=a[2], tc=a[3], cc=a[4])
                                           )))
  s <- summary(fmla)

  test_that(paste("it works for a directed network with a <- ", paste(deparse(a), collapse=""), collapse=" "), {
    expect_equal(ignore_attr=TRUE,
      0,
      as.vector( crossprod( c(-1, a), s ) ),
      info = paste(names(s), "=", s)
    )
  })
}

data(faux.mesa.high, package="ergm")

for(i in 1:10) {
  set.seed(i)
  a <- rnorm(2)
  a[rbinom(2,1,.5)==0] <- 0

  f <- as.formula(do.call(substitute, list(
                                        faux.mesa.high~netsize.adj(edges=ec, transitiveties=tc) + edges + gwesp(0, fix=TRUE),
                                        list(ec=a[1], tc=a[2])
                                      )))
  s <- summary(f)

  test_that(paste("it works for an undirected network with a <- ", paste(deparse(a), collapse=""), collapse=" "), {
    expect_equal(ignore_attr=TRUE,
      0,
      as.vector(crossprod(c(-1, a), s)),
      info = paste(names(s), "=", s)
    )
  })
}
