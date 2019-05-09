#  File tests/EgoStat.tests.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2015-2018 Statnet Commons
#######################################################################


context("Testing EgoStat")

test_that("egostats are close to complete network stats", {
  # Test data
  set.seed(0)
  n <- 100
  e <- 150
  ds <- c(10,15,5,20)
  y <- network.initialize(n, directed=FALSE)
  y %v% "a" <- sample(1:3+6,n,replace=TRUE)
  y %v% "b" <- sample(letters[1:4],n,replace=TRUE)
  y <- san(y~edges+degree(0:3), target.stats=c(e,ds))
  y.e <- as.egor(y)
  
  
  
  f <- ~ edges +
    nodecov("a") +
    
    nodefactor("a", 0) + nodefactor("a", 1) + nodefactor("a", 2) +
    
    nodematch("a") + nodematch("a", TRUE) + nodematch("a", TRUE, 2) +
    
    absdiff("a") + absdiff("a", 2) +
    
    degree(0) + degree(3) + degree(0:6) +
    degree(0, by="a") + degree(3, by="a") + degree(0:6, by="a") +
    degree(0, by="a", homophily=TRUE) + degree(3, by="a", homophily=TRUE) + degree(0:6, by="a", homophily=TRUE) +
    
    degrange(0) + degrange(3) + degrange(0:6) +
    degrange(0) + degrange(3) + degrange(0:6) +
    degrange(0, by="a") + degrange(3, by="a") + degrange(0:6, by="a") +
    degrange(0, by="a", homophily=TRUE) + degrange(3, by="a", homophily=TRUE) + degrange(0:6, by="a", homophily=TRUE) +
    
    degrange(0,2) + degrange(3,5) + degrange(0:6,7) +
    degrange(0,2) + degrange(3,5) + degrange(0:6,7) +
    degrange(0,2, by="a") + degrange(3,5, by="a") + degrange(0:6,7, by="a") +
    degrange(0,2, by="a", homophily=TRUE) + degrange(3,5, by="a", homophily=TRUE) + degrange(0:6,7, by="a", homophily=TRUE) +
    
    concurrent + concurrent("a") +
    
    concurrentties + concurrentties("a") +
    
    degree1.5 +
    
    nodemix("a") + nodemix("a", base=1) + nodemix("a", base=2) + nodemix("a", base=2:3) +
    
    transitiveties + transitiveties("a") + esp(0:6) + gwesp(fix=FALSE) + gwesp(0.5, fix=TRUE) +
    
    mm("a") + mm("a", levels2=~-1) + mm("a", levels2=-2) + mm("a", levels2=-(2:3)) + mm(~a>7) + mm(a~b) + mm(.~a) +
    
    gwdegree(fix=FALSE) + gwdegree(0.5, fix=TRUE)
  
  f.y <- statnet.common::nonsimp_update.formula(f, y~.)
  # environment(f.y) <- globalenv()
  f.y.e <- statnet.common::nonsimp_update.formula(f, y.e~.)
  # environment(f.y.e) <- globalenv()
  
  expect_equivalent(
    as.vector(summary(f.y)),
    as.vector(summary(f.y.e))
  )
})



test_that("egostats with alter missing data are close to complete network stats", {
  
  ergm.ego:::long_test()
  
  # Test data
  set.seed(0)
  n <- 100
  e <- 150
  ds <- c(10,15,5,20)
  y <- network.initialize(n, directed=FALSE)
  y %v% "a" <- sample(1:3+6,n,replace=TRUE)
  y %v% "b" <- sample(letters[1:4],n,replace=TRUE)
  y <- san(y~edges+degree(0:3), target.stats=c(e,ds))
  y.e <- as.egor(y)
  

  f <- ~ edges +
    nodecov("a") +
    
    nodefactor("a", 0) + nodefactor("a", 1) + nodefactor("a", 2) +
    
    nodematch("a") + nodematch("a", TRUE) + nodematch("a", TRUE, 2) +
    
    absdiff("a") + absdiff("a", 2) +
    
    degree(0) + degree(3) + degree(0:6) +
    degree(0, by="a") + degree(3, by="a") + degree(0:6, by="a") +
    
    degrange(0) + degrange(3) + degrange(0:6) +
    degrange(0) + degrange(3) + degrange(0:6) +
    degrange(0, by="a") + degrange(3, by="a") + degrange(0:6, by="a") +
    
    degrange(0,2) + degrange(3,5) + degrange(0:6,7) +
    degrange(0,2) + degrange(3,5) + degrange(0:6,7) +
    degrange(0,2, by="a") + degrange(3,5, by="a") + degrange(0:6,7, by="a") +
    
    concurrent + concurrent("a") +
    
    concurrentties + concurrentties("a") +
    
    degree1.5 +
    
    nodemix("a") + nodemix("a", base=1) + nodemix("a", base=2) + nodemix("a", base=2:3) +
    
    transitiveties + esp(0:6) + gwesp(fix=FALSE) + gwesp(0.5, fix=TRUE) +
    
    mm("a") + mm("a", levels2=~-1) + mm("a", levels2=-2) + mm("a", levels2=-(2:3)) + mm(~a>7) + mm(a~b) + mm(.~a) +
    
    gwdegree(fix=FALSE) + gwdegree(0.5, fix=TRUE)
  
  
  replicate(30,{
    y.em <- y.e
    y.em$.alts <- lapply(y.em$.alts, function(a){
      N <- nrow(a)
      if(N){
        am <- a
        am[-sample.int(N,1),] <- NA
        am$.altID <- a$.altID
        am
      }else a
    })
    f.y.em <- statnet.common::nonsimp_update.formula(f, y.em~., from.new="y.em")
    summary(f.y.em)
  })->s
  
  # Non-varying
  f.y <- statnet.common::nonsimp_update.formula(f, y~., from.new="y")
  d <- sweep(s, 1, summary(f.y))
  novar <- apply(d, 1, sd) < sqrt(.Machine$double.eps)
  ok <- abs(d[novar,]) < sqrt(.Machine$double.eps)
  expect_true(
    all(ok),
    info = paste("Non-varying missing alter data estimate is off:", 
                 paste(rownames(d[novar,])[!ok], collapse=", ")
    )
  )
  
  # Varying
  pvals <- sapply(apply(d[!novar,], 1, t.test), "[[", "p.value")
  pval <- pchisq(-2*sum(log(pvals)), 2*sum(!novar), lower.tail=FALSE)
  expect_true(
    pval > 0.001, # Not very safe, since this test is stochastic.
    info = paste("Varying missing alter data estimate is off, p-value =", pval)
  )
})
