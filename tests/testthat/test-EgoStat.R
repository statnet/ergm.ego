#  File tests/testthat/test-EgoStat.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2021 Statnet Commons
################################################################################
library(ergm.ego)
library(purrr)
library(dplyr)

# Test data
set.seed(0)
n <- 100
e <- 150
ds <- c(10,15,5,20)

y <- network.initialize(n, directed=FALSE)
y %v% "a" <- sample(1:3+6,n,replace=TRUE)
aM <- matrix(FALSE, 3, 3)
aM[1,1] <- aM[1,3] <- TRUE
y %v% "b" <- sample(letters[1:4],n,replace=TRUE)
y %v% "c" <- sample(runif(10),n,replace=TRUE)
y %v% "d" <- runif(n)
y <- san(y~edges+degree(0:3), target.stats=c(e,ds))

y.e <- as.egor(y)


test_that("egostats are close to complete network stats", {
  f <- ~ edges +
    nodecov("a") + nodecov(c("a", "c")) + nodecov(c("a", "c", "d")) + nodecov(~c^2 + sin(d)) + offset(nodecov(~c^2 + sin(d))) + 
    
    nodefactor("a", levels=NULL) + nodefactor("a", levels=-1) + nodefactor("a", levels=1) + nodefactor("a", levels=-2) + nodefactor("b", levels=c("c")) + 
    nodefactor("b", levels=c("c", "a")) + nodefactor("b", levels=c("c", "e", "a")) + nodefactor("b", levels=c("c", "e")) + nodefactor(~a, levels = 2:3) + 
    nodefactor(c("b", "c"), levels = -(5:10)) + offset(nodefactor(c("a", "b"), levels = -(5:10))) +
    
    nodematch("a") + nodematch("a", TRUE) + nodematch("a", TRUE, levels=2) + nodematch(c("a", "b"), TRUE, levels = 5:10) + nodematch(~2*a) + offset(nodematch(~a^2)) + 
    
    absdiff("a") + absdiff("a", 2) + absdiff(~a + 2*c - exp(d)) + 
    absdiff(function(x) if(is(x, "data.frame")) x[["a"]] + 2*x[["c"]] - exp(x[["d"]]) else (x %v% "a") + 2*(x %v% "c") - exp(x %v% "d")) + 
    
    degree(0) + degree(3) + degree(0:6) +
    degree(0, by="a") + degree(3, by="a") + degree(0:6, by="a") +
    degree(0, by="a", homophily=TRUE) + degree(3, by="a", homophily=TRUE) + degree(0:6, by="a", homophily=TRUE) + degree(0:6, by="a", homophily=FALSE, levels = -1) + 
    degree(0:6, by="a", homophily=FALSE, levels=1) + degree(0:6, by=c("a", "b"), homophily=FALSE, levels = -1) + degree(0:6, by=c("a", "b"), homophily=FALSE, levels=1) +
    
    degrange(0) + degrange(3) + degrange(0:6) +
    degrange(0) + degrange(3) + degrange(0:6) +
    degrange(0, by="a") + degrange(3, by="a") + degrange(0:6, by="a") + degrange(0:6, by=~b) + degrange(0:6, by= function(x) if(is(x, "data.frame")) x[["b"]] else x %v% "b") + 
    degrange(0, by="a", homophily=TRUE) + degrange(3, by="a", homophily=TRUE) + degrange(0:6, by="a", homophily=TRUE) + offset(degrange(0:6, by="a", homophily=TRUE)) + 
    degrange(0:6, by="a", homophily=FALSE, levels=2) + offset(degrange(0:6, by="a", homophily=FALSE, levels=-2)) + 
    degrange(0:6, by=c("a","b"), homophily=FALSE, levels=3) + offset(degrange(0:6, by=c("a","b"), homophily=FALSE, levels=-3)) + 
    
    degrange(0,2) + degrange(3,5) + degrange(0:6,7) +  
    degrange(0,2) + degrange(3,5) + degrange(0:6,7) +
    degrange(0,2, by="a") + degrange(3,5, by="a") + degrange(0:6,7, by="a") + degrange(0:6,7, by=~b) + 
    degrange(0,2, by="a", homophily=TRUE) + degrange(3,5, by="a", homophily=TRUE) + degrange(0:6,7, by="a", homophily=TRUE) +
    
    concurrent + concurrent("a") + concurrent("a", levels=2:3) + concurrent("a", levels=3) + concurrent(c("a", "b"), levels=2:3) + 
    concurrent(c("a", "b"), levels=4) + offset(concurrent(~a)) + 
    
    concurrentties + concurrentties("a") + offset(concurrentties(function(x) if(is(x, "data.frame")) x[["a"]]else (x %v% "a"))) + 
    concurrentties("a", levels=2:3) + concurrentties("a", levels=-1) + concurrentties("b", levels=c("a", "c")) + concurrentties("b", levels=c("a", "b", "e")) +
    concurrentties(c("a", "b"), levels=2:3) + concurrentties(c("a", "b"), levels=-3) +
    
    degree1.5 + offset(degree1.5) + 
    
    nodemix("a", levels2=TRUE) + nodemix("a", levels2=-1) + nodemix("a", levels2=1) + nodemix("b", levels = c("b"), levels2=1) + nodemix("a", levels2=-2) + nodemix("a", levels2=-(2:3)) +
    nodemix(c("a", "b"), levels = -(2:3)) + nodemix(c("a", "b"), levels = -(2:3), levels2 = -(3:4)) + nodemix("b", levels = -2, levels2 = -3) + nodemix("a",levels2=aM) +
    offset(nodemix("b", levels = -2, levels2 = -3)) + 
    
    mm("a", levels2=TRUE) + mm("a", levels2=~-1) + mm("a", levels2=-2) + mm("a", levels2=-(2:3)) + mm(~a>7) + mm(a~b) + mm(.~a) + offset(mm(.~a)) + mm("a", levels2 = 1) +
    mm("b", levels = c("a", "c", "e")) + mm("b", levels = c("a", "c", "e"), levels2 = 3) +
    esp(0:6) + gwesp(fix=FALSE) + gwesp(0.5, fix=TRUE) +
    gwdegree(0.5, fix=TRUE) +
    
    meandeg +

    transitiveties + transitiveties("a") + gwdegree(fix=FALSE)

  f.y <- statnet.common::nonsimp_update.formula(f, y~.)
  f.y.e <- statnet.common::nonsimp_update.formula(f, y.e~.)
  s.e <- summary(f.y.e)

  isMD <- names(s.e)=="meandeg"

  expect_equal(c(s.e), summary(f.y))
  
  expect_true(all(is.na(vcov(s.e)[!isMD, isMD])))
  expect_true(all(is.na(vcov(s.e)[isMD, !isMD])))
  expect_false(is.na(vcov(s.e)[isMD, isMD]))

  expect_false(anyNA(vcov(s.e)[!isMD, !isMD]))

  s.e2 <- s.e*2

  expect_equal(as.vector(s.e2)[!isMD], as.vector(s.e)[!isMD]*2)
  expect_equal(as.vector(s.e2)[isMD], as.vector(s.e)[isMD])

  expect_equal(vcov(s.e2)[!isMD, !isMD], vcov(s.e)[!isMD, !isMD] * 4)
  expect_equal(vcov(s.e2)[isMD, isMD], vcov(s.e)[isMD, isMD])

  expect_warning(s.e3 <- s.e*seq_along(s.e), ".*attempting to scale a nonscalable ego statistic.*")
  expect_equal(as.vector(s.e3)[!isMD], as.vector(s.e)[!isMD]*(seq_along(s.e))[!isMD])
  expect_equal(as.vector(s.e3)[isMD], as.vector(s.e)[isMD])

  expect_equal(vcov(s.e3)[!isMD, !isMD], t(vcov(s.e)[!isMD, !isMD] * (seq_along(s.e))[!isMD]) * (seq_along(s.e))[!isMD])
  expect_equal(vcov(s.e3)[isMD, isMD], vcov(s.e)[isMD, isMD])
})

test_that("scaling and nonscaling egostats are combined correctly", {
  expect_error(summary(y.e~edges + meandeg, individual=TRUE), ".onscaling statistic detected.*meaningless.*")
  expect_silent(summary(y.e~meandeg, scaleto=1))
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
    
    esp(0:6) + gwesp(fix=FALSE) + gwesp(0.5, fix=TRUE) +
    
    mm("a") + mm("a", levels2=~-1) + mm("a", levels2=-2) + mm("a", levels2=-(2:3)) + mm(~a>7) + mm(a~b) + mm(.~a) +
    
    gwdegree(0.5, fix=TRUE) +

    meandeg +

    transitiveties + gwdegree(fix=FALSE)

  alter_l <- alters_by_ego(y.e)
  replicate(30,{
    y.em <- y.e
    y.em$alter <- alter_l %>% map(function(a){
      N <- nrow(a)
      if(N){
        am <- a
        keep <- sample.int(N,1)
        am$a[-keep] <- NA_real_
        am$b[-keep] <- NA_character_
        am
      }else a
    }) %>% bind_rows
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
