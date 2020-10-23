#  File tests/EgoStat.tests.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################
library(ergm.ego)
library(ergm)

n <- 100
e <- 150
ds <- c(10,15,5,20)

y <- network.initialize(n, directed=FALSE)
y %v% "a" <- sample(1:3+6,n,replace=TRUE)
y %v% "b" <- sample(letters[1:4],n,replace=TRUE)
y %v% "c" <- sample(runif(10),n,replace=TRUE)
y %v% "d" <- runif(n)
y <- san(y~edges+degree(0:3), target.stats=c(e,ds))

y.e <- as.egodata(y)

f <- ~ edges +
  nodecov("a") + nodecov(c("a", "c")) + nodecov(c("a", "c", "d")) + nodecov(~c^2 + sin(d)) + offset(nodecov(~c^2 + sin(d))) + 
  
  nodefactor("a", levels=NULL) + nodefactor("a", levels=-1) + nodefactor("a", levels=1) + nodefactor("a", levels=-2) + nodefactor("b", levels=c("c")) + 
  nodefactor("b", levels=c("c", "a")) + nodefactor("b", levels=c("c", "e", "a")) + nodefactor("b", levels=c("c", "e")) + nodefactor(~a, levels = 2:3) + 
  nodefactor(c("b", "c"), levels = -(5:10)) + offset(nodefactor(c("b", "c"), levels = -(5:10))) + 
  
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
  
  nodemix("a") + nodemix("a", levels2=-1) + nodemix("a", levels2=1) + nodemix("b", levels = c("b"), levels2=1) + nodemix("a", levels2=-2) + nodemix("a", levels2=-(2:3)) + 
  nodemix(c("a", "b"), levels = -(2:3)) + nodemix(c("a", "b"), levels = -(2:3), levels2 = -(3:4)) + nodemix("b", levels = -2, levels2 = -3) + 
  offset(nodemix("b", levels = -2, levels2 = -3)) + 
  
  mm("a") + mm("a", levels2=~-1) + mm("a", levels2=-2) + mm("a", levels2=-(2:3)) + mm(~a>7) + mm(a~b) + mm(.~a) + offset(mm(.~a)) + mm("a", levels2 = 1) + 
  mm("b", levels = c("a", "c", "e")) + mm("b", levels = c("a", "c", "e"), levels2 = 3) +

  meandeg

f.y <- statnet.common::nonsimp_update.formula(f, y~.)
environment(f.y) <- globalenv()
f.y.e <- statnet.common::nonsimp_update.formula(f, y.e~.)
environment(f.y.e) <- globalenv()

stopifnot(all.equal(summary(f.y),summary(f.y.e)))

