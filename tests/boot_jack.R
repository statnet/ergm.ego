#  File tests/boot_jack.R in package ergm.ego, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2015-2020 Statnet Commons
#######################################################################
library(ergm.ego)
data(faux.mesa.high)
fmh.ego <- as.egodata(faux.mesa.high)

head(fmh.ego)

set.seed(0)

## Two parameters
egofit <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                          popsize=network.size(faux.mesa.high))

egofit.boot <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="bootstrap"))

egofit.jack <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="jackknife"))

stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.boot), tolerance=0.02)))
stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.jack), tolerance=0.02)))

stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.boot), tolerance=0.05)))
stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.jack), tolerance=0.05)))


## One parameter
egofit <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high))

egofit.boot <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="bootstrap"))

egofit.jack <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="jackknife"))

stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.boot), tolerance=0.02)))
stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.jack), tolerance=0.02)))

stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.boot), tolerance=0.05)))
stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.jack), tolerance=0.05)))

