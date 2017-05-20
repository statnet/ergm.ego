library(ergm.ego)
data(faux.mesa.high)
fmh.ego <- as.egor(faux.mesa.high)

head(fmh.ego)

## Two parameters
egofit <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                          popsize=network.size(faux.mesa.high))

egofit.boot <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="bootstrap"))

egofit.jack <- ergm.ego(fmh.ego~edges+nodematch("Sex"), 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="jackknife"))

stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.boot), tolerance=0.01)))
stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.jack), tolerance=0.01)))

stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.boot), tolerance=0.05)))
stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.jack), tolerance=0.05)))


## One parameter
egofit <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high))

egofit.boot <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="bootstrap"))

egofit.jack <- ergm.ego(fmh.ego~edges, 
                          popsize=network.size(faux.mesa.high), control=control.ergm.ego(stats.est="jackknife"))

stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.boot), tolerance=0.01)))
stopifnot(isTRUE(all.equal(coef(egofit),coef(egofit.jack), tolerance=0.01)))

stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.boot), tolerance=0.05)))
stopifnot(isTRUE(all.equal(vcov(egofit),vcov(egofit.jack), tolerance=0.05)))

