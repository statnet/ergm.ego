library(ergm.ego)
data(faux.mesa.high)
fmh.ego <- as.egor(faux.mesa.high)

set.seed(0)
ppop <- rbind(fmh.ego, fmh.ego)
ppop$.alters <- NULL
ppop$.alter_ties <- NULL
ppop <- as.data.frame(ppop)
egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
                   +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
                   popsize=network.size(faux.mesa.high),
                   control=control.ergm.ego(ppopsize=ppop))

stopifnot(isTRUE(all.equal(c(## -0.8407,
                             2.3393, 1.4686, 0.6323, 0.5287, -1.3603, -1.0454,
                             -2.4998, -0.7207, 0.833, -0.1823, 0.6357, -1.3513),
                           unname(coef(egofit))[-(1:2)],tolerance=.1)))

ppop <- fmh.ego[sample.int(nrow(fmh.ego), nrow(fmh.ego)*1.5, replace=TRUE),]
ppop$.alters <- NULL
ppop$.alter_ties <- NULL
ppop <- as.data.frame(ppop)
egosim <- simulate(egofit, popsize=ppop)
