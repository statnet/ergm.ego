#  File tests/table_ppop.R in package ergm.ego, part of the Statnet suite
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

set.seed(0)
ppop <- rbind(fmh.ego$egos, fmh.ego$egos)
egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
                   +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
                   popsize=network.size(faux.mesa.high),
                   control=control.ergm.ego(ppopsize=ppop))

stopifnot(isTRUE(all.equal(c(## -0.8407,
                             2.3393, 1.4686, 0.6323, 0.5287, -1.3603, -1.0454,
                             -2.4998, -0.7207, 0.833, -0.1823, 0.6357, -1.3513),
                           unname(coef(egofit))[-(1:2)],tolerance=.1)))

ppop <- fmh.ego$egos[sample.int(nrow(fmh.ego$egos), nrow(fmh.ego$egos)*1.5, replace=TRUE),]
egosim <- simulate(egofit, popsize=ppop)
