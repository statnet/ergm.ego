#  File R/reweight.egor.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################
## reweight.egor <- function(x, g, gw){
##   gw <- gw[order(gw[,1]),2]
  
##   w0 <- x$egoWt

##   # aggregate will sort by g.
##   gw0 <- aggregate(w0~g, FUN=sum)

##   gwm <- setNames(gw/gw0[,2],gw0[,1])
    
##   x$egoWt <- x$egoWt*gwm[g]

##   x
## }

## category.weights.egor <- function(x, pop, by){
##   x.g <- apply(x$ego[by],1,paste,collapse="\n")
##   pop.g <- apply(pop$ego[by],1,paste,collapse="\n")

##   pop.freq <- as.data.frame(table(pop.g))

##   list(g=x.g, gw=pop.freq)
## }
