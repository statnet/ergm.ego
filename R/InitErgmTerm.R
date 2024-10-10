#  File R/InitErgmTerm.R in package ergm.ego, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2015-2024 Statnet Commons
################################################################################

#' A custom term giving a linear combination of change statistics
#' useful for network size adjustment.
#' 
#' @useDynLib ergm.ego
#' @noRd
InitErgmTerm.netsize.adj<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("edges", "mutual", "transitiveties", "cyclicalties") ,
                      vartypes = c("numeric", "numeric", "numeric", "numeric"),
                      defaultvalues = list(+1, 0, 0, 0),
                      required = c(FALSE, FALSE, FALSE, FALSE))

  if(with(a, edges==1 && mutual==0 && transitiveties==0 && cyclicalties==0)){
    list(name="edges", coef.names="netsize.adj", dependence=FALSE,
         pkgname="ergm")
  }else{
    if(!is.directed(nw) && with(a, mutual!=0 || (transitiveties!=0 && cyclicalties!=0))) stop("Network is undirected: mutuality cannot be specified and transitiveties and cyclicalties are identical.")
    list(name="netsize_adj", coef.names="netsize.adj", dependence=TRUE,
         inputs=with(a, c(edges, mutual, transitiveties, cyclicalties)))
  }
}
