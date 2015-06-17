# This is just an alias for the number of edges, to distinguish it in the summary tables.

InitErgmTerm.netsize.adj<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = NULL,
                      vartypes = NULL,
                      defaultvalues = list(),
                      required = NULL)
  
  list(name="edges", coef.names="netsize.adj", dependence=FALSE,
       minval = 0, maxval = network.dyadcount(nw,FALSE), conflicts.constraints="edges", pkgname="ergm")
}
