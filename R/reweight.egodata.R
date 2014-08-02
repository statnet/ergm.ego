reweight.egodata <- function(x, g, gw){
  gw <- gw[order(gw[,1]),2]
  
  w0 <- x$egoWt

  # aggregate will sort by g.
  gw0 <- aggregate(w0~g, FUN=sum)

  gwm <- setNames(gw/gw0[,2],gw0[,1])
    
  x$egoWt <- x$egoWt*gwm[g]

  x
}

category.weights.egodata <- function(x, pop, by){
  x.g <- apply(x$egos[by],1,paste,collapse="\n")
  pop.g <- apply(pop$egos[by],1,paste,collapse="\n")

  pop.freq <- as.data.frame(table(pop.g))

  list(g=x.g, gw=pop.freq)
}
