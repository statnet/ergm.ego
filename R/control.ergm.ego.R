control.ergm.ego <- function(
  ppopsize = c("samp", "pop"),
  ppopsize.mul = 1,
  ppop.wt = c("round","greedy","sample"),                                      
  stats.wt = c("data","ppop"),
  stats.est = c("asymptotic", "bootstrap", "jackknife", "naive"),
  boot.R = 10000,
  ergm.control = control.ergm(),
  ...){
  match.arg.pars <- c("stats.est", "ppop.wt", "stats.wt", if(!is.numeric(ppopsize)) "ppopsize")

  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))
  if(length(list(...))) stop("Unrecognized control parameter: ",arg,".")

  for(arg in match.arg.pars)
    control[arg]<-list(match.arg(control[[arg]][1],eval(formal.args[[arg]])))

  set.control.class("control.simulate.ergm.ego")
}
