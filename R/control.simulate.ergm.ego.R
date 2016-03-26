control.simulate.ergm.ego <- function(
  ppop.wt = c("round","sample"),
  SAN.control = control.san(),
  simulate.control = control.simulate(),
  ...){
  match.arg.pars <- c("ppop.wt")

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
