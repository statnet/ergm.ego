control.simulate.ergm.ego <- function(SAN.control = control.san(),
                                      simulate.control = control.simulate(),
                                      ...){
  control<-list()
  formal.args<-formals(sys.function())
  formal.args[["..."]]<-NULL
  for(arg in names(formal.args))
    control[arg]<-list(get(arg))
  if(length(list(...))) stop("Unrecognized control parameter: ",arg,".")
  set.control.class("control.simulate.ergm.ego")
}
