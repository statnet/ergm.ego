library(ergm.ego)

options(error=sessioninfo::session_info)

data(faux.mesa.high)
fmh.ego <- as.egor(faux.mesa.high)
fmhfit <- ergm.ego(
  fmh.ego ~ edges + degree(0:3) + 
    nodefactor("Race") + nodematch("Race")
  + nodefactor("Sex")+nodematch("Sex")
  + absdiff("Grade") + transitiveties, 
  popsize = network.size(faux.mesa.high),
  control = control.ergm.ego(
    ergm.control = control.ergm(parallel=2)
  )
)








# save
save(fmhfit, file="../data/fmhfit.rda")













if(FALSE) {
  # Complete fit
  fit <- ergm(
    faux.mesa.high ~ edges + degree(0:3) + 
      nodefactor("Race") + nodematch("Race")
    + nodefactor("Sex")+nodematch("Sex")
    + absdiff("Grade") + transitiveties, 
    control = control.ergm(parallel=2)
  )

fmhfit2 <- ergm.ego(
  fmh.ego ~ edges + degree(0:3) + 
    nodefactor("Race") + nodematch("Race")
  + nodefactor("Sex")+nodematch("Sex")
  + absdiff("Grade") + transitiveties, 
  popsize = network.size(faux.mesa.high),
  control = control.ergm.ego(
    ergm.control = control.ergm(parallel=2, init=c(0,coef(fit)))
  )
)

}
