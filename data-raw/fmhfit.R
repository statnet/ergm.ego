library(ergm.ego)

data(faux.mesa.high)
fmh.ego <- egor::as.egor(faux.mesa.high)
fmhfit <- ergm.ego(
  fmh.ego ~ edges + degree(0:3) + 
    nodefactor("Race") + nodematch("Race")
  + nodefactor("Sex")+nodematch("Sex")
  + absdiff("Grade") + transitiveties, 
  popsize = network.size(faux.mesa.high),
  verbose = 3,
  control = control.ergm.ego(
    ergm.control = control.ergm(
      parallel=2,
      seed = 666
      )
  )
)


# save
usethis::use_data(fmhfit, overwrite=TRUE)








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
