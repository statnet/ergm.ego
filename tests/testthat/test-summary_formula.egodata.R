
# Testing summary_formula.egodata() ---------------------------------------


test_that("summary statistics for egodata are identical to the complete network data", {
  data(faux.mesa.high)
  fmh.ego <- as.egodata(faux.mesa.high)
  nw.summ <- summary(faux.mesa.high~edges+degree(0:3)+nodematch("Race")+
                        nodematch("Sex")+absdiff("Grade")+nodemix("Grade"))
  ego.summ <- summary(fmh.ego~edges+degree(0:3)+nodematch("Race")+nodematch("Sex")+
                         absdiff("Grade")+nodemix("Grade"),
                       scaleto=network.size(faux.mesa.high))
  expect_equal(ego.summ, nw.summ)
})

