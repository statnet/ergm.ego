%  File man/simulate.ergm.ego.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{simulate.ergm.ego}
\alias{simulate.ergm.ego}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Simulate from a \code{\link{ergm.ego}} fit.
}
\description{
  A wrapper around \code{\link[ergm]{simulate.formula}} to simulate
  networks from an ERGM fit using \code{\link{ergm.ego}}.
}
\usage{
\method{simulate}{ergm.ego}(object, nsim = 1, seed = NULL, 
               popsize=if(object$popsize==1) object$ppopsize else object$popsize, 
               control = control.simulate.ergm.ego(), ..., verbose = FALSE)
}
\arguments{
  \item{object}{
    An \code{\link{ergm.ego}} fit.
}
  \item{nsim}{
    Number of realizations to simulate.
}
  \item{seed}{
    Random seed.
}
  \item{popsize}{
    Network size to which to scale the model for simulation.
}
  \item{control}{
    A \code{\link{control.simulate.ergm.ego}} control list.
}
  \item{\dots}{
    Additional arguments passed to \code{\link[ergm]{san}} and \code{\link[ergm]{simulate.formula}}.
}
  \item{verbose}{
    Verbosity of output.
}
}
\value{
  The ouput has the same format (with the same options) as
  \code{\link[ergm]{simulate.formula}}. If \code{statsonly=TRUE} is
  passed, an additional attribute, \code{"ppopsize"} is set, giving the
  actual size of the network reconstructed, when the \code{pop.wt}
  control parameter is set to \code{"round"} and \code{"popsize"} is not
  a multiple of the egocentric sample size or the sampling weights.
}
\references{
  Pavel N. Krivitsky and Martina Morris. Inference for Social Network
  Models from Egocentrically-Sampled Data, with Application to
  Understanding Persistent Racial Disparities in HIV Prevalence in the
  US. Thechnical Report. National Institute for Applied Statistics
  Research Australia, University of Wollongong,
  2015(05-15). \url{http://niasra.uow.edu.au/publications/UOW190187.html}

  Pavel N. Krivitsky, Mark S. Handcock, and Martina Morris. Adjusting
  for Network Size and Composition Effects in Exponential-Family Random
  Graph Models. \emph{Statistical Methodology}, 2011, 8(4),
  319-339. \doi{10.1016/j.stamet.2011.01.005}
  
}
\author{
  Pavel N. Krivitsky
}
\seealso{
  \code{\link[ergm]{simulate.formula}}, \code{\link[ergm]{simulate.ergm}}
}

\examples{
data(faux.mesa.high)
fmh.ego <- as.egodata(faux.mesa.high)
egofit <- ergm.ego(fmh.ego~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
                         +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"), 
                          popsize=network.size(faux.mesa.high))
colMeans(egosim <- simulate(egofit, popsize=500,nsim=200,
statsonly=TRUE, control=control.simulate.ergm.ego(
                    simulate.control=control.simulate.formula(MCMC.burnin=2e6))))
colMeans(egosim)/attr(egosim,"ppopsize")*network.size(faux.mesa.high)
summary(faux.mesa.high~edges+degree(0:3)+nodefactor("Race")+nodematch("Race")
                         +nodefactor("Sex")+nodematch("Sex")+absdiff("Grade"))
}
\keyword{models}
