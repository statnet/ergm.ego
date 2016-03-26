%  File man/control.simulate.ergm.ego.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{control.simulate.ergm.ego}
\alias{control.simulate.ergm.ego}
\title{
  Control parameters for \code{\link{simulate.ergm.ego}}.
}
\description{
  Constructs and checks the list of control parameters for simulation by \code{\link{simulate.ergm.ego}}.
}
\usage{
control.simulate.ergm.ego(ppop.wt = c("round","sample"),
                          SAN.control = control.san(),
                          simulate.control = control.simulate(),
                          ...)
}
\arguments{
    \item{ppop.wt}{
    Because each ego must be represented in the pseuodopopulation
    network an integral number of times, if the sample is weighted (or
    the target \eqn{|N'|} calculated from \code{ppopsize} and \code{ppopsize.mul} is not a
    multiple of the sample size), it may not be possible, for a finite
    \eqn{|N'|} to represent each ego exactly according to its
    relative weight, and \code{ppop.wt} controls how the fractional egos
    are allocated:
    \describe{
      \item{\code{"round"} (default)}{Rather than treating
	\code{ppopsize} as a hard setting, calculate \eqn{|N'|
	  w_i / w_\cdot} for each ego \eqn{i} and round it to the
	nearest integer. Then, the \eqn{|N'|} actually used will be the
	sum of these rounded freqencies.}
      \item{\code{"sample"}}{Resample in proportion to \eqn{w_i}.} 
    }

  }
  
  \item{SAN.control}{
    A list of control parameters for \code{\link[ergm]{san}} constructed
    by \code{\link[ergm]{control.ergm}}, called to construct a
    pseudopopulation network consistent with the data.
}
  \item{simulate.control}{
    A list of control parameters for \code{\link[ergm]{simulate.formula}} constructed
    by \code{\link[ergm]{control.simulate}}, called to simulate from the
    model fit.
}
  \item{\dots}{
    Not used at this time.
}
}
\value{
  A list with arguments as components.
}
\author{
  Pavel N. Krivitsky
}

\seealso{
  control.simulate, control.san
}
\keyword{ models }
