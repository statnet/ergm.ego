%  File man/as.network.egodata.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{as.network.egodata}
\alias{as.network.egodata}
\title{
  Construct an Empty ``Template'' Network Consistent with an Egocentric Sample
}
\description{
  Taking a \code{\link{egodata}} object, constructs a
  \code{\link[network]{network}} object with no edges whose
  vertices have the attributes of the egos in the dataset, replicating
  the egos as needed, and taking into accounts their sampling weights.
}
\usage{
\method{as.network}{egodata}(x, N, scaling = c("round", "sample"), ...)
}
\arguments{
  \item{x}{
    A \code{\link{egodata}} object.
  }
  \item{N}{
    The target number of vertices the output network should have.
  }
  \item{scaling}{
    If \code{\link{egodata}} contains weights or 
    \code{N} is not a multiple of number of egos in the sample, it may not be possible, for a finite
    \code{N} to represent each ego exactly according to its
    relative weight, and \code{scaling} controls how the fractional egos
    are allocated:
    \describe{
      \item{\code{"round"} (default)}{Rather than treating
	\code{N} as a hard setting, calculate \eqn{N
	  w_i / w_\cdot} for each ego \eqn{i} and round it to the
	nearest integer. Then, the \code{N} actually used will be the
	sum of these rounded freqencies.}
      \item{\code{"sample"}}{Resample in proportion to \eqn{w_i}.} 
    }
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
}
\value{
  A \code{\link[network]{network}} object.
}
\author{
  Pavel N. Krivitsky
}

\seealso{
  \code{\link{as.egodata.network}}, which performs the inverse operation.
}
\examples{

data(faux.mesa.high)
summary(faux.mesa.high, print.adj = FALSE)

fmh.ego <- as.egodata(faux.mesa.high)

# Same actor attributes
fmh.template <- as.network(fmh.ego, N=network.size(faux.mesa.high))
summary(fmh.template, print.adj = FALSE)

# Twice the actors, same distribution
fmh2.template <- as.network(fmh.ego, N=2*network.size(faux.mesa.high))
summary(fmh2.template, print.adj = FALSE)
}
\keyword{ manip }
