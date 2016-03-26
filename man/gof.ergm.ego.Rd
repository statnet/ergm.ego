%  File man/gof.ergm.ego.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\alias{gof.ergm.ego}
\title{Conduct Goodness-of-Fit Diagnostics on a Exponential Family
  Random Graph Model fit to Egocentrically Sampled Data}
\description{
  \code{\link{gof.ergm.ego}} implements the \code{\link[ergm]{gof}}
  method for \code{\link{ergm.ego}} fit objects.
}
\usage{
\method{gof}{ergm.ego}(object, ...,
                       GOF = c("model", "degree"),
                       control = control.gof.ergm(),
                       verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An \code{\link{ergm.ego}} fit.
  }
  \item{\dots}{
    Additional arguments, currently unused.
  }
  \item{GOF}{
    A string specifying the statistics whose goodness of fit is to be
    evaluated. Currently, only \dQuote{degree} and \dQuote{model} are
    implemented; see \code{\link[ergm]{gof}} documentation for details.
  }
  \item{control}{
    A list to control parameters, constructed using
    \code{\link{control.gof.formula}} or \code{\link{control.gof.ergm}}
    (which have different defaults).
  }
  \item{verbose}{
    Provide verbose information on the progress of the
    simulation.
  }
}
\value{
  An object of class \code{\link[ergm:gof]{gofobject}}.
}
\author{
  Pavel N. Krivitsky
}
\seealso{For examples, see \code{\link{ergm.ego}}.}
\keyword{ models }
