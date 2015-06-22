\name{gof.ergm.ego}
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
  An object of class \code{\link{gofobject}}.
}
\author{
  Pavel N. Krivitsky
}
\seealso{For examples, see \code{\link{ergm.ego}}.}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
