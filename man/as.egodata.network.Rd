%  File man/as.egodata.network.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{as.egodata.network}
\alias{as.egodata.network}
\title{
  Construct an Egocentric View of a Network
}
\description{
  Given a \code{\link[network]{network}} object, construct an
  \code{\link{egodata}} object representing a census of all the actors in
  the network. Used mainly for testing.
}
\usage{
\method{as.egodata}{network}(object, special.cols = c("na", "vertex.names"), ..., egoIDcol = "vertex.names")
}
\arguments{
  \item{object}{
    A \code{\link[network]{network}} object.
}
  \item{special.cols}{
    Vertex attributes that should not be copied to the \code{egos} and
    \code{alters} tables. Defaults to attributes special to the
    \code{\link[network]{network}} objects.
}
  \item{\dots}{
    Additional arguments, currently unused.
}
  \item{egoIDcol}{
    The name of the vertex attribute containg unique ego IDs. Defaults
    to \code{"vertex.names"}.
}
}
\value{
  An \code{\link{egodata}} object.
}
\author{
  Pavel N. Krivitsky
}
\seealso{
  \code{\link{as.network.egodata}}, which performs the inverse operation
  (though drops the ties).
}
\examples{
# See example(ergm.ego) and example(as.network.egodata).
}
\keyword{ datagen }
\keyword{ manip }
