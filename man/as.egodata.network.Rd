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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An \code{\link{egodata}} object.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{as.network.egodata}}, which performs the inverse operation
  (though drops the ties).
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
