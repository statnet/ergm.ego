\name{as.egodata}
\alias{as.egodata}
\alias{as.egodata.data.frame}
\alias{as.egodata.network}
\alias{egodata}
\title{
  Convert to or Construct \link[=egodata.object]{\code{egodata} Objects}
}
\description{
  \code{\link{as.egodata}} is a generic function to construct
  \code{\link[=egodata.object]{egodata}} objects from a variety of
  sources. \code{\link{egodata}} function is the standard constructor,
  taking two data frames.
}
\usage{
as.egodata(object, \dots, egoIDcol = "egoID")

\method{as.egodata}{data.frame}(object, alters, egoWt = 1, ..., egoIDcol="egoID")

egodata(egos, alters, egoWt=1, ..., egoIDcol="egoID")

}
\arguments{
  \item{object, egos}{
    The object from which the egocentric data should be constructed. For
    the \code{\link{data.frame}} methods and \code{\link{egodata}}
    itself, a data frame containing the egos, as specified in
    \link[=egodata.object]{\code{egodata} object} documentation.
  }
  \item{\dots}{
    Additional arguments; currently unused.
}
  \item{egoIDcol}{
    Name of the column in the ego table containing the unique ego identifier.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An \link[=egodata.object]{\code{egodata} Object}.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Pavel N. Krivitsky
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
