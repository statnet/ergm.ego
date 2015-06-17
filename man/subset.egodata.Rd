\name{subset.egodata}
\alias{subset.egodata}
\alias{[.egodata}
\title{
  Subsetting \code{\link{egodata}} Objects
}
\description{
  Returns subsets of \code{\link{egodata}} objects that meet conditions.
}
\usage{
\method{subset}{egodata}(x, subset, select, ...,
                         dup.action = c("make.unique", "fail", "number"))

\method{[}{egodata}(x, i, j, ..., dup.action=c("make.unique", "fail", "number"))
}

\arguments{
  \item{x}{
An \code{\link{egodata}} object.
}
  \item{subset, i}{
An expression (evaluated in the context of the \code{egos} table of
  \code{x} producing a logical, integer, or character vector indicating
  which egos to select (and, for the latter two, how many times).
}
  \item{select, j}{
A numeric or character vector specifying the columns of \code{egos} and \code{alters} to select.
}
  \item{\dots}{
Additional arguments, currently unused.
}
  \item{dup.action}{
What to do when an ego is referenced multiple times:
\describe{
\item{\code{"make.unique"}}{Construct new unique ego IDs using the
  \code{\link[base]{make.unique}} function
}
\item{\code{"fail"}}{Exit with an error.}
\item{\code{"number"}}{Number the egos consecutively in the order they
  were selected
}
}
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
Pavel N. Krivitsky
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
