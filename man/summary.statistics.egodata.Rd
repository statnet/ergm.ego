\name{summary.statistics.egodata}
\alias{summary.statistics.egodata}
\alias{summary}
\title{
  Claculation of ERGM-style summary statistics for \code{\link{egodata}} objects.
}
\description{
  Used to calculate the specified network statistics inferred from a
  \code{\link{egodata}} object.
}
\usage{
summary.statistics.egodata(object, ..., basis = NULL, individual = FALSE, scaleto = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    An \code{\link[ergm]{ergm}}-style formula with a \code{\link{egodata}}
    object as the LHS.
  }
  \item{\dots}{
    Not used at this time.
  }
  \item{basis}{
    An optional \code{\link{egodata}} object relative to which
    the statistics should be calculated.
  }
  \item{individual}{
    If \code{FALSE} (the default), calculate the estimated per-capita
    statistics, weighted according to the ego weights, then scale them
    up to a network of size \code{scaleto}.

    If \code{TRUE}, calculate each ego's individual contribution to the
    specified network statistics.
  }
  \item{scaleto}{
    Size of a hypothetical network to which to scale the
    statistics. Defaults to the number of egos in the dataset.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  If \code{individual==FALSE}, a named vector of statistics. If
    \code{individual==TRUE}, a matrix with a row for each ego, giving
    that ego's contribution to the network statistic.
}
\references{
  Pavel N. Krivitsky and Martina Morris. Inference for Social Network
  Models from Egocentrically-Sampled Data, with Application to
  Understanding Persistent Racial Disparities in HIV Prevalence in the
  US. Thechnical Report. National Institute for Applied Statistics
  Research Australia, University of Wollongong,
  2015(05-15). \url{http://niasra.uow.edu.au/publications/UOW190187.html}
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
