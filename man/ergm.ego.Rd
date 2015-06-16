\name{ergm.ego}
\alias{ergm.ego}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit an ERGM to a \code{\link{egodata}} dataset.
}
\description{
  A wrapper to infer ERGM sufficient statistics, construct the
  pseudopopulation network, fit the model, and evaluate the standard
  errors.
}
\usage{
ergm.ego(formula, popsize = 1, offset.coef = NULL, ..., control = control.ergm.ego(), na.action = na.fail, do.fit = TRUE)
}
\arguments{
  \item{formula}{
    An \code{\link[ergm]{ergm}}-style formula, whose LHS is a
    \code{\link{egodata}} object.
}
  \item{popsize}{
    Size of the population from which the data were drawn; ERGM
    parameters would be shifted via the network size offset to estimate
    the parameters that would 
    have been gotten had the whole population been observed. Defaults to
    1, so that the shift is a function only of the pseudopopulation size.
}
  \item{offset.coef,\dots}{
    Additional arguments passed to \code{\link[ergm]{ergm}}.
  }
  \item{control}{
    The control list, constructed by \code{\link{control.ergm.ego}}.
  }
  \item{na.action}{
    What to do about missing covariate attributes.
  }
  \item{do.fit}{
    Whether to actually perform the fit. If set to \code{FALSE}, will
    estimate the statistics and their variance--covariance matrix, but
    not call \code{\link[ergm]{ergm}}.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class \code{ergm.ego} inheriting from \code{\link{ergm}},
  overriding the \code{ergm} estimates with \code{ergm.ego}.
}
\references{
Pavel N. Krivitsky and Martina Morris. Inference for Social Network
Models from Egocentrically-Sampled Data, with Application to
Understanding Persistent Racial Disparities in HIV Prevalence in the US.
Working paper. National Institute for Applied Statistics Research Australia, University of Wollongong, 2015(05-15). \url{http://niasra.uow.edu.au/publications/UOW190187.html}
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
\keyword{models}
