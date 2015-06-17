\name{as.network.egodata}
\alias{as.network.egodata}
%- Also NEED an '\alias' for EACH other topic documented here.
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
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{\link[network]{network}} object.
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
  \code{\link{as.egodata.network}}, which performs the inverse operation.
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
