%  File man/summary.statistics.egodata.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{summary.statistics.egodata}
\alias{summary.statistics.egodata}
\alias{summary}
\alias{summary.formula}
\title{
  Calculation of ERGM-style summary statistics for \code{\link{egodata}} objects.
}
\description{
  Used to calculate the specified network statistics inferred from a
  \code{\link{egodata}} object.
}
\usage{
\method{summary.statistics}{egodata}(object, ..., basis = NULL, 
                            individual = FALSE, scaleto = NULL)
}
\arguments{
  \item{object}{
    An \code{\link[ergm]{ergm}}-style formula with a \code{\link{egodata}}
    object as the LHS.

    For a list of currently implemented egocentric terms for the RHS, see
    \code{\link{ergm.ego-terms}}. 
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

  Pavel N. Krivitsky, Mark S. Handcock, and Martina Morris. Adjusting
  for Network Size and Composition Effects in Exponential-Family Random
  Graph Models. \emph{Statistical Methodology}, 2011, 8(4),
  319-339. \doi{10.1016/j.stamet.2011.01.005}

}
\author{
  Pavel N. Krivitsky
}

\seealso{
\code{\link[ergm]{summary.statistics}}, \code{\link[ergm]{summary.statistics.ergm}}
}
\examples{
data(faux.mesa.high)
fmh.ego <- as.egodata(faux.mesa.high)
(nw.summ <- summary(faux.mesa.high~edges+degree(0:3)+nodematch("Race")+
                    nodematch("Sex")+absdiff("Grade")+nodemix("Grade")))

(ego.summ <- summary(fmh.ego~edges+degree(0:3)+nodematch("Race")+nodematch("Sex")+
                     absdiff("Grade")+nodemix("Grade"),
                     scaleto=network.size(faux.mesa.high)))

stopifnot(isTRUE(all.equal(nw.summ,ego.summ)))
}
