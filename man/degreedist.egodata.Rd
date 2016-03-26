%  File man/degreedist.egodata.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{degreedist.egodata}
\alias{degreedist.egodata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting the degree distribution of an egocentric dataset
}
\description{
A function to plot a histogram of the degree distribution of actors in
the egocentric dataset, optionally broken down by group and/or compared
with a Bernoulli graph.
}
\usage{
degreedist.egodata(egodata, freq = FALSE, prob = !freq, by = NULL, brgmod = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egodata}{
    A \code{\link{egodata}} object.
  }
  \item{freq, prob}{
    Whether to plot the raw frequencies or the conditional proportions
    of the degree values. Defaults to the latter.
  }
  \item{by}{
    A character vector giving the name of a vertex attribute; if given,
    plots the frequences broken down by that attribute.
  }
  \item{brgmod}{
    Plot the range of predicted frequencies/probabilities according to a
    Bernoulli graph having the same expected density as the observed.
  }
}
\seealso{
  \code{\link{degreedist}}, \code{\link[ergm:summary.formula]{summary}}
}
\examples{
data(faux.mesa.high)
fmh.ego <- as.egodata(faux.mesa.high)

degreedist.egodata(fmh.ego,by="Grade",brgmod=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
