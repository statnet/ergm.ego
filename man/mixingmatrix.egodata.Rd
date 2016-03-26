%  File man/mixingmatrix.egodata.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{mixingmatrix.egodata}
\alias{mixingmatrix.egodata}
\title{
Summarizing the mixing among groups in an egocentric dataset
}
\description{
A function to return counts of how often a ego of each group nominates
an alter of each group.
}
\usage{
mixingmatrix.egodata(egodata, attrname, rowprob = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{egodata}{
    A \code{\link{egodata}} object.
  }
  \item{attrname}{
    A character vector containing the name of the network attribute
    whose mixing matrix is wanted.
  }
  \item{rowprob}{
    Whether the counts should be normalized by row sums. That is,
    whether they should be proportions conditional on the ego's group.
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  A matrix with a row and a column for each level of \code{attrname}.

  Note that, unlike \code{\link[network]{mixingmatrix}}, what is counted
  are \emph{nominations}, not ties. This means that under an egocentric
  census, the diagonal of \code{mixingmatrix.egodata} will be twice that
  returned by \code{\link[network]{mixingmatrix}} for the original undirected
  network.
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }

%% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link[network]{mixingmatrix}}, \code{\link[ergm]{nodemix}},
  \code{\link[ergm.ego]{summary}} method for egocentric data
}
\examples{
data(faux.mesa.high)
fmh.ego <- as.egodata(faux.mesa.high)

(mm <- mixingmatrix(faux.mesa.high,"Grade"))
(mm.ego <- mixingmatrix.egodata(fmh.ego,"Grade"))

stopifnot(isTRUE(all.equal({tmp<-unclass(mm$matrix); diag(tmp) <- diag(tmp)*2;
tmp}, mm.ego, check.attributes=FALSE)))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%% \keyword{ ~kwd1 }
%% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
