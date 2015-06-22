\name{ergm.ego-terms}
\alias{ergm.ego-terms}
\alias{terms-ergm.ego}
\alias{ergm.ego.terms}
\alias{terms.ergm.ego}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{EgoStat}
\alias{EgoStat.edges}
\alias{EgoStat.nodecov}
\alias{EgoStat.nodefactor}
\alias{EgoStat.nodematch}
\alias{EgoStat.nodemix}
\alias{EgoStat.absdiff}
\alias{EgoStat.degree}
\alias{EgoStat.degrange}
\alias{EgoStat.concurrent}
\alias{EgoStat.concurrentties}
\alias{EgoStat.degreepopularity}

\docType{methods}
\title{\code{\link[ergm]{ergm}} Terms Implemented for \code{\link[=egodata.object]{egodata}}}
\description{
  This page describes the \code{\link[ergm]{ergm}} terms (and hence network statistics)
  for which inference based on egocentrically sampled data is
  implemented in \code{ergm.ego} package. Other packages may add their own terms.

  The current recommendation for any package implementing additional
  egocentric calculator terms is to create a help file with a name or
  alias \code{ergm.egodata-terms}, so that
  \code{help("ergm.egodata-terms")} will list egocentric ERGM terms
  available from all loaded packages.
}

\section{Currently implemented egocentric statistics}{
  For each of these, please see \code{\link[ergm]{ergm-terms}} for
  meaning and parameters.

  \itemize{
    \item edges
    \item nodecov
    \item nodefactor
    \item nodematch
    \item nodemix
    \item absdiff
    \item degree
    \item degrange
    \item concurrent
    \item concurrentties
    \item degreepopularity
  }
}  

\seealso{\code{\link[ergm]{ergm-terms}}}
\keyword{models}
