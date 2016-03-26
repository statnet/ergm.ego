%  File man/ergm.ego-terms.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2016 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
\alias{EgoStat.mean.age}
\alias{netsize.adj}
\alias{InitErgmTerm.netsize.adj}

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
  For each of these, please see their respective package's
  \code{ergm-terms} help for meaning and parameters. The simplest way to
  do this is usually via \code{? TERM}.
  
  \describe{
    \item{Special-purpose terms:}{
      \describe{
	\item{\code{netsize.adj}}{A special-purpose term equivalent to
	  \code{\link[ergm]{edges}}, to house the network-size
	  adjustment offset. This term is added to the model
	  automatically and should not be used in the model formula
	  directly.
	}
      }
    }
      
    \item{\code{\link[ergm:ergm-terms]{ergm}}:}{
      \itemize{
	\item \code{edges}
	\item \code{nodecov}
	\item \code{nodefactor}
	\item \code{nodematch}
	\item \code{nodemix}
	\item \code{absdiff}
	\item \code{degree}
	\item \code{degrange}
	\item \code{concurrent}
	\item \code{concurrentties}
	\item \code{degreepopularity}
      }
    }
  }
}

\seealso{\code{\link[ergm]{ergm-terms}}}
\keyword{models}
