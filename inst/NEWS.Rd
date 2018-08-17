%  File inst/NEWS.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2015-2018 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm.ego-package]{ergm.ego}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.

\section{Changes in version 0.4.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      
      \item \code{\link{mm}} term has been implemented.
      \code{\link[=degreedist.egodata]{degreedist}} now has an option to
      not plot, and returns the calculated degree distribution
      (invisibly, if plotting).
      
      \item \code{\link{offset}} terms are now handled.
      
      \item More \code{EgoStat} now handle more options that their
      \code{ergm} counterparts do.

      \item \code{ergm.ego}'s \code{ppopsize} control parameter and
      \code{simulate} method for \code{ergm.ego}'s \code{popsize}
      argument now take a data frame of egos to use as the
      pseudopopulation.

    }
  }
   
  
   \subsection{BUG FIXES}{
     \itemize{

       \item Package now works with \code{\link[=ergm-package]{ergm}}
       3.9.

       \item \code{\link[=degreedist.egodata]{degreedist}} now handles
       sampling weights correctly, and has been fixed in other ways.

       \item Bootstrap and jackknife now handle one-dimentional stats
       correctly.

       \item \code{mixingmatrix.egodata} now handles ego ID column names
       other than \code{vertex.names}. Thanks to Deven Hamilton for
       reporting this bug. Non-numeric ego IDs are also handled correctly.

       \item \code{mixingmatrix.egodata} no longer rounds the row
       probabilities before returning when called \code{rowprob=TRUE}.


     }
   }

  
   \subsection{OTHER USER-VISIBLE CHANGES}{
     \itemize{
       
       \item \code{\link{degreedist.egodata}} is now an
       \code{\link{egodata}} method of \code{\link[ergm]{degreedist}}.
       
     }
   }

  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }
}

\section{Changes in version 0.3.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item This is the initial public release.
    }
  }
   
  
   %% \subsection{BUG FIXES}{
   %%   \itemize{
   %%   }
   %% }

  
   %% \subsection{OTHER USER-VISIBLE CHANGES}{
   %%   \itemize{
   %%   }
   %% }

  % \subsection{C-LEVEL FACILITIES}{
  %   \itemize{
  %   }
  % }

  % \subsection{UTILITIES}{
  %   \itemize{
  %     \item Likewise, "UTILITIES" is probably unnecessary.
  %   }
  % }

  % \subsection{INSTALLATION}{
  %   \itemize{
  %     \item ... as is "INSTALLATION".
  %   }
  % }

  % \subsection{PACKAGE INSTALLATION}{
  %   \itemize{
  %     \item ditto.
  %   }
  % }
}

