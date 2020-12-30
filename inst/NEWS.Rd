%  File inst/NEWS.Rd in package ergm.ego, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2015-2020 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{NEWS}
\title{NEWS file for the \code{\link[=ergm.ego-package]{ergm.ego}} package}
\encoding{UTF-8}

% When a release is forked, copy the sections to a new release and
% delete unused sections.

%% \section{Changes in version X.Y.Z}{

%%   \subsection{NEW FEATURES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }


%%   \subsection{BUG FIXES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }


%%   \subsection{OTHER USER-VISIBLE CHANGES}{
%%     \itemize{
%%       \item{
%%       }
%%     }
%%   }

  %% \subsection{C-LEVEL FACILITIES}{
  %%   \itemize{
  %%   }
  %% }

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

%% }

\section{Changes in version 0.6.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        \code{\link{predict.ergm.ego}}, a \code{\link{predict}} method for \code{\link{ergm.ego}} has been implemented. (Thanks, Michał Bojanowski.)
      }
      \item{
        Nonscaling statistic \code{\link{meandeg}} has been added.
      }
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        \code{EgoStat.*} functions no longer need to be exported, reducing namespace pollution.
      }
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item{
        \code{\link{ergm.ego}} now detects when a coefficient has been
        dropped by \code{\link[ergm]{ergm}} due to the statistic having
        an extreme value and subsets the variance matrices accordingly.
      }
      \item{
        \code{\link{control.ergm.ego}} now calls \code{match.arg} on 
        \code{ppopsize} only if \code{ppopsize} is of class 
        \code{character}.  This allows \code{ppopsize} to be of class 
        \code{\link[network]{network}} when calling 
        \code{\link{control.ergm.ego}}.
      }
      \item{
        A more thorough search mechanism for \code{EgoStat.*} functions no longer requires them to be exported.
      }
    }
  }
}

\section{Changes in version 0.5.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
	\code{\link[ergm]{ergm}}'s new nodal attributes user interface has been extended to \code{ergm.ego}.
      }
      \item{
	\code{\link{mixingmatrix.egodata}} and \code{\link{degreedist.egodata}} now have an option to ignore sampling weights.
      }
      \item{
	Simulation frmo an \code{\link{ergm.ego}} fit now inherints the constraints.
      }
      \item{
	It is now possible to specify the (pseudo)population network temlate directly by passing it to \code{control$ppopsize}.
      }
      \item{
	It is now possible to infer main effects (\code{nodefactor} and \code{nodecov}) when the attribute has only been obseved on the egos.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
	A wide variety of minor bugs has been fixed. See commit log and issue tracker for details.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
	A number of robustifications have been made.
      }
      \item{
	\code{\link{ergm.ego}} now produces sensible error messages when terms have alter categories that egos do not.
      }
      \item{
	Chad Klumb has been added as a contributor.
      }
      \item{
	\code{\link{gof.ergm.ego}}'s default MCMC.interval is now the MCMC.interval of the ergm fit \emph{scaled} by the ratio between the fit's \code{MCMC.samplesize} and GoF control's \code{nsim}.
      }
      \item{
	\code{\link{gof.ergm.ego}} now only calculates GOF for degree values up to twice the highest observed in the data or 6, whichever is higher with an additional category to catch the higher values.
      }
    }
  }
}

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

