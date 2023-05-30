%  File inst/NEWS.Rd in package ergm.ego, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2015-2023 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


\section{Changes in version 1.1.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item \code{\link{ergm.ego}()} now has a \code{basis=} argument. So does \code{\link{simulate.ergm.ego}()}, for consistency (as an alias for the \code{popsize=} argument).

      \item \code{\link{simulate.ergm.ego}()}'s \code{popsize=} argument can now be a network object, enabling simulation from any starting network.
    }
  }

  \subsection{BUG FIXES}{
    \itemize{
      \item{
        \item \code{\link{gof.ergm.ego}(GOF="degree")} now handles the case in which the observed or simulated degree distribution is dense and the LHS network is small more gracefully.

        \item \code{\link{gof.ergm.ego}()} was scrambling the the order of ESP terms.

        \item \code{\link{simulate.ergm.ego}()} is now more robust to models with offsets and extreme \dQuote{dropped} statistics.

        \item \code{\link{ergm.ego}()} (via \code{\link{control.ergm.ego}(ppopsize=)}) and \code{\link{simulate.ergm.ego}(popsize=)} can once again take \code{\link{data.frame}}s and \code{\link{tibble}}s to specify the pseudopopulation network composition directly.
      }
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item \code{\link{simulate.ergm.ego}()} now preserves some of the attributes attached by \code{\link{simulate.ergm}()} to the statistics matrix, including \code{"monitored"}.

      \item \code{\link{simulate.ergm.ego}()} no longer supports \code{\link{ergm.ego}} objects fit under under \CRANpkg{ergm} < 4.
    }
  }
}


\section{Changes in version 1.0.1}{

  \subsection{BUG FIXES}{
    \itemize{
      \item{
        Documentation fixes, particularly for compatibility with \pkg{ergm} 4.2.
      }
      \item{
        Summary for \code{\link{ergm.ego}} fits now displays the original call rather than the instrumental \code{\link{ergm}()} call.
      }
    }
  }

  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        \code{\link{control.ergm.ego}()} praameter \code{ignore.max.alters=} now defaults to \code{TRUE}, since simulation studies (Krivitsky, et al. 2020) showed that they did more harm than good.
      }
    }
  }
}


\section{Changes in version 1.0.0}{

  \subsection{NEW FEATURES}{
    \itemize{
      \item{
        This package now uses the \pkg{\link[egor:egor-package]{egor}} package's \code{egor} class for data storage and manipulation. A converter \code{\link{as.egor.egodata}()} is provided.
      }
      \item{
        \code{\link{ergm.ego}()} now supports complex survey designs set on \code{\link[egor:egor-object]{egor}} objects.
      }
      \item{
        \code{\link{ergm.ego}()} and the summary methods can now fit triadic effects (\code{gwesp}, \code{esp}, \code{transitiveties}) when alter-alter ties are available.
      }
      \item{
        \code{\link{ergm.ego}()} can now handle missing alter attributes in some circumstances, and provided they are missing completely at random.
      }
      \item{
        A number of new egostats have been implemented, including \code{gwdegree}
      }
      \item{
        A number of improvements to the goodness-of-fit routines.
      }
      \item{
        \code{\link{snctrl}()} UI for specifying control parameters is supported.
      }
      \item{
        Curved ERGMs are now supported; this capability should be considered experimental, as uncertainty estimates have not been rigorously derived.
      }
      \item{
        For nonscaling statistics such as \code{\link{meandeg}}, standard errors can now be computed.
      }
      \item{
        Network size adjustment can now be disabled during fitting.
      }
    }
  }


  \subsection{BUG FIXES}{
    \itemize{
      \item{
        Various fixes to \code{\link[=degreedist.egor]{degreedist}()}, \code{\link[=mixingmatrix.egor]{mixingmatrix}()}, and other methods.
      }
    }
  }


  \subsection{OTHER USER-VISIBLE CHANGES}{
    \itemize{
      \item{
        The function that was previously \code{as.network.egodata()} for constructing an empty network having the same composition as the egocentric dataset has been superseded by \code{\link{template_network}()}.
      }
      \item{
        Manually specified pseudo-population is handled better.
      }
      \item{
        \code{\link[=degreedist.egor]{degreedist}()} method for egocentric data now defauts to not making plots.
      }
      \item{
        \code{\link[=mixingmatrix.egor]{mixingmatrix}()} method for egocentric data now returns a \code{\link{table}}.
      }
    }
  }
}


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

