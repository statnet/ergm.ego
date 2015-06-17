\name{as.egodata}
\alias{as.egodata}
\alias{as.egodata.data.frame}
\alias{as.egodata.network}
\alias{egodata.object}
\alias{egodata}
\title{
  Convert to or Construct \code{\link{egodata}} Objects
}
\description{
  \code{\link{as.egodata}} is a generic function to construct
  \code{\link{egodata}} objects from a variety of
  sources. \code{\link{egodata}} function is the standard constructor,
  taking two data frames.
}
\usage{
as.egodata(object, \dots, egoIDcol = "egoID")

\method{as.egodata}{data.frame}(object, alters, egoWt = 1, ..., egoIDcol="egoID")

egodata(egos, alters, egoWt=1, ..., egoIDcol="egoID")

}
\arguments{
  \item{object, egos}{
    The object from which the egocentric data should be constructed. For
    the \code{\link{data.frame}} methods and \code{\link{egodata}}
    itself, a data frame containing at least the column named in
    \code{egoIDcol}, whose values must all be unique. Other columns
    contain information about the egos.
  }

  \item{alters}{ A data frame containing at least the column named in
    \code{egoIDcol}, whose values do not have to be unique, and not
    every ego must be represented. Other columns contain information
    about the alters.  }

  \item{egoWt}{A vector of the same length as number of rows in
    \code{egos} or \code{object}, containing the relative sampling
    weight of each ego.}

  \item{\dots}{
    Additional arguments; currently unused.
}
  \item{egoIDcol}{
    Name of the column in the ego table containing the unique ego identifier.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An \code{\link{egodata}} object. The object is a list containing the
  following elements:

    \item{egos}{
      A data frame with one row for each ego, containing at least the column named in
      \code{egoIDcol}, and other columns
      containing attributes of the egos.
    }

    \item{alters}{ A data frame containing at least the column named in
      \code{egoIDcol}, and other columns containing attributes of the alters.  }
    
    \item{egoWt}{A vector of the same length as the number of egos,
      containing the relative sampling weight of each ego.}
    
    \item{egoIDcol}{
      Name of the column in the ego table containing the unique ego identifier.
    }
  
}
\references{
%% ~put references to the literature/web site here ~
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
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
