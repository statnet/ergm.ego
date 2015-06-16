\name{egodata.object}
\alias{egodata.object}
\title{
  Description of Egocentric Sample Data Objects
}
\description{
  An object of class \code{\link{egodata}} is a container for the ego
  and alter information comprising an egocentric sample, as well as metadata
}
\value{
  The object is a list containing the following elements:

  \item{egoIDcol}{The name of the column containing unique identifiers
    of respondents (egos) in the sample.}
  
  \item{egos}{A data frame containing at least the column named in
    \code{egoIDcol}, whose values must all be unique. Other columns
    contain information about the egos.}
  
  \item{alters}{A data frame containing at least the column named in
    \code{egoIDcol}, whose values do not have to be unique, and not
    every ego must be represented. Other columns contain information
    about the alters.}

  \item{egoWt}{A vector of the same length as number of rows in
    \code{egos}, containing the relative sampling weight of each ego.}

  The object has class \code{"egodata"}.
}
\seealso{
  \code{\link{egodata}}, \code{\link{summary.statistics.egodata}}, 
}
