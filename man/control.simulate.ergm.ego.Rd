\name{control.simulate.ergm.ego}
\alias{control.simulate.ergm.ego}
\title{
  Control parameters for \code{\link{simulate.ergm.ego}}.
}
\description{
  Constructs and checks the list of control parameters for simulation by \code{\link{simulate.ergm.ego}}.
}
\usage{
control.simulate.ergm.ego(SAN.control = control.san(), simulate.control = control.simulate(), ...)
}
\arguments{
  \item{SAN.control}{
    A list of control parameters for \code{\link[ergm]{san}} constructed
    by \code{\link[ergm]{control.ergm}}, called to construct a
    pseudopopulation network consistent with the data.
}
  \item{simulate.control}{
    A list of control parameters for \code{\link[ergm]{simulate.formula}} constructed
    by \code{\link[ergm]{control.simulate}}, called to simulate from the
    model fit.
}
  \item{\dots}{
    Not used at this time.
}
}
\value{
  A list with arguments as components.
}
\author{
  Pavel N. Krivitsky
}

\seealso{
  control.simulate, control.san
}
\keyword{ models }
