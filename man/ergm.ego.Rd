\name{ergm.ego}
\alias{ergm.ego}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Inference for Exponential-Family Random Graph Models based on
  Egocentrically Sampled Data
}
\description{
  A wrapper around the \code{\link[ergm]{ergm}} to fit an ERGM to an \code{\link[=egodata-object]{egodata}}.
}
\usage{
ergm.ego(formula, popsize, ppopsize = popsize, offset.coef = NULL, na.action = na.fail, stats.est = c("asymptotic", "bootstrap", "jackknife", "naive"), stats.wt = c("data", "ppop"), ppop.wt = c("round", "greedy", "sample"), R = 10000, ..., control = control.ergm(), do.fit = TRUE)
}
\arguments{
  \item{formula}{
    An \R \code{\link{formula}} object, of the form
    \code{e ~ <model terms>},
    where \code{e} is a \code{\link{egodata}} object. See
    \code{\link[ergm]{ergm}} for details and examples. 
  }
  \item{popsize}{
    The size \eqn{|N|} of the finite population network from which the egocentric
    sample was taken; only affects the shift in the coefficients of the
    terms modeling the overall propensity to have ties.
  }
  \item{ppopsize, ppop.wt}{
    Parameters for constructing the pseudopopulation
    network. \code{ppopsize} \eqn{|N'|} is the target size of the pseudopopulation network used for the
    \code{\link[ergm]{ergm}} call. It must be at least equal to the sample
    size, and, particularly if there are sampling weights in the data,
    must be bigger.

    Because each ego must be represented in the pseuodopopulation
    network an integral number of times, if the sample is weighted (or \code{ppopsize} is not a
    multiple of the sample size), it may not be possible, for a finite
    \code{ppopsize} to represent each ego exactly according to its
    relative weight, and \code{ppop.wt} controls how the fractional egos
    are allocated:
    \description{
      \item{\code{"round"} (default)}{Rather than treating
	\code{ppopsize} as a hard setting, calculate \eqn{|N'|
	  w_i / w_\cdot} for each ego \eqn{i} and round it to the
	nearest integer. Then, the \eqn{|N'|} actually used will be the
	sum of these rounded freqencies.}
      \item{\code{"sample"}}{Resample in proportion to \eqn{w_i}.} 
    }
  }
  \item{offset.coef}{
    A vector of coefficients for the offset terms.
  }
  \item{na.action}{
    How to handle missing actor attributes in egos or alters, when the
    terms need them.
  }
  \item{stats.est, R}{
    Method to be used to estimate the ERGM's sufficient statistics and their
    variance:
    \description{
      \item{\code{"asymptotic"} (default)}{Delta method, as derived by Krivitsky and
	Morris (2015), assuming the ego weights are sampled alongside the egos.}
      \item{\code{"bootstrap"}}{Nonparametric bootstrap with bias correction, resampling
	egos, using \code{R} replications.}
      \item{\code{"jackknife"}}{Jackknife.}
      \item{\code{"naive"}}{"Naive" estimator, assuming that weights are fixed.}
      }
}
  \item{stats.wt}{
    Weight assigned to each ego's contribution to the ERGM's sufficient
    statistic:
    \description{
      \item{\code{"data"} (default)}{Use weights \eqn{|N'|
	  w_i / w_\cdot} for each ego \eqn{i} as in the data.}
      \item{\code{"ppop"}}{Use weights ultimately used in the
	pseudopopulation network.} 
    }
}
  \item{\dots}{
    Additional arguments passed to \code{\link[ergm]{ergm}}.
}
  \item{control}{
    A \code{\link[ergm]{control.ergm}} control list.
}
  \item{do.fit}{
    Whether to actually call \code{\link[ergm]{ergm}}
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An object of class \code{ergm.ego} inheriting from
    \code{\link[ergm]{ergm}}, with the following additional or
    overridden elements:
    \item{\code{v}}{Estimated variance-covariance matrix of the estimate
      of the sufficient statistics}
    \item{\code{m}}{Per-capita estimate of the sufficient statistics}
    \item{\code{egodata}}{The egodata object passed}
    \item{\code{popsize}, \code{ppopsize}}{Population network size and
      pseudopopulation size used, respectively}
    \item{\code{coef}}{The coefficients adjusted for network size.}
    \item{\code{covar}}{Pseudo-MLE estimate of the variance-covariance
      matrix of the parameter estimates under repeated egocentric
      sampling}
    \item{\code{ergm.coef}, \code{ergm.covar},
      \code{ergm.est.cov}, \code{est.cov}}{Respectively, the
      unadjusted ERGM parameter estimates, the variance-covariance matrix of
      parameter estimates under the ERGM superpopulation process, and two
      copies of the variance-covariance matrix of uncertainty due to MCMLE
      approximation.}
    \item{\code{DtDe}}{Estimated Jacobian of the expectation of the
      sufficient statistics with respect to the model parameters}
}
\references{
  
Pavel N. Krivitsky and Martina Morris. Inference for Social Network
    Models from Egocentrically-Sampled Data, with Application to
    Understanding Persistent Racial Disparities in HIV Prevalence in the
    US. Thechnical Report. National Institute for Applied Statistics
    Research Australia, University of Wollongong,
    2015(05-15). \url{http://niasra.uow.edu.au/publications/UOW190187.html}

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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (formula, popsize, ppopsize = popsize, offset.coef = NULL, 
    na.action = na.fail, stats.est = c("asymptotic", "bootstrap", 
        "jackknife", "naive"), stats.wt = c("data", "ppop"), 
    ppop.wt = c("round", "greedy", "sample"), R = 10000, ..., 
    control = control.ergm(), do.fit = TRUE) 
{
    stats.est <- match.arg(stats.est)
    stats.wt <- match.arg(stats.wt)
    egodata <- get(as.character(formula[[2]]), envir = environment(formula))
    message("Constructiong pseudopopulation network.")
    popnw <- as.network(egodata, ppopsize, scaling = ppop.wt)
    if (network.size(popnw) != ppopsize) {
        message("Note: Constructed network has size ", network.size(popnw), 
            ", different from requested ", ppopsize, ". Estimation should not be meaningfully affected.")
        ppopsize <- network.size(popnw)
    }
    w <- switch(stats.wt, data = egodata$egoWt, ppop = tabulate(popnw \%n\% 
        "ego.inds", nbins = nrow(egodata)))
    stats <- summary(remove.offset.formula(formula), individual = TRUE)
    tmp <- na.action(cbind(w, stats))
    w <- tmp[, 1]
    stats <- tmp[, -1, drop = FALSE]
    n <- length(w)
    wmean <- function(w, s) {
        w <- w/sum(w)
        colSums(s * w)
    }
    m <- wmean(w, stats)
    if (stats.est == "bootstrap") {
        m.b <- t(replicate(R, {
            i <- sample.int(length(w), replace = TRUE)
            wmean(w[i], stats[i, , drop = FALSE])
        }))
        m <- m - (colMeans(m.b) - m)
    }
    else if (stats.est == "jackknife") {
        m.j <- t(sapply(seq_len(n), function(i) {
            wmean(w[-i], stats[-i, , drop = FALSE])
        }))
        m <- n * m - (n - 1) * colMeans(m.j)
    }
    v <- switch(stats.est, bootstrap = cov(m.b), jackknife = (n - 
        1)/n * crossprod(sweep(m.j, 2, colMeans(m.j))), naive = {
        w <- w/sum(w)
        crossprod(sweep(stats, 2, m, "-") * sqrt(w))/(1 - sum(w^2)) * 
            sum(w^2)
    }, asymptotic = .asymptotic.var(stats, w)/length(w))
    out <- list(v = v, m = m, formula = formula, egodata = egodata, 
        ppopsize = ppopsize, popsize = popsize)
    if (do.fit) {
        ergm.formula <- ergm.update.formula(formula, popnw ~ 
            offset(edges) + ., from.new = "popnw")
        ergm.fit <- ergm(ergm.formula, target.stats = m * ppopsize, 
            offset.coef = c(-log(ppopsize/popsize), offset.coef), 
            ..., eval.loglik = FALSE, control = control)
        coef <- coef(ergm.fit)
        oi <- offset.info.formula(ergm.formula)
        DtDe <- -ergm.fit$hessian[!oi$theta, !oi$theta, drop = FALSE]
        vcov <- matrix(NA, length(coef), length(coef))
        vcov[!oi$theta, !oi$theta] <- solve(DtDe/ppopsize) \%*\% 
            v \%*\% solve(DtDe/ppopsize)
        rownames(vcov) <- colnames(vcov) <- names(coef)
        out <- c(out, list(coef = coef, covar = vcov, est.cov = ergm.fit$est.cov, 
            ergm.coef = ergm.fit$coef, ergm.covar = ergm.fit$covar, 
            ergm.est.cov = ergm.fit$est.cov, DtDe = DtDe))
        for (name in names(out)) ergm.fit[[name]] <- out[[name]]
        out <- ergm.fit
    }
    class(out) <- c("ergm.ego", "ergm")
    out
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
