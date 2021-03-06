% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prcbr.R
\name{prcbr}
\alias{prcbr}
\title{Profiling estimation methods for RCBR models}
\usage{
prcbr(
  formula,
  b0,
  data,
  logL = TRUE,
  omethod = "BFGS",
  lo = -Inf,
  up = Inf,
  ...
)
}
\arguments{
\item{formula}{is of the extended form enabled by the \pkg{Formula} package.
In the Cosslett, or current status, model the formula takes the form 
\code{y ~ v | z} where \code{v} is the covariate designated to have coefficient
one, and \code{z} is another covariate or group of covariates that are assumed
fixed coefficients that are to be estimated.}

\item{b0}{is either an initial value of the parameter for the Z covariates
or a matrix of such values, in which case optimization occurs over this discrete
set, when there is only one covariate then b0 is either scalar, or a vector.}

\item{data}{data frame for formula variables}

\item{logL}{if logL is TRUE the log likelihood is optimized, otherwise
a GEE score criterion is minimized.}

\item{omethod}{optimization method for \code{optim}, default "BFGS".}

\item{lo}{lower bound(s) for the parameter domain}

\item{up}{upper bound(s) for the parameter domain}

\item{\dots}{other arguments to be passed to \code{rcbr.fit} to control fitting.}
}
\value{
a list comprising the components:
\describe{
\item{bopt}{output of the optimizer for the profiled parameters beta}
\item{fopt}{output of the optimizer for the random coefficients eta}
}
}
\description{
Profile likelihood and (GEE) score methods for estimation of random coefficient binary
response models.  This function is a wrapper for \code{rcbr} that uses the offset
argument to implement estimation of additional fixed parameters.  It may be useful to
restrict the domain of the optimization over the profiled parameters, this can be
accomplished, at least for box constraints by setting \code{omethod = "L-BFGS-B"}
and specifying the \code{lo} and \code{up} accordingly.
}
