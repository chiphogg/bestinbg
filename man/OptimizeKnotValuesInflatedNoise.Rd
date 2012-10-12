\name{OptimizeKnotValuesInflatedNoise}
\alias{OptimizeKnotValuesInflatedNoise}
\title{Optimize knot values with inflated noise level.}
\usage{
  OptimizeKnotValuesInflatedNoise(bgr, data, sigma)
}
\arguments{
  \item{bgr}{(\code{\link{background}} object) The
  estimated background (and associated hyperparameters)}

  \item{data}{(\code{\link{dataset}} object) The data whose
  background to estimate}
}
\value{
  (\code{\link{background}} object) A copy of \code{bgr},
  but with \code{Y(bgr)} optimized.
}
\description{
  We artificially inflate the noise level to avoid getting
  stuck in local optima.
}

