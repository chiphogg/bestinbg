\name{OptimizeCC}
\alias{OptimizeCC}
\title{Optimizes cc (i.e., the spline values; 'c' in Fischer et al.) with a given
level of experimental noise.}
\usage{
  OptimizeCC(x, y, cc = (0 * xi), xi, p.bgr, sigma, lambda,
    eta, Phi = spline.matrix(x = x, spline.x = xi))
}
\arguments{
  \item{x}{The x-values of the measured datapoints.}

  \item{y}{The y-values of the measured datapoints.}

  \item{cc}{Initial guess for the spline knot y-values ('c'
  in Fischer et al.).}

  \item{xi}{The spline knot x-values.}

  \item{p.bgr}{The probability that a single pixel contains
  "only" background.}

  \item{sigma}{The Gaussian noise level for background-only
  points}

  \item{lambda}{Expected value (over background) of
  signal-containing pixels}

  \item{eta}{A limit to the straightness sensitivity}

  \item{Phi}{The spline matrix taking cc into the
  background curve

  Returns: The optimized values of cc.}
}
\description{
  Optimizes cc (i.e., the spline values; 'c' in Fischer et
  al.) with a given level of experimental noise.
}

