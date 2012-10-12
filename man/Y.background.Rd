\name{Y.background}
\alias{Y.background}
\alias{Y<-.background}
\title{Estimate of background function}
\usage{
  Y.background(b, x = NA)

  Y.background(b, value) <- value
}
\arguments{
  \item{x}{(optional; numeric vector) X-values where the
  background estimate should be evaluated.}

  \item{value}{(numeric vector) New Y-values for the
  background estimate}
}
\description{
  By default, the spline knot values.  However, if the user
  supplies a numeric vector of \code{x}-values, the spline
  will be evaluated at those points.
}

