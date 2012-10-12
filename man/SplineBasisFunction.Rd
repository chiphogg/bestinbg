\name{SplineBasisFunction}
\alias{SplineBasisFunction}
\title{Spline basis functions.}
\usage{
  SplineBasisFunction(x, x.out, i, deriv = 0)
}
\arguments{
  \item{x}{(numeric vector) Spline knot positions}

  \item{x.out}{(numeric vector) x-values where we evaluate
  the basis functions.}

  \item{i}{(integer from 1 to \code{length(x.out)}): which
  basis function to return}

  \item{deriv}{(integer from 0 to 3) Which derivative of
  the spline function (0 means return the function itself)}
}
\value{
  (numeric vector of length \code{length(x.out)}) The
  requested derivative of the requested basis function,
  evaluated at the points in \code{x.out}.
}
\description{
  Compute the i'th spline basis function for taking
  \code{x} into \code{x.out}.
}

