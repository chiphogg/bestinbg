\name{SplineMatrix}
\alias{SplineMatrix}
\title{Matrix of spline basis functions}
\usage{
  SplineMatrix(x, x.out, deriv = 0)
}
\arguments{
  \item{x}{(numeric vector) Spline knot positions}

  \item{x.out}{(numeric vector) x-values where we evaluate
  the basis functions.}

  \item{deriv}{(integer from 0 to 3) Which derivative of
  the spline function (0 means return the function itself)}
}
\value{
  (numeric matrix) Matrix whose columns are spline basis
  functions taking \code{x} into \code{x.out}.
}
\description{
  Computes the matrix whose columns are spline basis
  functions taking \code{x} into \code{x.out}.
}

