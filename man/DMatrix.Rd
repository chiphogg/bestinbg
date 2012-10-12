\name{DMatrix}
\alias{DMatrix}
\title{Spline second-derivative overlap matrix}
\usage{
  DMatrix(x, only.trace = FALSE, robust.factor = 1e-12)
}
\arguments{
  \item{x}{Numeric vector of spline knot positions}

  \item{only.trace}{If true, return trace(D) (a single
  number) instead of D.}

  \item{robust.factor}{We add a constant to all eigenvalues
  to promote stability; this is the ratio of that constant
  to the smallest eigenvalue which is SUPPOSED to be
  nonzero.}
}
\value{
  A matrix M where M[i, j] gives the integral of the
  product of the second derivatives of spline basis
  functions i and j.
}
\description{
  Computes the second derivative overlap matrix ('D' in
  Fischer et al. (1999)) of the spline basis functions with
  knots at \code{x}.
}
\details{
  The second derivative of a natural cubic spline is a
  piecewise linear function, whose kinks coincide with the
  spline knots.  This is very convenient for an analytical
  treatment.
}

