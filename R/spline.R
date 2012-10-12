# file:     R/spline.R
# author:   Charles R. Hogg III
# date:     2012-10-09
# purpose:  Helper functions for splines

#' Spline basis functions.
#'
#' Compute the i'th spline basis function for taking \code{x} into \code{x.out}.
#'
#' @param x  (numeric vector) Spline knot positions
#' @param x.out  (numeric vector) x-values where we evaluate the basis functions.
#' @param i  (integer from 1 to \code{length(x.out)}): which basis function to
#'     return
#' @param deriv (integer from 0 to 3) Which derivative of the spline function
#'     (0 means return the function itself)
#'
#' @return (numeric vector of length \code{length(x.out)}) The requested
#'     derivative of the requested basis function, evaluated at the points in
#'     \code{x.out}.
SplineBasisFunction <- function(x, x.out, i, deriv=0) {
  y <- 0 * x
  y[i] <- 1
  spline.fun <- splinefun(x=x, y=y, method="natural")
  return (spline.fun(x=x.out, deriv=deriv))
}

#' Matrix of spline basis functions
#'
#' Computes the matrix whose columns are spline basis functions taking \code{x}
#' into \code{x.out}.
#'
#' @param x  (numeric vector) Spline knot positions
#' @param x.out  (numeric vector) x-values where we evaluate the basis functions.
#' @param deriv (integer from 0 to 3) Which derivative of the spline function
#'     (0 means return the function itself)
#'
#' @return (numeric matrix) Matrix whose columns are spline basis functions
#'     taking \code{x} into \code{x.out}.
SplineMatrix <- function(x, x.out, deriv=0) {
  n <- length(x)
  n.out <- length(x.out)
  M <- matrix(NA, nrow=n.out, ncol=n)
  for (i in 1:n) {
    M[, i] <- SplineBasisFunction(x=x, x.out=x.out, i=i, deriv=deriv)
  }
  return (M)
}

#' Spline second-derivative overlap matrix
#'
#' Computes the second derivative overlap matrix ('D' in Fischer et al.
#' (1999)) of the spline basis functions with knots at \code{x}.
#'
#' The second derivative of a natural cubic spline is a piecewise linear
#' function, whose kinks coincide with the spline knots.  This is very
#' convenient for an analytical treatment.
#'
#' @param x  Numeric vector of spline knot positions
#' @param only.trace  If true, return trace(D) (a single number) instead of D.
#' @param robust.factor  We add a constant to all eigenvalues to promote
#'     stability; this is the ratio of that constant to the smallest eigenvalue
#'     which is SUPPOSED to be nonzero.
#'
#' @return A matrix M where M[i, j] gives the integral of the product of the
#'     second derivatives of spline basis functions i and j.
DMatrix <- function(x, only.trace=FALSE, robust.factor=1e-12) {
  n <- length(x)
  d.d.M <- SplineMatrix(x=x, x.out=x, deriv=2)
  x.mat <- matrix(rep(x, n), nrow=n)
  d.y <- apply(d.d.M, 2, diff)
  d.x <- apply(x.mat, 2, diff)
  # slope[i, j] is the slope of the i'th segment of the j'th basis function;
  # intercept[i, j] is its y-intercept
  slope <- d.y / d.x
  intercept <- d.d.M[-1, ] - slope * x.mat[-1, ]
  d.x3 <- diff(x ^ 3)
  d.x2 <- diff(x ^ 2)
  # Computation is simpler if we only care about the trace.
  if (only.trace) {
    return (sum((slope ^ 2) * d.x3 / 3.0 + (slope * intercept) * d.x2
        + intercept ^ 2 * diff(x)))
  }
  # If we're this far, we need to compute the whole matrix.
  DD <- matrix(NA, nrow=n, ncol=n)
  for (i in 1:n) {
    for (j in i:n) {
      prod.ss <- slope[, i] * slope[, j] / 3.0
      prod.si <- 0.5 * (
        slope[, i] * intercept[, j] + slope[, j] * intercept[, i])
      prod.ii <- intercept[, i] * intercept[, j]
      integrals <- prod.ss * d.x3 + prod.si * d.x2 + prod.ii * diff(x)
      DD[i, j] <- DD[j, i] <- sum(integrals)
    }
  }
  # Make sure we avoid negative eigenvalues!
  eigenvals <- eigen(x=DD, symmetric=TRUE, only.values=TRUE)$values
  lowest.true.eigenval <- sort(eigenvals)[3]
  return (DD + lowest.true.eigenval * robust.factor * diag(nrow(DD)))
}

#' Computes the *modified* determinant (i.e., product of highest (N-2)
#' eigenvalues) for a D-matrix (see Fischer et al. for an explanation of why we
#' would want this).
#'
#' @param m  (numeric matrix) Overlap integrals for second derivatives of
#'      spline basis functions
#'
#' @return The product of the largest all-but-two eigenvalues of m (the
#'     smallest two should be identically zero, although they won't due to
#'     roundoff errors).
DetD <- function(m) {
  eigenvals <- eigen(x=m, symmetric=TRUE, only.values=TRUE)$values
  return (prod(rev(eigenvals)[-(1:2)]))
}

