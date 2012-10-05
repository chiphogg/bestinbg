##############################################################################
# SECTION: Spline Knot Representations
#
# Spline knots can be moved left or right, but they can't approach closer than
# some minimum distance ("min.dx").  There are two ways to represent spline
# knot positions (assume there are 'E' knots; the first is fixed at 'a' and the
# last at 'b'):
#
#   1) Natural Form
#       A sorted vector of 'E' numeric values which starts at 'a' and ends at
#       'b'.  This is easy to interpret, but probabilistic calculations and
#       proposals become very complicated.
#
#   2) Fractions of Available Space
#       A sorted vector of ('E' - 1) numeric values 
##############################################################################


#' Available space for moving spline knots
#'
#' Spline knots cannot be closer than a minimum spacing \code{min.dx}.  Thus,
#' the total \emph{available} range (the \dQuote{wiggle room}) is less than
#' \code{x.range}.
#'
#' @param E (positive integer) The number of spline knots.
#' @param min.dx (positive numeric) The minimum spacing of the spline knots.
#' @param x.range (positive numeric) The width of the region where we're trying
#'     to estimate the background.
#' 
#' @seealso \code{\link{WiggleRoomFractions}}
#'
#' @return The total amount of free space for moving the spline knots.
WiggleRoom <- function(E, min.dx, x.range) {
  x.range - (E - 1) * min.dx
}

#' Normalized distances between consecutive spline knots
#'
#' @param xi (sorted numeric vector) The x-coordinates of the spline knots.
#' @param min.dx (positive numeric) The minimum spacing of the spline knots.
#'
#' @seealso \code{\link{WiggleRoom}}
#' @seealso \code{\link{Xi}}
#'
#' @return A numeric vector whose entries are non-negative and sum to 1.
WiggleRoomFractions <- function(xi, min.dx) {
  shares <- diff(xi) - min.dx
  fractions <- shares / sum(shares)
  return (fractions)
}

#' Compute spline knot positions
#'
#' Spline knot positions for background functions must have several properties:
#' the endpoints are fixed, they can't be closer than some minimum distance,
#' etc.  It's more convenient to store them in a format which *automatically*
#' fulfills these requirements: here, we use the fraction of the *actually
#' available* space (i.e., the space leftover after all the "minimum distance"
#' buffers are taken into account).  This function goes from this abstract
#' representation back to actual positions.
#'
#' @param fracts (numeric vector) Fractional distances between consecutive
#'     spline-knots.  Non-negative; sums to 1.
#' @param min.dx (positive numeric) The minimum spacing of the spline knots.
#' @param xi.min (numeric) The (fixed) position of the leftmost knot.
#' @param xi.max (numeric) The (fixed) position of the rightmost knot.
#'
#' @seealso \code{\link{WiggleRoomFractions}}
Xi <- function(fracts, min.dx, xi.min, xi.max) {
  x.range <- xi.max - xi.min
  room <- WiggleRoom(E=(length(fracts) + 1), min.dx=min.dx, x.range=x.range)
  spaces <- fracts * room
  return (xi.min + cumsum(c(0, (spaces + min.dx))))
}

##############################################################################
# SECTION: "Core" spline functions

#' Spline basis functions.
#'
#' Compute the i'th spline basis function for taking spline.x into x.
#'
#' @param x  Numeric vector giving x-values where we evaluate the basis functions.
#' @param spline.x  Numeric vector of spline knot positions
#' @param i  Integer from 1 to length(spline.x): which basis function to return?
#' @param deriv  Which derivative of the spline function? (0 to 3; 0 is
#'     function itself)
#'
#' @return A numeric vector of length length(x): the requested derivative of
#'     the requested basis function, evaluated at the points in x.
spline.basis.function <- function(x, spline.x, i, deriv=0) {
  y <- 0 * spline.x
  y[i] <- 1
  spline.fun <- splinefun(x=spline.x, y=y, method="natural")
  return (spline.fun(x=x, deriv=deriv))
}

#' Matrix of spline basis functions
#'
#' Computes the matrix whose columns are spline basis functions taking spline.x
#' into x.
#'
#' @param x  Numeric vector giving x-values where we evaluate the basis functions.
#' @param spline.x  Numeric vector of spline knot positions
#' @param deriv  Which derivative of the spline function? (0 to 3; 0 is
#'     function itself)
#'
#' @return The matrix whose columns are spline basis functions taking spline.x
#'     into x.
spline.matrix <- function(x, spline.x, deriv=0) {
  N <- length(spline.x)
  N.out <- length(x)
  M <- matrix(NA, nrow=N.out, ncol=N)
  for (i in 1:N) {
    M[, i] <- spline.basis.function(x=x, spline.x=spline.x, i=i, deriv=deriv)
  }
  return (M)
}

#' Computes the second derivative overlap matrix ('D' in Fischer et al.
#' (1999)) of the spline basis functions with knots at 'spline.x'.
#'
#' @param spline.x  Numeric vector of spline knot positions
#' @param only.trace  If true, return trace(D) (a single number) instead of D.
#' @param robust.factor  We add a constant to all eigenvalues to promote
#'     stability; this is the ratio of that constant to the smallest eigenvalue
#'     which is SUPPOSED to be nonzero.
#'
#' @return A matrix DD where DD(i, j) gives the integral of the product of the
#'     second derivatives of spline basis functions i and j.
D.Matrix <- function(spline.x, only.trace=FALSE, robust.factor=1e-12) {
  N <- length(spline.x)
  ddM <- spline.matrix(x=spline.x, spline.x=spline.x, deriv=2)
  x.mat <- matrix(rep(spline.x, N), nrow=N)
  dy <- apply(ddM, 2, diff)
  dx <- apply(x.mat, 2, diff)
  # slope[i, j] is the slope of the i'th segment of the j'th basis function;
  # intercept[i, j] is its y-intercept
  slope <- dy / dx
  intercept <- ddM[-1, ] - slope * x.mat[-1, ]
  d.x3 <- diff(spline.x ^ 3)
  d.x2 <- diff(spline.x ^ 2)
  # The for-loop code below reduces to this if we only care about diagonal
  # elements (which is the case for computing the trace).
  if (only.trace) {
    return (sum((slope ^ 2) * d.x3 / 3.0 + (slope * intercept) * d.x2
        + intercept ^ 2 * diff(spline.x)))
  }
  # If we're this far, we need to compute the whole matrix.
  DD <- matrix(NA, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in i:N) {
      prod.ss <- slope[, i] * slope[, j] / 3.0
      prod.si <- 0.5 * (
        slope[, i] * intercept[, j] + slope[, j] * intercept[, i])
      prod.ii <- intercept[, i] * intercept[, j]
      integrals <- prod.ss * d.x3 + prod.si * d.x2 + prod.ii * diff(spline.x)
      DD[i, j] <- DD[j, i] <- sum(integrals)
    }
  }
  # Make sure we avoid negative eigenvalues!
  eigenvals <- eigen(x=DD, symmetric=TRUE, only.values=TRUE)$values
  lowest.true.eigenval <- sort(eigenvals)[3]
  return (DD + lowest.true.eigenval * robust.factor * diag(nrow(DD)))
}


