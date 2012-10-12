# file:     R/background.R
# author:   Charles R. Hogg III
# date:     2012-10-05
# purpose:  'background' class -- a set of spline knots defining a background
#    estimate, together with hyperparameters encoding underlying assumptions.

#' Background estimate
#'
#' A set of spline knots defining a background estimate, together with
#' hyperparameters encoding underlying assumptions.
#'
#' If \code{epsilon} is zero, there is a \dQuote{straight line trap}: the model
#' assigns infinite probability to perfectly straight lines.  A small but
#' nonzero \code{epsilon} circumvents this: the model won't discriminate
#' among straight lines which vary (vertically) by this amount or less.
#'
#' @rdname background
#' @aliases background new_background
#'
#' @param delta.x (numeric) The minimum allowed separation distance between two
#'     spline knots (larger \code{delta.x} yields smoother backgrounds)
#' @param lambda (numeric) Estimated average signal value (above the
#'     background)
#' @param beta (numeric) Probability for a point to contain only background
#' @param epsilon (numeric) \dQuote{Vertical} precision: we don't worry if
#'     spline knots fluctuate by this amount or less. 
#' @param E (numeric) Number of spline knots to use
#' @param data (\code{\link{dataset}}) The dataset to analyze
#' @export
new_background <- function(delta.x, lambda, beta, epsilon, E, data) {
  x.range <- range(X(data))
  knots.x <- seq(from=x.range[1], to=x.range[2], length.out=E)
  knots.y <- rep(mean(Y(data)), E)
  b <- list(
    hypers=list(
      delta.x=delta.x,
      lambda=lambda,
      beta=beta,
      epsilon=epsilon,
      eta.universal=EtaUniversal(diff(range(X(data))) / delta.x),
      x.bounds=x.range),
    knots=data.frame(
      x=knots.x,
      y=knots.y))
  class(b) <- "background"
  return (b)
}

#' Is this object a background?
#'
#' @param x Object to test
#'
#' @export
#' @return \code{TRUE} if this is a \code{background} object.
is.background <- function(x) inherits(x, "background")

###############################################################################
# background: other functions

DeltaX <- function(b) {
  return (b$hypers$delta.x)
}

FreeXSpace <- function(b) {
  return (diff(b$hypers$x.bounds) - (E(b) - 1) * DeltaX(b))
}

#' Spline knot positions as fractions of available space
#'
#' @param b (\code{\link{background}} object) An estimated background (and
#'     associated hyperparameters)
#'
#' @export
XFractions <- function(b) {
  return ((diff(X(b)) - DeltaX(b)) / FreeXSpace(b))
}

#' @rdname XFractions
#' @param value (numeric vector) New fractions to 
`XFractions<-` <- function(b, value) {
  new.fracts <- NA
  if (is.vector(value)) {
    if (length(value) < E(b) - 2 || length(value) > E(b) - 1) {
      warn("'value' has the wrong length")
    } else if (min(value) < 0 || max(value) > 1) {
      warn("'value' can only contain values between 0 and 1")
    } else {
      new.fracts <- value
    }
  } else if (is.null(value)) {
    new.fracts <- sort(runif(n=E(b) - 2, min=0, max=1))
  }
  # At this point, if new.fracts is not NA, it's probably legit.
  if (!identical(new.fracts, NA)) {
    # First knot is always 0; last is always 1; fractions are always sorted.
    fracts.processed <- c(0, sort(new.fracts)[1:(E(b) - 2)], 1)
    distances <- diff(fracts.processed * FreeXSpace(b)) + DeltaX(b)
    b$knots$x <- b$knots$x[1] + cumsum(c(0, distances))
  }
  return (invisible(b))
}

#' Spline knot positions for background function
#'
#' @S3method X background
#' @export
X.background <- function(b) {
  return (b$knots$x)
}

#' @rdname X.background
#' @param value (numeric vector) New X-values for the background spline knots
#' @export
`X<-.background` <- function(b, value) {
  # This idiom prevents changing the length of X
  b$knots$x[1:length(b$knots$x)] <- value
  return (invisible(b))
}

#' Estimate of background function
#'
#' By default, the spline knot values.  However, if the user supplies a numeric
#' vector of \code{x}-values, the spline will be evaluated at those points.
#'
#' @param x (optional; numeric vector) X-values where the background estimate
#'     should be evaluated.
#'
#' @S3method Y background
#' @export
Y.background <- function(b, x=NA) {
  if (identical(x, NA)) {
    return (b$knots$y)
  }
  return (spline(x=b$knots$x, y=b$knots$y, xout=x)$y)
}

#' @rdname Y.background
#' @S3method Y<- background
#' @param value (numeric vector) New Y-values for the background estimate
#' @export
`Y<-.background` <- function(b, value) {
  # This idiom prevents changing the length of Y
  b$knots$y[1:length(b$knots$y)] <- value
  return (invisible(b))
}

#' Spline matrix
#'
#' The matrix for estimating background values given spline values.
#'
#' Referred to as \code{Phi} in Fischer et al. (1999); hence the function name.
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#'
#' @return Matrix which takes spline knot values into background estimates at
#'     the datapoints
Phi <- function(bgr, data) {
  return (SplineMatrix(x=X(bgr), x.out=X(data)))
}

#' Number of spline knots
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#'
#' @return The number of spline knots in \code{bgr}
E <- function(bgr) {
  return (length(X(bgr)))
}

#' Mean signal magnitude.
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' 
Lambda <- function(bgr) {
  return (bgr$hypers$lambda)
}

#' Probability for pixel to contain only background
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' 
Beta <- function(bgr) {
  return (bgr$hypers$beta)
}

#' Maximum straightness
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' 
Eta <- function(bgr) {
  h <- bgr$hypers
  return ((h$epsilon ^ 2) * h$eta.universal / (h$delta.x ^ 3))
}

#' Add a background object to a plot
#'
#' @param x \code{\link{background}} object to plot
#' @param x.vals (numeric vector) x-values to evaluate the background function
#'
#' @export
points.background <- function(x, x.vals=X(x), ...) {
  points(x.vals, Y(x, x.vals), ...)
}
