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
#' @param delta.x (numeric) The minimum allowed separation distance between two
#'     spline knots (larger \code{delta.x} yields smoother backgrounds)
#' @param lambda (numeric) Estimated average signal value (above the
#'     background)
#' @param beta (numeric) Probability for a point to contain only background
#' @param eta (numeric) Maximum straightness threshhold; should be computed
#'     using \code{\link{Eta}}
#' @param E (numeric) Number of spline knots to use
#' @param data (\code{\link{dataset}}) The dataset to analyze
new_background <- function(delta.x, lambda, beta, E, data) {
  x.range <- range(x(data))
  knots.x <- seq(from=x.range[1], to=x.range[2], length.out=E)
  knots.y <- rep(mean(y(data)), E)
  b <- list(
    hypers=list(
      delta.x=delta.x,
      lambda=lambda,
      beta=beta,
      eta=eta),
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



