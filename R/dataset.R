# file:     R/dataset.R
# author:   Charles R. Hogg III
# date:     2012-10-05
# purpose:  'dataset' class -- A set of (x, y) pairs, and a noise level (which
#   can be different for each datapoint).

###############################################################################
# dataset: constructors

#' Data whose background we want to estimate
#'
#' @export
new_dataset <- function(x, ...) {
  UseMethod("new_dataset", x)
}

#' New dataset from numeric data
#'
#' @S3method new_dataset numeric
#'
#' @param x (numeric vector) X-values where we have data
#' @param y (numeric vector) Measured datapoints
#' @param sigma (numeric) Gaussian uncertainty (optionally different for each
#'     datapoint)
new_dataset.numeric <- function(x, y, sigma) {
  d <- data.frame(x=x, y=y, sigma=sigma)
  class(d) <- "dataset"
  return(d)
}

#' New dataset from datafile
#'
#' @S3method new_dataset character
#'
#' @param f.name (character) The name of a datafile which has 2 columns (3 if
#'     sigma is included)
#' @param sigma (numeric) Gaussian uncertainty (optionally different for each
#'     datapoint)
new_dataset.character <- function(f.name, sigma=NA) {
  d <- data.frame(x=1, y=2, sigma=3)
  class(d) <- "dataset"
  return (d)
}

###############################################################################
# dataset: other functions

#' Is this object a dataset?
#'
#' @param x Object to test
#'
#' @export
#' @return \code{TRUE} if this is a \code{dataset} object.
is.dataset <- function(x) inherits(x, "dataset")

#' @S3method x dataset
x.dataset <- function(d) {
  return (d$x)
}

