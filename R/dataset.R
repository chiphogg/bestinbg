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
  class(d) <- c("dataset", "data.frame")
  return(d)
}

#' New dataset from datafile
#'
#' Create a new \code{dataset} object by reading from a datafile.
#'
#' This function expects at least two columns (three if \code{sigma} is not
#' passed).  \code{columns} is a numeric vector of length 2 (or 3) indicating
#' which datafile columns contain the data.
#' \describe{
#'   \item{\code{columns[1]}}{The x-values}
#'   \item{\code{columns[2]}}{The y-values}
#'   \item{\code{columns[3]}}{The sigma-values (unless \code{sigma} is supplied)}
#' }
#'
#' @S3method new_dataset character
#'
#' @param x (character) The name of a datafile which has 2 columns (3 if
#'     sigma is included)
#' @param sigma (numeric) Gaussian uncertainty (optionally different for each
#'     datapoint)
#' @param columns (numeric vector) The column numbers corresponding to x, y,
#'     and sigma, respectively (the column for sigma may be omitted if
#'     \code{sigma} is supplied as an argument).
#' @param ... Other parameters for \code{\link{read.table}}
new_dataset.character <- function(x, sigma=NA, columns=1:3, ...) {
  d <- read.table(file=x, header=TRUE, ...)
  colnames(d)[columns[1:2]] <- c("x", "y")
  if (is.na(sigma)) {
    colnames(d)[columns[3]] <- "sigma"
  } else {
    d$sigma <- sigma
  }
  class(d) <- c("dataset", "data.frame")
  return (d)
}

###############################################################################
# dataset: other functions

#' Is this object a dataset?
#'
#' @param x Object to test
#' TODO
#' @export
#' @return \code{TRUE} if this is a \code{dataset} object.
is.dataset <- function(x) inherits(x, "dataset")

#' @S3method X dataset
X.dataset <- function(d) {
  return (d$x)
}

#' Y-values for a dataset
#'
#' @S3method Y dataset
Y.dataset <- function(d) {
  return (d$y)
}

#' @rdname Y.dataset
#' @S3method Y<- dataset
#' @param value (numeric vector) New Y-values for the dataset
#' @export
`Y<-.dataset` <- function(d, value) {
  # This idiom prevents changing the length of Y
  d$y[1:length(d$y)] <- value
  return (invisible(d))
}

#' Datapoint deviation from background
#'
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#'
#' @return (numeric vector) The difference between the measured data and the
#'     estimated background
#' @export
Deviation <- function(data, bgr) {
  return (Y(data) - Y(bgr, x=X(data)))
}

#' Pointwise experimental uncertainty
#'
#' @param data (\code{\link{dataset}} object) The measured data
#'
#' @return (numeric vector) The pointwise experimental uncertainty (assuming
#'     Gaussian noise)
#' @export
Sigma <- function(data) {
  return (data$sigma)
}

#' @rdname Sigma
#' @export
`Sigma<-` <- function(data, value) {
  # This idiom prevents changing the length of Sigma
  data$sigma[1:length(data$sigma)] <- value
  return (data)
}

#' Plot a dataset object
#'
#' @export
plot.dataset <- function(data, ...) {
  plot(data$x, data$y, ...)
}
