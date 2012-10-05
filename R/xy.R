# file:     R/xy.R
# author:   Charles R. Hogg III
# date:     2012-10-05
# purpose:  Generic methods for x and y coordinates

#' X-values for a function
x <- function(obj, ...) {
  UseMethod("x", obj)
}

#' Y-values for a function
y <- function(obj, ...) {
  UseMethod("y", obj)
}

