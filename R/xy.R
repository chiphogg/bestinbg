# file:     R/xy.R
# author:   Charles R. Hogg III
# date:     2012-10-05
# purpose:  Generic methods for x and y coordinates

#' X-values for a function
#'
#' @export
X <- function(obj, ...) {
  UseMethod("X", obj)
}
#' @export
`X<-` <- function(obj, ...) {
  UseMethod("X<-", obj)
}

#' Y-values for a function
#'
#' @export
Y <- function(obj, ...) {
  UseMethod("Y", obj)
}
#' @export
`Y<-` <- function(obj, ...) {
  UseMethod("Y<-", obj)
}

