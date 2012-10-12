# file:    R/model.R
# author:  Charles R. Hogg III
# date:    2012-10-09
# purpose: Functions relating to the probability model in Fischer et al. (2000)

###############################################################################
# SECTION: Probability model functions

#' Background-only contribution to likelihood
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#' @param calc.grad  If TRUE, also calculate the gradient
#' @param calc.hess  If TRUE, also calculate the Hessian
#'
#' @seealso \code{\link{PsiHSig}}
#' @return (named list) The background-only contribution to the likelihood (and
#'     any requested derivatives with respect to the spline values):
#' \item{funct}{(numeric) The function itself}
#' \item{grad}{(numeric vector, or NA) The gradient of the function with
#'     respect to the spline values.}
#' \item{hess}{(numeric matrix, or NA) The hessian of f with respect to the
#'     spline values}
PsiFBgr <- function(bgr, data, calc.grad, calc.hess) {
  sigma <- Sigma(data)
  diff.abs <- Deviation(data, bgr)
  diff.norm <- diff.abs / sigma
  value <- (
    log(Beta(bgr))
    - 0.5 * log(2 * pi)
    - log(sigma)
    - 0.5 * (diff.norm ^ 2)
    )
  grad <- hess <- NA
  if (calc.grad) {  # Calculate the gradient
    grad <- as.vector(diff.abs / (sigma ^ 2)) * Phi(bgr, data)
  }
  if (calc.hess) {  # Calculate the Hessian
    PPT <- RowOuterProduct(Phi(bgr, data))
    sigma.rep <- array(rep(sigma, each=E(bgr) ^ 2), dim=dim(PPT))
    hess <- -PPT / (sigma.rep ^ 2)
  }
  return (list(funct=value, grad=grad, hess=hess))
}

#' Helper function for \code{\link{PsiHSig}}
#'
#' This combination of variables appears frequently in calculations for the
#' gradient and Hessian.
#'
#' @param qq (numeric vector) A convenient variable (see Charles Hogg's HTML
#'     notes for elaboration).
#' @param rho (numeric vector) The experimental noise at each datapoint, in
#'     units of the expected average signal magnitude.
GammaQ <- function(qq, rho) {
  gq <- exp(-0.5 * qq ^ 2 - pnorm(q=qq, log.p=TRUE))
  return (gq / (rho * sqrt(2 * pi)))
}

#' Signal-including contribution to likelihood
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#' @param calc.grad  If TRUE, also calculate the gradient
#' @param calc.hess  If TRUE, also calculate the Hessian
#'
#' @seealso \code{\link{PsiFBgr}}
#' @return (named list) The contribution to the likelihood when a datapoint
#'     includes signal as well as background (and any requested derivatives
#'     with respect to the spline values):
#' \item{funct}{(numeric) The function itself}
#' \item{grad}{(numeric vector, or NA) The gradient of the function with
#'     respect to the spline values.}
#' \item{hess}{(numeric matrix, or NA) The hessian of f with respect to the
#'     spline values}
PsiHSig <- function(bgr, data, calc.grad, calc.hess) {
  # Some of these variable names -- 'rho', 'z', 'qq' -- are based on my notes.
  # I broke the likelihood from Fischer et al. into pieces, and these variables
  # are for combinations which frequently occur.  The notes are in HTML format;
  # I left one copy with Igor Levin and one with Steve Lund.
  sigma <- Sigma(data)
  rho <- sigma / Lambda(bgr)
  diff.abs <- Deviation(data, bgr)
  z <- diff.abs / Lambda(bgr)
  qq <- z / rho - rho
  funct <- (
    log(1 - Beta(bgr))
    - log(Lambda(bgr))
    + pnorm(log.p=TRUE, q=qq)
    - z
    + 0.5 * (rho ^ 2)
    )
  grad <- hess <- NA
  if (calc.grad || calc.hess) {  # Calculate the gradient
    gamma.q <- GammaQ(qq=qq, rho=rho)
    grad <- as.vector((1 - gamma.q) / Lambda(bgr)) * Phi(bgr, data)
  }
  if (calc.hess) {               # Calculate the Hessian
    prefactor <- as.vector(-gamma.q * (gamma.q + qq / rho) / (Lambda(bgr) ^ 2))
    PPT <- RowOuterProduct(Phi(bgr, data))
    prefactor.rep <- array(rep(prefactor, each=E(bgr) ^ 2), dim=dim(PPT))
    hess <- prefactor.rep * PPT
  }
  return (list(funct=funct, grad=grad, hess=hess))
}

###############################################################################
# SECTION: Matrix bookkeeping functions
#
# The model calculations sometimes take us beyond 2D matrices.  For example, we
# need to calculate a Hessian (2 indices) for every datapoint (1 index): a
# total of 3 indices.  Simple-minded conventions for multiplying matrices and
# vectors no longer suffice; we have to be very careful.  This section holds
# the functions which handle these computations.

#' Outer product of matrix rows
#'
#' Computes the (3d) array consisting of the outer product of each row of m
#' with itself.
#'
#' @param m (numeric matrix) Matrix whose rows we desire to expand.  Usually a
#'     spline matrix, where the columns correspond to the knots, and the rows
#'     correspond to evaluation points).
#'
#' @return A 3-index numeric array, \code{value}, such that \code{value[, , i]}
#'     is the outer product of \code{m[i, ]} with itself.
RowOuterProduct <- function(m) {
  rows <- nrow(m)
  cols <- ncol(m)
  ppt <- array(apply(m, 1, function(x) x %*% t(x)), dim=c(cols, cols, rows))
  return (ppt)
}
