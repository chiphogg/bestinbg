# file:    R/optimizers.R
# author:  Charles R. Hogg III
# date:    2012-10-05
# purpose: Optimization code which shouldn't be exposed to users.

require("trust")

#' Posterior warpper for trust region algorithm
#'
#' Wrap the posterior function for use with \code{\link{trust}}.
#'
#' @param params (numeric vector) The initial parameters.
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#'
#' @return As required by \code{trust()}, a list with three elements:
#'  \item{value}{(numeric) The value of the function}
#'  \item{gradient}{(numeric vector) The gradient of the optimization function
#'      w.r.t. the parameters}
#'  \item{hessian}{(numeric matrix) The matrix of the optimization function
#'      w.r.t. the parameters}
PsiWrapped <- function(params, bgr, data) {
  p <- Psi(bgr=bgr, data=data, calc.grad=TRUE, calc.hess=TRUE)
  return (list(value=p$funct, gradient=as.vector(p$grad), hessian=p$hess))
}

#' Optimize knot values with inflated noise level.
#'
#' We artificially inflate the noise level to avoid getting stuck in local
#' optima.
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#'
#' @return (\code{\link{background}} object) A copy of \code{bgr}, but with
#'     \code{Y(bgr)} optimized.
OptimizeKnotValuesInflatedNoise <- function(bgr, data, sigma) {
  Sigma(data) <- pmax(Sigma(data), sigma)
  opt <- trust(objfun=PsiWrapped, parinit=Y(bgr),
    rinit=1e-2, rmax=1e6, blather=FALSE,
    bgr=bgr, data=data)
  Y(bgr) <- opt$argument
  return (bgr)
}

#' Optimizes cc (i.e., the spline values; 'c' in Fischer et al.) with a given
#' level of experimental noise.
#'
#' @param x  The x-values of the measured datapoints.
#' @param y  The y-values of the measured datapoints.
#' @param cc  Initial guess for the spline knot y-values ('c' in Fischer et al.).
#' @param xi  The spline knot x-values.
#' @param p.bgr  The probability that a single pixel contains "only" background.
#' @param sigma  The Gaussian noise level for background-only points
#' @param lambda  Expected value (over background) of signal-containing pixels
#' @param eta  A limit to the straightness sensitivity
#' @param Phi  The spline matrix taking cc into the background curve
#'
#' Returns:
#'   The optimized values of cc.
OptimizeCC <- function(x, y, cc=(0 * xi), xi, p.bgr, sigma, lambda, eta,
  Phi=spline.matrix(x=x, spline.x=xi)) {
  opt <- trust(objfun=PsiWrapped, parinit=cc, rinit=1e-2, rmax=1e6,
    blather=FALSE,
    x=x, y=y, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, Phi=Phi)
  return (opt$argument)
}

