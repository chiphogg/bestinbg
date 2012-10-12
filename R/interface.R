# file:    R/interface.R
# author:  Charles R. Hogg III
# date:    2012-10-05
# purpose: Functions which end users of the package are likely to call.

###############################################################################
# SECTION: Probability model functions

#' Negative log of posterior
#'
#' Computes the optimization function from Fischer et al.: the negative log of
#' the posterior.
#'
#' @param bgr (\code{\link{background}} object) The estimated background (and
#'     associated hyperparameters)
#' @param data (\code{\link{dataset}} object) The data whose background to
#'     estimate
#' @param calc.grad  If TRUE, also calculate the gradient
#' @param calc.hess  If TRUE, also calculate the Hessian
#'
#' @export
#' @return (named list) The negative logarithm of the posterior (and any
#' requested derivatives with respect to the spline values):
#' \item{funct}{(numeric) The function itself}
#' \item{grad}{(numeric vector, or NA) The gradient of the function with
#'     respect to the spline values.}
#' \item{hess}{(numeric matrix, or NA) The hessian of Psi with respect to the
#'     spline values}
Psi <- function(bgr, data, calc.grad=TRUE, calc.hess=TRUE) {
  # Some notes on notation:
  # QUANTITY              | FISCHER | US
  # ----------------------+---------+-------------
  # No. of spline knots   | E       | E
  # Spline knot locations | xi      | X(bgr)
  # Spline knot values    | c       | Y(bgr)

  E <- length(X(bgr))
  d.mat <- DMatrix(x=X(bgr))
  cDc <- as.vector(t(Y(bgr)) %*% d.mat %*% Y(bgr)) + Eta(bgr)
  # Prior: separate terms according to which variables contribute
  contr.E <- (0.5 * E * log(pi) - log(gamma(0.5 * E)) 
    - sum(log(1:(E - 2))))  # "Depends" on xi, but constant as long as xi OK.
  contr.E.xi <- (-0.5 * log(DetD(d.mat)))
  contr.E.xi.cc <- 0.5 * E * log(cDc)
  psi.prior <- contr.E + contr.E.xi + contr.E.xi.cc
  # Likelihood:
  # (Includes my idiosyncratic notation:
  #   'f' is for points containing only background;
  #   'h' is for points which have signal, too;
  #   'g' combines both 'f' and 'h'.)
  calc.grad <- calc.grad || calc.hess  # Hessian requires gradient!
  f <- PsiFBgr(bgr=bgr, data=data, calc.grad=calc.grad, calc.hess=calc.hess)
  h <- PsiHSig(bgr=bgr, data=data, calc.grad=calc.grad, calc.hess=calc.hess)
  # log(exp(a) + exp(b)) has some SERIOUS pitfalls for the unwary...
  max.log <- pmax(f$funct, h$funct)
  f.fract <- as.vector(1 / (1 + exp(h$funct - f$funct)))
  g <- max.log + log(exp(f$funct - max.log) + exp(h$funct - max.log))
  psi.likelihood <- -sum(g)
  # Calculate higher derivatives, if requested
  grad <- hess <- NA
  if (calc.grad) {  # Calculate the gradient
    grad.g <- f.fract * f$grad + (1 - f.fract) * h$grad
    grad <- -colSums(grad.g) + (t(Y(bgr)) %*% d.mat) * E / cDc
  }
  if (calc.hess) {  # Calculate the Hessian
    f.contr <- (    f.fract) * (f$hess + RowOuterProduct(f$grad))
    h.contr <- (1 - f.fract) * (h$hess + RowOuterProduct(h$grad))
    hess.g <- f.contr + h.contr - RowOuterProduct(grad.g)
    hess <- apply(hess.g, c(1, 2), sum) + E * (
      (d.mat / cDc) - (2 * d.mat %*% Y(bgr) %*% t(Y(bgr)) %*% d.mat) / (cDc ^ 2)
      )
  }
  return (list(funct=psi.prior + psi.likelihood, grad=grad, hess=hess))
}

###############################################################################
# SECTION: Optimization functions

#' Compute background values at spline knots
#'
#' @param bgr (background object) The current background estimate (includes all
#'     relevant hyperparameters).
#' @param data (dataset object) The data whose background to estimate.
#'
#' @export
#' @return (background object) The most probable spline knot values
OptimizeKnotValues <- function(bgr, data) {
  # We follow the strategy outlined in Fischer et al. (2000):
  #   1) Fit using artificially inflated noise (comparable to the mean *signal*
  #      level, "lambda")
  #   2) Fit with progressively smaller noise (down to the true value), always
  #      using the previous fit as a starting point.
  # This helps avoid getting stuck in local optima.

  # We start doubled, because we halve *before* optimizing.  This lets me use a
  # "while" loop rather than "repeat/break", which I find less readable.
  sigma <- 2.0 * Lambda(bgr)

  sigma.min <- min(Sigma(data))
  while (sigma > sigma.min) {
    sigma <- sigma * 0.5
    bgr <- OptimizeKnotValuesInflatedNoise(bgr, data, sigma)
  }
}

#' Random function with smooth background
#'
#' Generates a random function with a spline background, random peaks, and
#' Gaussian noise.
#'
#' @param x (numeric vector) The x-points where data should be generated
#' @param lambda (numeric; positive) The mean signal magnitude (approximate!)
#' @param sigma (numeric; positive) The noise level to generate
#' @param delta.x (numeric; positive) The minimum spacing allowed between
#'     spline knots
#' @param E (numeric; positive integer) The number of spline knots
#' @param n.peaks (numeric; positive integer) The number of peaks to generate
#'
#' There are some obvious constraints on these parameters.  For instance,
#' \code{(E - 1) * delta.x < diff(range(x))} should always be satisfied.
#'
#' @export
NoisyTestFunction <- function(x, lambda, sigma, delta.x, E, n.peaks) {
  # 1.  Generate the signal (random peaks).
  # Assume mean height roughly half peak height (hence "2/lambda"):
  amplitude <- rexp(n=n.peaks, rate=(2 / lambda))
  # Peaks should fit (within +-2 sigma) than any 2 spline knots
  sigma.range <- c(0.4, 1.0) * delta.x / 4
  # Keep peak centres at least 2 sigma away from the edge:
  width <- runif(n=n.peaks, min=sigma.range[1], max=sigma.range[2])
  centre <- runif(n=n.peaks,
    min=min(x) + (2 * width),
    max=max(x) - (2 * width))
  separate.peaks <- mapply(
    FUN=function(amp, wid, cen) amp * exp (-0.5 * ((x - cen) / wid) ^ 2),
    amp=amplitude, wid=width, cen=centre)
  peaks <- rowSums(separate.peaks)

  # 2. Add some noise
  noise <- rnorm(sd=sigma, n=length(x))

  # 3. Construct dataset and background objects
  data <- new_dataset(x=x, y=peaks + noise, sigma=sigma)
  bgr <- new_background(delta.x=delta.x, lambda=lambda, beta=0.5, 
    epsilon=sigma * 0.001, E=E, data=data)

  # 4. Construct random background and add to data
  Y(bgr) <- rnorm(n=E(bgr))
  XFractions(bgr) <- NULL  # Choose knot positions randomly
  Y(data) <- Y(data) + Y(bgr, X(data))

  return (list(data=data, bgr=bgr))
}

