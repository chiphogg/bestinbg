\name{Psi}
\alias{Psi}
\title{Negative log of posterior}
\usage{
  Psi(bgr, data, calc.grad = TRUE, calc.hess = TRUE)
}
\arguments{
  \item{bgr}{(\code{\link{background}} object) The
  estimated background (and associated hyperparameters)}

  \item{data}{(\code{\link{dataset}} object) The data whose
  background to estimate}

  \item{calc.grad}{If TRUE, also calculate the gradient}

  \item{calc.hess}{If TRUE, also calculate the Hessian}
}
\value{
  (named list) The negative logarithm of the posterior (and
  any requested derivatives with respect to the spline
  values): \item{funct}{(numeric) The function itself}
  \item{grad}{(numeric vector, or NA) The gradient of the
  function with respect to the spline values.}
  \item{hess}{(numeric matrix, or NA) The hessian of Psi
  with respect to the spline values}
}
\description{
  Computes the optimization function from Fischer et al.:
  the negative log of the posterior.
}

