\name{PsiHSig}
\alias{PsiHSig}
\title{Signal-including contribution to likelihood}
\usage{
  PsiHSig(bgr, data, calc.grad, calc.hess)
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
  (named list) The contribution to the likelihood when a
  datapoint includes signal as well as background (and any
  requested derivatives with respect to the spline values):
  \item{funct}{(numeric) The function itself}
  \item{grad}{(numeric vector, or NA) The gradient of the
  function with respect to the spline values.}
  \item{hess}{(numeric matrix, or NA) The hessian of f with
  respect to the spline values}
}
\description{
  Signal-including contribution to likelihood
}
\seealso{
  \code{\link{PsiFBgr}}
}

