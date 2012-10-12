\name{PsiWrapped}
\alias{PsiWrapped}
\title{Posterior warpper for trust region algorithm}
\usage{
  PsiWrapped(params, bgr, data)
}
\arguments{
  \item{params}{(numeric vector) The initial parameters.}

  \item{bgr}{(\code{\link{background}} object) The
  estimated background (and associated hyperparameters)}

  \item{data}{(\code{\link{dataset}} object) The data whose
  background to estimate}
}
\value{
  As required by \code{trust()}, a list with three
  elements: \item{value}{(numeric) The value of the
  function} \item{gradient}{(numeric vector) The gradient
  of the optimization function w.r.t. the parameters}
  \item{hessian}{(numeric matrix) The matrix of the
  optimization function w.r.t. the parameters}
}
\description{
  Wrap the posterior function for use with
  \code{\link{trust}}.
}

