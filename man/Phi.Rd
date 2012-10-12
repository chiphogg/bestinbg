\name{Phi}
\alias{Phi}
\title{Spline matrix}
\usage{
  Phi(bgr, data)
}
\arguments{
  \item{bgr}{(\code{\link{background}} object) The
  estimated background (and associated hyperparameters)}

  \item{data}{(\code{\link{dataset}} object) The data whose
  background to estimate}
}
\value{
  Matrix which takes spline knot values into background
  estimates at the datapoints
}
\description{
  The matrix for estimating background values given spline
  values.
}
\details{
  Referred to as \code{Phi} in Fischer et al. (1999); hence
  the function name.
}

