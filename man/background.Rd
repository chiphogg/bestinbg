\name{new_background}
\alias{background}
\alias{new_background}
\title{Background estimate}
\usage{
  new_background(delta.x, lambda, beta, epsilon, E, data)
}
\arguments{
  \item{delta.x}{(numeric) The minimum allowed separation
  distance between two spline knots (larger \code{delta.x}
  yields smoother backgrounds)}

  \item{lambda}{(numeric) Estimated average signal value
  (above the background)}

  \item{beta}{(numeric) Probability for a point to contain
  only background}

  \item{epsilon}{(numeric) \dQuote{Vertical} precision: we
  don't worry if spline knots fluctuate by this amount or
  less.}

  \item{E}{(numeric) Number of spline knots to use}

  \item{data}{(\code{\link{dataset}}) The dataset to
  analyze}
}
\description{
  A set of spline knots defining a background estimate,
  together with hyperparameters encoding underlying
  assumptions.
}
\details{
  If \code{epsilon} is zero, there is a \dQuote{straight
  line trap}: the model assigns infinite probability to
  perfectly straight lines.  A small but nonzero
  \code{epsilon} circumvents this: the model won't
  discriminate among straight lines which vary (vertically)
  by this amount or less.
}

