\name{NoisyTestFunction}
\alias{NoisyTestFunction}
\title{Random function with smooth background}
\usage{
  NoisyTestFunction(x, lambda, sigma, delta.x, E, n.peaks)
}
\arguments{
  \item{x}{(numeric vector) The x-points where data should
  be generated}

  \item{lambda}{(numeric; positive) The mean signal
  magnitude (approximate!)}

  \item{sigma}{(numeric; positive) The noise level to
  generate}

  \item{delta.x}{(numeric; positive) The minimum spacing
  allowed between spline knots}

  \item{E}{(numeric; positive integer) The number of spline
  knots}

  \item{n.peaks}{(numeric; positive integer) The number of
  peaks to generate

  There are some obvious constraints on these parameters.
  For instance, \code{(E - 1) * delta.x < diff(range(x))}
  should always be satisfied.}
}
\description{
  Generates a random function with a spline background,
  random peaks, and Gaussian noise.
}

