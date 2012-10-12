\name{OptimizeKnotValues}
\alias{OptimizeKnotValues}
\title{Compute background values at spline knots}
\usage{
  OptimizeKnotValues(bgr, data)
}
\arguments{
  \item{bgr}{(background object) The current background
  estimate (includes all relevant hyperparameters).}

  \item{data}{(dataset object) The data whose background to
  estimate.}
}
\value{
  (background object) The most probable spline knot values
}
\description{
  Compute background values at spline knots
}

