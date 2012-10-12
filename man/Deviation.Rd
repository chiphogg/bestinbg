\name{Deviation}
\alias{Deviation}
\title{Datapoint deviation from background}
\usage{
  Deviation(data, bgr)
}
\arguments{
  \item{data}{(\code{\link{dataset}} object) The data whose
  background to estimate}

  \item{bgr}{(\code{\link{background}} object) The
  estimated background (and associated hyperparameters)}
}
\value{
  (numeric vector) The difference between the measured data
  and the estimated background
}
\description{
  Datapoint deviation from background
}

