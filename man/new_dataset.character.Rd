\name{new_dataset.character}
\alias{new_dataset.character}
\title{New dataset from datafile}
\usage{
  new_dataset.character(x, sigma = NA, columns = 1:3, ...)
}
\arguments{
  \item{x}{(character) The name of a datafile which has 2
  columns (3 if sigma is included)}

  \item{sigma}{(numeric) Gaussian uncertainty (optionally
  different for each datapoint)}

  \item{columns}{(numeric vector) The column numbers
  corresponding to x, y, and sigma, respectively (the
  column for sigma may be omitted if \code{sigma} is
  supplied as an argument).}

  \item{...}{Other parameters for \code{\link{read.table}}}
}
\description{
  Create a new \code{dataset} object by reading from a
  datafile.
}
\details{
  This function expects at least two columns (three if
  \code{sigma} is not passed).  \code{columns} is a numeric
  vector of length 2 (or 3) indicating which datafile
  columns contain the data. \describe{
  \item{\code{columns[1]}}{The x-values}
  \item{\code{columns[2]}}{The y-values}
  \item{\code{columns[3]}}{The sigma-values (unless
  \code{sigma} is supplied)} }
}

