\name{RowOuterProduct}
\alias{RowOuterProduct}
\title{Outer product of matrix rows}
\usage{
  RowOuterProduct(m)
}
\arguments{
  \item{m}{(numeric matrix) Matrix whose rows we desire to
  expand.  Usually a spline matrix, where the columns
  correspond to the knots, and the rows correspond to
  evaluation points).}
}
\value{
  A 3-index numeric array, \code{value}, such that
  \code{value[, , i]} is the outer product of \code{m[i, ]}
  with itself.
}
\description{
  Computes the (3d) array consisting of the outer product
  of each row of m with itself.
}

