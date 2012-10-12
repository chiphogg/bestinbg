\name{DetD}
\alias{DetD}
\title{Computes the *modified* determinant (i.e., product of highest (N-2)
eigenvalues) for a D-matrix (see Fischer et al. for an explanation of why we
would want this).}
\usage{
  DetD(m)
}
\arguments{
  \item{m}{(numeric matrix) Overlap integrals for second
  derivatives of spline basis functions}
}
\value{
  The product of the largest all-but-two eigenvalues of m
  (the smallest two should be identically zero, although
  they won't due to roundoff errors).
}
\description{
  Computes the *modified* determinant (i.e., product of
  highest (N-2) eigenvalues) for a D-matrix (see Fischer et
  al. for an explanation of why we would want this).
}

