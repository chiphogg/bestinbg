
#' Eta: maximum straightness parameter
#'
#' This is basically a \dQuote{patch} to the background correction approach of
#' Fischer et al. (2000).  Their stated prior implicitly assigns infinite
#' probability to perfectly-straight backgrounds: a \dQuote{straight-line
#' trap}.  To fix this, we use a \dQuote{maximum straightness} parameter
#' \code{eta}.  It turns out to depend only on three parameters: the domain of
#' the background function, the distance of closest approach for spline knots,
#' and a vertical lengthscale.
#'
#' The vertical lengthscale \code{sigma} must be chosen manually by the user.
#' However, \dQuote{vertical lengthscales} are more intuitively accessible for
#' researchers than \dQuote{integrated second derivatives}.
#'
#' The other parameters, \code{dx.min} and \code{dx.total}, are already fixed
#' by the Fischer model.
#'
#' Simple dimensional analysis reveals that \code{eta} can be related to a
#' "universal" function, which depends only on the \dQuote{width} (basically,
#' the maximum number of spline knots which could be used).
#'
#' @rdname eta
#'
#' @param E (positive integer) The number of spline knots.
#' @param width (positive numeric) The width of the background function's
#'    domain, in units of the minimum knot spacing.  (e.g., if the background
#'    function goes from x=1 to x=9, and spline knots must be separated by at
#'    least 2, then \code{width} is 4.)
#'
#' @references  R. Fischer, K.M. Hanson, V. Dose, and W. von der Linden,
#' \dQuote{Background estimation in experimental spectra,} Physical Review E,  vol.
#' 61, Feb. 2000, pp. 1152-1160.
#'
#' @examples
#' \dontrun{
#'    # Compute the maximum straightness parameter for a background function
#'    # whose domain is 15 units wide,
#'    # whose spline knots can't be within 2 units of each other,
#'    # and whose straightness is relevant *only* to within *vertical*
#'    #   perturbabions of 0.1 units:
#'    eta <- Eta(sigma=0.1, dx.min=2, dx.total=15)
#' }
RandomSplinePositions <- function(E, width) {
  num.regions <- E - 1
  free.space <- width - num.regions
  xi <- 0:num.regions
  frac.free.space.added <- c(0, sort(runif(n=(E - 2), min=0, max=1)), 1)
  return (xi + frac.free.space.added * free.space)
}

#' @rdname eta
#' @param n (positive integer) The number of samples to use when estimating the
#'    universal eta function.
EtaUniversal <- function(width, n=2000) {
  if (width < 2) {
    return (0)
  }
  E.vals <- floor(runif(n=n, min=3, max=(width + 1)))
  traces <- 0 * (1:n)
  for (i in 1:n) {
    xi <- RandomSplinePositions(E=E.vals[i], width=width)
    traces[i] <- DMatrix(x=xi, only.trace=TRUE)
  }
  return (mean(traces))
}

