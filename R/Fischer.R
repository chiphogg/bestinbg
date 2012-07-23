#!/usr/bin/R

library("debug")
library("ggplot2")
library("Cairo")
library("trust")
library("colorspace")

##############################################################################
# UTILITY FUNCTIONS                                                          #
##############################################################################

invert.order <- function(i) {
  # Inverts the result of the 'order' function.
  #
  # Args:
  #   i:  Numeric vector with a permutation of the integers from 1:length(i).
  #
  # Returns:
  #   Numeric vector r such that r[i] == 1:length(i).
  r <- i + NA
  for (j in 1:length(i)) {
    r[i[j]] <- j
  }
  return (r)
}

Dx <- function(x) {
  # Compute the width for each x (specifically, its Voronoi cell size) to aid
  # in numerical integration.
  #
  # Args:
  #   x:  A sorted numeric vector of x-values
  n <- length(x)
  i <- order(x)
  x.sort <- x[i]
  dx <- diff(c(x.sort[1], 0.5 * (x.sort[-1] + x.sort[-n]), x.sort[n]))
  return (dx[invert.order(i)])
}

LayoutNewGridPage <- function(Layout, ...) {
  # Setup a new grid page with the specified Layout.  
  # (Adapted from http://bit.ly/9m4zyD)
  #
  # Args:
  #   Layout:  Result of calling grid.layout function with the desired
  #      arrangement.
  #
  # Returns:
  #   Used for its side-effect.
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}

Subplot <- function(x, y) {
  # Terse wrapper to return a viewport into the specified row (x) and 
  # column (y) of a viewport.
  # (Adapted from http://bit.ly/9m4zyD)
  #
  # Args:
  #   x:  Numeric; the row(s) of the layout.
  #   y:  Numeric; the column(s) of the layout.
  #
  # Returns:
  #   Used for its side-effect.
  viewport(layout.pos.row=x, layout.pos.col=y)
}

PlotBackgroundGuess <- function(curves, xi, cc, label="Background guess") {
  # Plots noisy data, the true background curve, and an estimate for the
  # background curve.
  #
  # Args:
  #   curves:  Something like the output of NoisyTestFunction()$curves: a
  #      data.frame having columns $x, $total, and $bgr (as well as $peaks and
  #      $noise).
  #   xi:  The spline knot x-values.
  #   cc:  The spline knot y-values ('c' in Fischer et al.).
  #   label:  The title for the plot.
  #
  # Returns:
  #   a ggplot() object with the desired plot.
  curve.names <- c("total", "bgr")
  knots.df <- data.frame(x=xi, value=cc, variable="bgr.guess")
  melted.curves <- rbind(
    melt(curves[, c("x", curve.names)], id.vars="x"),
    data.frame(x=curves$x, value=spline(x=xi, y=cc, xout=curves$x)$y,
      variable="bgr.guess"))
  p <- (ggplot(data=melted.curves, aes(x=x, y=value, colour=variable))
    + geom_line()
    + geom_point(data=knots.df)
    + scale_colour_manual(values=c(
        total="black", bgr="red", bgr.guess="blue"))
    + opts(title=label)
    )
  return (p)
}

##############################################################################
# FUNCTIONS RELATING TO FISCHER ET AL.                                       #
##############################################################################

WiggleRoom <- function(E, min.dx, x.range) {
  x.range - (E - 1) * min.dx
}

WiggleRoomFractions <- function(xi, min.dx) {
  shares <- diff(xi) - min.dx
  fractions <- shares / sum(shares)
  return (fractions)
}

Xi <- function(fracts, min.dx, xi.min, xi.max) {
  x.range <- xi.max - xi.min
  room <- WiggleRoom(E=(length(fracts) + 1), min.dx=min.dx, x.range=x.range)
  spaces <- fracts * room
  return (xi.min + cumsum(c(0, (spaces + min.dx))))
}

spline.basis.function <- function(x, spline.x, i, deriv=0) {
  # Computes the i'th spline basis function for taking spline.x into x.
  #
  # Args:
  #   x:  Numeric vector giving x-values where we evaluate the basis functions.
  #   spline.x:  Numeric vector of spline knot positions
  #   i:  Integer from 1 to length(spline.x): which basis function to return?
  #   deriv:  Which derivative of the spline function? (0 to 3; 0 is function
  #     itself)
  #
  # Returns:
  #   A numeric vector of length length(x): the requested derivative of the
  #   requested basis function, evaluated at the points in x.
  y <- 0 * spline.x
  y[i] <- 1
  spline.fun <- splinefun(x=spline.x, y=y, method="natural")
  return (spline.fun(x=x, deriv=deriv))
}

spline.matrix <- function(x, spline.x, deriv=0) {
  # Computes the matrix whose columns are spline basis functions taking
  # spline.x into x.
  #
  # Args:
  #   x:  Numeric vector giving x-values where we evaluate the basis functions.
  #   spline.x:  Numeric vector of spline knot positions
  #   deriv:  Which derivative of the spline function? (0 to 3; 0 is function
  #     itself)
  #
  # Returns:
  #   The matrix whose columns are spline basis functions taking spline.x into
  #   x.
  N <- length(spline.x)
  N.out <- length(x)
  M <- matrix(NA, nrow=N.out, ncol=N)
  for (i in 1:N) {
    M[, i] <- spline.basis.function(x=x, spline.x=spline.x, i=i, deriv=deriv)
  }
  return (M)
}

D.matrix <- function(x, spline.x, robust.factor=1e-15) {
  # Computes the second derivative overlap matrix ('D' in Fischer et al.
  # (1999)) of the spline basis functions taking spline.x into x.
  #
  # Args:
  #   x:  Numeric vector giving x-values where we evaluate the basis functions.
  #   spline.x:  Numeric vector of spline knot positions
  #   robust.factor:  We add a constant to all eigenvalues to promote
  #      stability; this is the ratio of that constant to the smallest
  #      eigenvalue which is SUPPOSED to be nonzero.
  #
  # Returns:
  #   A matrix DD where DD(i, j) gives the integral of the product of the
  #   second derivatives of spline basis functions i and j.
  N <- length(spline.x)
  N.out <- length(x)
  ddM <- spline.matrix(x=x, spline.x=spline.x, deriv=2)
  d.x <- Dx(x)
  DD <- matrix(NA, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in 1:N) {
      DD[i, j] <- sum(ddM[, i] * ddM[, j] * d.x)
    }
  }
  eigenvals <- eigen(x=DD, symmetric=TRUE, only.values=TRUE)$values
  lowest.true.eigenval <- min(sort(eigenvals)[-(1:2)])
  return (DD + lowest.true.eigenval * robust.factor * diag(nrow(DD)))
}

DMatrix <- function(spline.x, only.trace=FALSE, robust.factor=1e-12) {
  # Computes the second derivative overlap matrix ('D' in Fischer et al.
  # (1999)) of the spline basis functions with knots at 'spline.x'.
  #
  # Args:
  #   spline.x:  Numeric vector of spline knot positions
  #   only.trace:  If true, return trace(D) (a single number) instead of D.
  #   robust.factor:  We add a constant to all eigenvalues to promote
  #      stability; this is the ratio of that constant to the smallest
  #      eigenvalue which is SUPPOSED to be nonzero.
  #
  # Returns:
  #   A matrix DD where DD(i, j) gives the integral of the product of the
  #   second derivatives of spline basis functions i and j.
  N <- length(spline.x)
  ddM <- spline.matrix(x=spline.x, spline.x=spline.x, deriv=2)
  x.mat <- matrix(rep(spline.x, N), nrow=N)
  dy <- apply(ddM, 2, diff)
  dx <- apply(x.mat, 2, diff)
  # slope[i, j] is the slope of the i'th segment of the j'th basis function;
  # intercept[i, j] is its y-intercept
  slope <- dy / dx
  intercept <- ddM[-1, ] - slope * x.mat[-1, ]
  d.x3 <- diff(spline.x ^ 3)
  d.x2 <- diff(spline.x ^ 2)
  # The for-loop code below reduces to this if we only care about diagonal
  # elements (which is the case for computing the trace).
  if (only.trace) {
    return (sum((slope ^ 2) * d.x3 / 3.0 + (slope * intercept) * d.x2
        + intercept ^ 2 * diff(spline.x)))
  }
  # If we're this far, we need to compute the whole matrix.
  DD <- matrix(NA, nrow=N, ncol=N)
  for (i in 1:N) {
    for (j in i:N) {
      prod.ss <- slope[, i] * slope[, j] / 3.0
      prod.si <- 0.5 * (
        slope[, i] * intercept[, j] + slope[, j] * intercept[, i])
      prod.ii <- intercept[, i] * intercept[, j]
      integrals <- prod.ss * d.x3 + prod.si * d.x2 + prod.ii * diff(spline.x)
      DD[i, j] <- DD[j, i] <- sum(integrals)
    }
  }
  # Make sure we avoid negative eigenvalues!
  eigenvals <- eigen(x=DD, symmetric=TRUE, only.values=TRUE)$values
  lowest.true.eigenval <- sort(eigenvals)[3]
  return (DD + lowest.true.eigenval * robust.factor * diag(nrow(DD)))
}

RandomSplinePositions <- function(E, width) {
  num.regions <- E - 1
  free.space <- width - num.regions
  xi <- 0:num.regions
  frac.free.space.added <- c(0, sort(runif(n=(E - 2), min=0, max=1)), 1)
  return (xi + frac.free.space.added * free.space)
}

EtaUniversal <- function(width, n=2000) {
  if (width < 2) {
    return (0)
  }
  E.vals <- floor(runif(n=n, min=3, max=(width + 1)))
  traces <- 0 * (1:n)
  for (i in 1:n) {
    xi <- RandomSplinePositions(E=E.vals[i], width=width)
    traces[i] <- DMatrix(spline.x=xi, only.trace=TRUE)
  }
  return (mean(traces))
}

Eta <- function(sigma, dx.min, dx.total) {
  return (sigma ^ 2 * EtaUniversal(width=(dx.total / dx.min)) / (dx.min ^ 3))
}

RowOuterProduct <- function(Phi) {
  # Computes the (3d) array consisting of the outer product of each row of Phi
  # with itself.
  #
  # Args:
  #   Phi:  A (R x C) numeric matrix (usually a spline matrix, where the
  #      columns correspond to the knots, and the rows correspond to evaluation
  #      points).
  #
  # Returns:
  #   A 3-index numeric array, whose dimensions are c(C, C, R), such that the
  #   matrix at [*, *, i] is the outer product of Phi[i, ] with itself.
  rows <- nrow(Phi)
  cols <- ncol(Phi)
  ppt <- array(apply(Phi, 1, function(x) x %*% t(x)), dim=c(cols, cols, rows))
  return (ppt)
}

ReduceFrom3D <- function(H, v, premultiply=FALSE) {
  # Reduces 3D array H to a 2D matrix by multiplying with vector v.
  #
  # Args:
  #   H:  a 3D numeric array; the 3rd index indicates which 2D matrix we're
  #      dealing with.
  #   v:  A vector for reducing H to a matrix.
  #
  # Returns:
  #   If v is a vector, returns the matrix M such that:
  #      M[k, i] = sum_j H[i, j, k] v[j]
  if (premultiply) {
    result <- apply(H, 3, function(x) t(v) %*% x)
  }
  else {
    result <- apply(H, 3, function(x) x %*% v)
  }
  return (t(result))
}

det.D <- function(DD) {
  # Computes the *modified* determinant (i.e., product of highest (N-2)
  # eigenvalues) for a D-matrix (see Fischer et al. for an explanation of why
  # we would want this).
  #
  # Args:
  #   DD:  A D-matrix (giving overlap integrals for second derivatives of
  #      spline basis functions)
  #
  # Returns:
  #   The product of the largest all-but-two eigenvalues of DD (the smallest
  #   two should be identically zero, although they won't due to roundoff
  #   errors).
  eigenvals <- eigen(x=DD, symmetric=TRUE, only.values=TRUE)$values
  return (prod(rev(eigenvals)[-(1:2)]))
}

SafetyNetEpsilon <- function(safety.net) {
  # Computes the fraction of the noise mixture model which is exponential
  # rather than Gaussian.
  #
  # Explanation:
  #   A purely-Gaussian noise model has trouble with outliers.  In log-space,
  #   it drops off quadratically.  Instead, define a mixture model which is
  #   mostly Gaussian but has exponential tails.  This makes it linear in
  #   log-space, which is much gentler (while still encouraging progress
  #   towards the peak).
  #
  # Args:
  #   safety.net:  (numeric) After how many sigma's should the exponential
  #      tails "take over"?
  #
  # Returns:
  #   Epsilon.
  return (1 / (1 + exp(0.5 * (safety.net - 1) ^ 2) * sqrt(0.5 * pi / exp(1))))
}

PsiFBgr <- function(cc, d, Phi, p.bgr, sigma, calc.grad, calc.hess, ...) {
  # The 'f'-function: i.e., the background-only contribution to the likelihood.
  #
  # Args:
  #   cc:  'c' in Fischer et al.: vector of spline knot y-values
  #   d:  vector of datapoints
  #   Phi:  spline matrix taking 'c' into 'd'
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  vector of experimental uncertainties
  #   calc.grad:  If TRUE, also calculate the gradient
  #   calc.hess:  If TRUE, also calculate the Hessian
  #   safety.net:  After this many sigmas, an exponential distribution takes
  #      over from the Gaussian distribution, to help the numerical optimizer.
  #
  # Returns:
  #   A list() with three elements:
  #     $funct:  vector of values for the (background-only) likelihood of each
  #        datapoint.
  #     $grad:  2-index array, where $grad[i, j] is the jth component of the
  #        gradient (w.r.t. cc) of $funct[i]
  #     $hess:  3-index array, where $hess[j, k, i] gives the (j, k)th
  #        component of the Hessian (w.r.t. cc) of $funct[i] (NOTE different
  #        ordering of index i; it's unfortunate, but this is a more "natural"
  #        ordering in R for multi-D arrays!)
  deviation <- (d - Phi %*% cc)
  norm.dev <- deviation / sigma
  funct <- (log(p.bgr) - 0.5 * log(2 * pi) - log(sigma) - 0.5 * norm.dev ^ 2)
  grad <- hess <- NA
  if (calc.grad) {  # Calculate the gradient
    grad <- as.vector(deviation / (sigma ^ 2)) * Phi
  }
  if (calc.hess) {  # Calculate the Hessian
    PPT <- RowOuterProduct(Phi)
    sigma.rep <- array(rep(sigma, each=length(cc) ^ 2), dim=dim(PPT))
    hess <- -PPT / (sigma.rep ^ 2)
  }
  return (list(funct=funct, grad=grad, hess=hess))
}

PsiHSig <- function(cc, d, Phi, p.bgr, sigma, lambda, calc.grad, calc.hess, ...) {
  # The 'h'-function: i.e., the signal-containing contribution to the
  # likelihood.
  #
  # Args:
  #   cc:  'c' in Fischer et al.: vector of spline knot y-values
  #   d:  vector of datapoints
  #   Phi:  spline matrix taking 'c' into 'd'
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  vector of experimental uncertainties
  #   lambda:  mean signal magnitude for signal-containing pixels
  #   calc.grad:  If TRUE, also calculate the gradient
  #   calc.hess:  If TRUE, also calculate the Hessian
  #
  # Returns:
  #   A vector of values for the (signal-containing) likelihood of each
  #   datapoint.
  rho <- sigma / lambda
  z <- (d - (Phi %*% cc)) / lambda
  qq <- z / rho - rho
  funct <- (log(1 - p.bgr) - log(lambda) + pnorm(log.p=TRUE, q=qq) - z
    + 0.5 * (rho ^ 2))
  grad <- hess <- gamma.q <- NA
  GammaQ <- function(qq, rho) {
    # This combination of variables appears frequently in calculations for the
    # gradient and Hessian.
    gq <- exp(-0.5 * qq ^ 2 - pnorm(q=qq, log.p=TRUE))
    return (gq / (rho * sqrt(2 * pi)))
  }
  if (calc.grad) {  # Calculate the gradient
    gamma.q <- GammaQ(qq=qq, rho=rho)
    grad <- as.vector((1 - gamma.q) / lambda) * Phi
  }
  if (calc.hess) {  # Calculate the Hessian
    prefactor <- as.vector(-gamma.q * (gamma.q + qq / rho) / (lambda ^ 2))
    PPT <- RowOuterProduct(Phi)
    prefactor.rep <- array(rep(prefactor, each=length(cc) ^ 2), dim=dim(PPT))
    hess <- prefactor.rep * PPT
  }
  return (list(funct=funct, grad=grad, hess=hess))
}

Psi <- function(x, y, cc, xi, p.bgr=0.5, sigma, lambda, eta,
  Phi=spline.matrix(x=x, spline.x=xi), calc.grad, calc.hess) {
  # Computes the optimization function from Fischer et al.: the negative log of
  # the posterior.
  #
  # Args:
  #   x:  The x-values of the measured datapoints.
  #   y:  The y-values of the measured datapoints.
  #   cc:  The spline knot y-values ('c' in Fischer et al.).
  #   xi:  The spline knot x-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve
  #   calc.grad:  If TRUE, also calculate the gradient
  #   calc.hess:  If TRUE, also calculate the Hessian
  #
  # Returns:
  #   A list() with three elements:
  #     $funct:  The negative log of the posterior.
  #     $grad:  Row vector (or NA); the gradient of funct with respect to cc
  #     $hess:  Matrix (or NA); the Hessian of funct with respect to cc
  E <- length(cc)
  if (length(xi) != E) {
    stop("Inconsistent lengths of parameters cc and xi!")
  }
  DD <- DMatrix(spline.x=xi)
  cDc <- as.vector(t(cc) %*% DD %*% cc) + eta
  # Prior: separate terms according to which variables contribute
  contr.E       <- (0.5 * E * log(pi) - log(gamma(0.5 * E)) 
    - sum(log(1:(E - 2))))  # "Depends" on xi, but constant as long as xi OK.
  contr.E.xi    <- (-0.5 * log(det.D(DD)))
  contr.E.xi.cc <- 0.5 * E * log(cDc)
  psi.prior <- contr.E + contr.E.xi + contr.E.xi.cc
  # Likelihood:
  calc.grad <- calc.grad || calc.hess  # Hessian requires gradient!
  f <- PsiFBgr(cc=cc, d=y, Phi=Phi, p.bgr=p.bgr, sigma=sigma,
    calc.grad=calc.grad, calc.hess=calc.hess)
  h <- PsiHSig(cc=cc, d=y, Phi=Phi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
    calc.grad=calc.grad, calc.hess=calc.hess)
  max.log <- pmax(f$funct, h$funct)
  f.fract <- as.vector(1 / (1 + exp(h$funct - f$funct)))
  g <- max.log + log(exp(f$funct - max.log) + exp(h$funct - max.log))
  psi.likelihood <- -sum(g)
  # Calculate higher derivatives, if requested
  grad <- hess <- NA
  if (calc.grad) {  # Calculate the gradient
    grad.g <- f.fract * f$grad + (1 - f.fract) * h$grad
    grad <- -colSums(grad.g) + (t(cc) %*% DD) * E / cDc
  }
  if (calc.hess) {  # Calculate the Hessian
    f.contr <- (    f.fract) * (f$hess + RowOuterProduct(Phi=f$grad))
    h.contr <- (1 - f.fract) * (h$hess + RowOuterProduct(Phi=h$grad))
    hess.g <- f.contr + h.contr - RowOuterProduct(Phi=grad.g)
    hess <- apply(hess.g, c(1, 2), sum) + E * (
      (DD / cDc) - (2 * DD %*% cc %*% t(cc) %*% DD) / (cDc ^ 2)
      )
  }
  return (list(funct=psi.prior + psi.likelihood, grad=grad, hess=hess))
}

PsiWrapped <- function(params, x, y, xi, p.bgr=0.5, sigma, lambda, eta,
  Phi=spline.matrix(x=x, spline.x=xi), ...) {
  # Wrap Psi(cc, xi) for use with trust so that it returns gradient and
  # hessian.
  #
  # Args:
  #   params:  The spline knot y-values ('c' in Fischer et al.).
  #   x:  The x-values of the measured datapoints.
  #   y:  The y-values of the measured datapoints.
  #   xi:  The spline knot x-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve
  #
  # Returns:
  #   A list() with elements:
  #     $value:  (numeric) The value of psi
  #     $gradient:  (numeric vector) The gradient of psi with respect to cc
  #     $hessian:  (numeric matrix) The Hessian of psi with respect to cc
  f <- Psi(x=x, y=y, cc=params, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
    eta=eta, Phi=Phi, calc.grad=TRUE, calc.hess=TRUE)
  return (list(value=f$funct, gradient=as.vector(f$grad), hessian=f$hess))
}

OptimizeCC <- function(x, y, cc=(0 * xi), xi, p.bgr, sigma, lambda, eta,
  Phi=spline.matrix(x=x, spline.x=xi)) {
  # Optimizes cc (i.e., the spline values; 'c' in Fischer et al.) with a given
  # level of experimental noise.
  #
  # Args:
  #   x:  The x-values of the measured datapoints.
  #   y:  The y-values of the measured datapoints.
  #   cc:  Initial guess for the spline knot y-values ('c' in Fischer et al.).
  #   xi:  The spline knot x-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve
  #
  # Returns:
  #   The optimized values of cc.
  opt <- trust(objfun=PsiWrapped, parinit=cc, rinit=1e-2, rmax=1e6,
    blather=FALSE,
    x=x, y=y, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, Phi=Phi)
  return (opt$argument)
}

OptimizeCCShrinkingSigma <- function(curves, cc=(0 * xi), xi, p.bgr, sigma,
  lambda, eta, Phi=spline.matrix(x=x, spline.x=xi)) {
  # Optimizes cc (i.e., the spline values; 'c' in Fischer et al.) by gradually
  # decreasing the experimental noise, from lambda down to its true value of
  # sigma.
  #
  # Args:
  #   curves:  Something like the output of NoisyTestFunction()$curves: a
  #      data.frame having columns $x, $total, and $bgr (as well as $peaks and
  #      $noise).
  #   cc:  Initial guess for the spline knot y-values ('c' in Fischer et al.).
  #   xi:  The spline knot x-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve
  #
  # Returns:
  #   A numeric vector containing optimal values for cc.
  sigma.current <- pmax(sigma, max(lambda, sigma))
  i <- 0
  repeat {
    cc <- OptimizeCC(x=curves$x, y=curves$total, cc=cc, xi=xi, p.bgr=p.bgr,
      sigma=sigma.current, lambda=lambda, eta=eta, Phi=Phi)
    if (length(which(sigma.current > sigma)) <= 0) {
      break
    }
    sigma.current <- pmax(sigma.current / 2, sigma)
    i <- i + 1
  }
  return (cc)
}

DirichletDraw <- function(alpha) {
  gamma.draws <- rgamma(n=length(alpha), shape=alpha)
  return (gamma.draws / sum(gamma.draws))
}

LogDDirichlet <- function(x, alpha) {
  # Computes the log of the probability density for x, given a Dirichlet
  # distribution parameterized by alpha.
  #
  # Args:
  #   x:  The values where we want to evaluate the pdf.
  #   alpha:  The parameters of the Dirichlet distribution.
  #
  # Returns:
  #   (numeric): log[p.gamma(x | alpha)]
  return (sum((alpha - 1) * log(x) - lgamma(alpha)) + lgamma(sum(alpha)))
}

DrawNewXi <- function(xi, min.dx, concentration) {
  # Draws new values for xi, given the old values of xi, the closest allowable
  # distance, and the 'concentration' (higher means smaller jumps).
  #
  # Returns:
  #   A list() of two elements:
  #     $xi:  (numeric) New xi values.
  #     $d.log.p:  (numeric) Log of proposal probability ratio,
  #       log[ p(new->old) / p(old->new) ]
  fracts <- WiggleRoomFractions(xi=xi, min.dx=min.dx)
  alpha <- 1 + concentration * fracts
  new.fracts <- DirichletDraw(alpha)
  new.alpha <- 1 + concentration * new.fracts
  xi.new <- Xi(fracts=new.fracts, min.dx=min.dx,
    xi.min=min(xi), xi.max=max(xi))
  log.p.goto.old <- LogDDirichlet(x=fracts, alpha=new.alpha)
  log.p.goto.new <- LogDDirichlet(x=new.fracts, alpha=alpha)
  return (list(xi=xi.new, d.log.p=(log.p.goto.old - log.p.goto.new)))
}

OptimizeXi <- function(curves, xi, min.dx, p.bgr, sigma, lambda, eta, 
  temp.init=500, reduce.new.best=0.95, reduce=0.995, c.at.1=1000) {
  # Optimizes xi (i.e., the spline locations) by simulated annealing Markov
  # Chain Monte Carlo (MCMC).
  #
  # Args:
  #   curves:  Something like the output of NoisyTestFunction()$curves: a
  #      data.frame having columns $x, $total, and $bgr (as well as $peaks and
  #      $noise).
  #   xi:  The starting spline knot x-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   temp.init:  Initial "temperature" for simulated annealing.
  #   reduce.new.best:  The factor to reduce the temperature when we find a new
  #      best value of psi
  #   reduce:  The factor to reduce the temperature when we accept a move which
  #      is *not* the best-so-far.
  #   c.at.1:  The (c)oncentration at temperature=1, i.e., a number which
  #      describes how closely the new values should stick to the present ones
  #      (higher means tighter; any positive number will do).
  #
  # Returns:
  #   The optimal xi values.
  temp <- temp.init
  cc <- OptimizeCCShrinkingSigma(curves=curves, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta)
  psi.best <- Psi(x=curves$x, y=curves$total, cc=cc, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, calc.grad=FALSE, calc.hess=FALSE)
  xi.traces <- data.frame(xi=matrix(xi, nrow=1), cc=matrix(cc, nrow=1),
    temp=temp, p=1, psi=psi.best$funct)
  psi.new <- psi.best
  while (temp > 1) {
    # Draw new xi-values (and calculate proposal probability ratio!)
    xi.new <- DrawNewXi(xi=xi, min.dx=min.dx, concentration=(c.at.1 / temp))
    # Optimize c-values for the new xi
    cc <- OptimizeCCShrinkingSigma(curves=curves, cc=cc, xi=xi.new$xi,
      p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)
    psi.prev <- psi.new
    psi.new <- Psi(x=curves$x, y=curves$total, cc=cc, xi=xi.new$xi,
      p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, calc.grad=FALSE,
      calc.hess=FALSE)
    # Calculate acceptance probability
    prob <- exp(((psi.new$funct - psi.prev$funct) / temp) + xi.new$d.log.p)
    accepted <- (prob > runif(1))
    # Set new value of xi and adjust the temperature
    temp.prev <- temp
    if (accepted) {
      xi <- xi.new$xi
      if (psi.new$funct < psi.best$funct) {
        psi.best <- psi.new
        temp <- temp * reduce.new.best
      } else {
        temp <- temp * reduce
      }
    }
    xi.traces <- rbind(xi.traces, data.frame(xi=matrix(xi, nrow=1),
        cc=matrix(cc, nrow=1), temp=temp.prev, p=prob, psi=psi.new$funct))
    save(file="xi_traces.ROBJ", xi.traces, curves)
  }
  return (xi.traces)
}

QuenchedXiOptimization <- function(curves, xi, min.dx, p.bgr, sigma, lambda, eta, 
  temp.init=500, reduce.new.best=0.95, reduce=0.995, c.at.1=1000) {
  xi.traces <- OptimizeXi(curves=curves, xi=xi, min.dx=min.dx, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, temp.init=temp.init,
    reduce.new.best=reduce.new.best, reduce=reduce, c.at.1=c.at.1)
  iter <- 1
  all.xi.traces <- data.frame()
  repeat {
    xi.traces$iter <- iter
    all.xi.traces <- rbind(all.xi.traces, xi.traces)
    save(file="all_xi_traces.ROBJ", all.xi.traces, curves)
    i.best <- which(all.xi.traces$psi == min(all.xi.traces$psi))
    temp.quench <- all.xi.traces$temp[i.best] * 0.5
    if (temp.quench < 1) {
      break;
    }
    iter <- iter + 1
    i.xi <- grep("^xi\\.\\d+$", names(all.xi.traces))
    xi <- unname(unlist(all.xi.traces[i.best, i.xi]))
    n.knots <- length(xi)
    xi.traces <- OptimizeXi(curves=curves, xi=xi, min.dx=min.dx, p.bgr=p.bgr,
      sigma=sigma, lambda=lambda, eta=eta, temp.init=temp.quench,
      reduce.new.best=reduce.new.best, reduce=reduce, c.at.1=c.at.1)
  }
  return (all.xi.traces)
}

ExtractKnotPaths <- function(xi.traces) {
  # Creates a new data.frame with the timetrace for each knot, using the output
  # of OptimizeXi().
  #
  # Args:
  #   xi.traces:  (data.frame) The output of OptimizeXi().
  #
  # Returns:
  #   data.frame() with columns temp (temperature), knot (which knot), p
  #   (acceptance probability), xi (knot x-values), cc (knot y-values).
  n.knots <- length(grep("^xi\\.\\d+$", names(xi.traces)))
  knot.paths <- data.frame(temp=c(), knot=c(), p=c(), xi=c(), cc=c())
  for (k in 1:n.knots) {
    point.names <- paste(sep="", c("xi", "cc"), ".", k)
    this.knot.names <- c("temp", "log.p", "psi", point.names)
    this.knot <- xi.traces[, this.knot.names]
    this.knot$knot <- k
    this.knot <- this.knot[, c("temp", "knot", "log.p", "psi", point.names)]
    colnames(this.knot)[5:6] <- c("xi", "cc")
    knot.paths <- rbind(knot.paths, this.knot)
  }
  return (knot.paths)
}

VisualizeKnotPaths <- function(xi.traces, hilite.within=2) {
  ekp <- ExtractKnotPaths(xi.traces)
  ekp$size <- -pmin(hilite.within, ekp$psi - min(ekp$psi))
  i.best <- which(ekp$size >= 0)
  hi.breaks <- seq(from=0, to=-hilite.within, length.out=3)
  hi.labels <- hi.breaks
  hi.labels[1] <- "0 (best)"
  hi.labels[3] <- paste(hi.labels[3], "\n(or worse)")
  p <- (ggplot(data=ekp, aes(x=xi, y=log(temp), group=knot,
        colour=factor(knot), size=size))
    + geom_line(inherit.aes=FALSE, data=ekp[i.best, ], size=1,
      aes(x=xi, y=log(temp), group=temp))
    + geom_point()
    + scale_colour_brewer("Knot")
    + scale_size_continuous("Comparison\nto best", to=c(0.8, 8),
      labels=hi.labels, breaks=hi.breaks))
  return (p)
}

VisualizeAllKnotPaths <- function(all.xi.traces, label, hilite.within=2) {
  require("Cairo")
  y.range <- log(range(all.xi.traces$temp))
  for (iter in unique(all.xi.traces$iter)) {
    i.iter <- which(all.xi.traces$iter == iter)
    p <- VisualizeKnotPaths(xi.traces=all.xi.traces[i.iter, ],
      hilite.within=hilite.within)
    CairoPNG(filename=sprintf("quenched_%s_%02d.png", label, iter),
      width=800, height=800)
    print(p + ylim(y.range))
    dev.off()
  }
  p <- VisualizeKnotPaths(xi.traces=all.xi.traces,
    hilite.within=hilite.within)
  CairoPNG(filename=sprintf("quenched_%s_all.png", label),
    width=800, height=800)
  print(p + ylim(y.range))
  dev.off()
}

PlotAllXiGuesses <- function(xi.traces, curves, label, max.steps.shown=25) {
  knot.paths <- ExtractKnotPaths(xi.traces)
  p <- (ggplot(data=knot.paths, aes(x=xi, y=cc, group=as.factor(knot),
        colour=log(temp)))
    + geom_path()
    + geom_line(data=curves, aes(x=x, y=total), inherit.aes=FALSE)
    + scale_colour_gradientn(colour=heat_hcl(10)))
  # Now, plot the evolution of the background curve
  n.steps <- min(nrow(xi.traces), max.steps.shown)
  i.row <- round(seq(from=1, to=nrow(xi.traces), length.out=n.steps))
  i.xi <- grep("^xi\\.\\d+$", names(xi.traces))
  i.cc <- grep("^cc\\.\\d+$", names(xi.traces))
  y.range <- range(c(xi.traces[, i.cc], curves$total))
  log.temp.range <- log(range(xi.traces$temp))
  log.t.min <- min(log.temp.range)
  log.t.max <- max(log.temp.range)
  for (i in i.row) {
    knot.data <- data.frame(x=unlist(matrix(xi.traces[i, i.xi])),
      y=unlist(matrix(xi.traces[i, i.cc])))
    curve.data <- with(knot.data, data.frame(spline(x=x, y=y, xout=curves$x)))
    knot.data$temp <- curve.data$temp <- xi.traces$temp[i]
    CairoPNG(filename=sprintf("opt-xi_%s_%05d.png", label, i),
      width=800, height=800)
    p.1curve <- (ggplot(data=knot.data, aes(x=x, y=y, colour=log(temp)))
      + geom_point()
      + geom_line(data=curve.data)
      + geom_line(data=curves, aes(x=x, y=total), inherit.aes=FALSE)
      #+ geom_line(data=curves, aes(x=x, y=bgr), inherit.aes=FALSE,
      #  colour='blue')
      + ylim(y.range)
      + scale_colour_gradientn(colours=heat_hcl(10)[1:8], rescale=FALSE,
        limits=log.temp.range,
        values=seq(from=log.t.min, to=log.t.max, length.out=10))
      )
    print(p.1curve)
    dev.off()
  }
  return (p)
}

##############################################################################
# PDF-RELATED FUNCTIONS                                                      #
##############################################################################

NoiseCovarianceGFromSSingleArg <- function(r1, r2, Q, sigma) {
  # Computes the covariance between two points in G(r) due to noise in S(Q)
  # (the latter is assumed i.i.d. Gaussian).
  #
  # Args:
  #   r1:  (numeric) One r-value to consider
  #   r2:  (numeric) The other r-value to consider
  #   Q:  (numeric vector) The Q-values where we have data.
  #   sigma:  (numeric) The standard deviation of the noise in S.
  #
  # Returns:
  #   (numeric) The covariance between G(r1) and G(r2) due to noise in S(Q).
  delta <- diff(Q)
  delta <- c(delta[1], delta)
  f.sum <- sum((2 * sigma * Q * delta / pi) ^ 2 * sin(Q * r1) * sin(Q * r2))
  return (f.sum)
}
NoiseCovarianceGFromS <- Vectorize(NoiseCovarianceGFromSSingleArg,
  vectorize.args=c("r1", "r2"))

NoiseCovarianceMatrixGFromS <- function(r, Q, sigma) {
  # Computes the covariance matrix in G(r) due to noise in S(Q) (the latter is
  # assumed i.i.d. Gaussian).
  #
  # Args:
  #   r:  (numeric vector) The r-values where we evaluate this covariance.
  #   Q:  (numeric vector) The Q-values where we have data.
  #   sigma:  (numeric) The standard deviation of the noise in S.
  #
  # Returns:
  #   An (N x N) matrix M, where N is the length of r, such that M[i, j] gives
  #   the covariance between r[i] and r[j].
  N <- length(r)
  M <- matrix(nrow=N,
    NoiseCovarianceGFromS(r1=rep(r, each=N), r2=rep(r, N), Q=Q, sigma=sigma))
  return (M)
}

FourierSineTransformMatrixQR <- function(Q, r) {
  # Computes the matrix which converts a function at the given Q-points into
  # r-space using a Fourier sine transform.
  #
  # Args:
  #   Q:  (numeric vector) The Q-values where the function is evaluated.
  #   r:  (numeric vector) The r-values where the function is evaluated.
  #
  # Returns:
  #   A (N.r x N.Q) numeric matrix which effects the Q-to-r sine Fourier
  #   transform.
  dQ <- Dx(Q)
  N.Q <- length(Q)
  return (sapply(X=1:N.Q, Q=Q, dQ=dQ, r=r, FUN=function(i, Q, dQ, r) {
        sin(Q[i] * r) * dQ[i]
      }))
}

FourierSineTransformQR <- function(f.Q, Q, r) {
  # Computes the matrix which converts a function at the given Q-points into
  # r-space using a Fourier sine transform.
  #
  # Args:
  #   Q:  (numeric vector) The Q-values where the function is evaluated.
  #   r:  (numeric vector) The r-values where the function is evaluated.
  #
  # Returns:
  #   A (N.r x N.Q) numeric matrix which effects the Q-to-r sine Fourier
  #   transform.
  N.Q <- length(Q)
  N.r <- length(r)
  dQ <- Dx(Q)
  f.r <- 0 * r
  for (i in 1:N.Q) {
    for (j in 1:N.r) {
      f.r[j] <- f.r[j] + sin(Q[i] * r[j]) * f.Q[i] * dQ[i]
    }
  }
  return (f.r)
}

SofQ.fitness <- function(params, Q, F.exp, xi, p.bgr=0.5, sigma, lambda, eta,
  Phi, J, mu, ...) {
  # Computes the total fitness function for the spline values: including
  # contributions from the Q(S(Q) - 1) fit (i.e., Fischer et al.), as well as
  # the low-r G(r) constraints.
  #
  # Args:
  #   params:  The spline knot y-values ('c' in Fischer et al.).
  #   Q:  The Q-values of the measured datapoints.
  #   F.exp:  The y-values of the measured datapoints.
  #   xi:  The spline knot Q-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve.
  #   J:  (numeric matrix) Matrix for the quadratic contribution of cc to the
  #      fitness function.
  #   mu:  (numeric vector) Vector for the linear contribution of cc to the
  #      fitness function.
  #
  # Returns:
  #   A list() with elements:
  #     $value:  (numeric) The value of psi
  #     $gradient:  (numeric vector) The gradient of psi with respect to cc
  #     $hessian:  (numeric matrix) The Hessian of psi with respect to cc
  fit <- PsiWrapped(params=params, x=Q, y=F.exp, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, ...)
  cT.J <- t(params) %*% J
  fit$value <- fit$value + (0.5 * cT.J %*% params) - (t(mu) %*% params)
  fit$gradient <- fit$gradient + t(cT.J) - mu
  fit$hessian <- fit$hessian + J
  if (any(is.nan(fit$hessian))) {
    fit$hessian <- diag(length(xi))
  }
  return (fit)
}

SofQ.optim.fn <- function(params, Q, F.exp, xi, p.bgr=0.5, sigma, lambda,
  eta, Phi, J, mu, ...) {
  fit <- SofQ.fitness(params=params, Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu, ...)
  return (fit$value)
}

SofQ.optim.gr <- function(params, Q, F.exp, xi, p.bgr=0.5, sigma, lambda,
  eta, Phi, J, mu, ...) {
  fit <- SofQ.fitness(params=params, Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu, ...)
  return (fit$gradient)
}

SofQ.OptimizeCC <- function(Q, F.exp, xi, cc=(0 * xi), p.bgr, sigma, lambda, eta,
  Phi=spline.matrix(x=x, spline.x=xi), J, mu, ...) {
  # Optimizes cc (i.e., the spline values; 'c' in Fischer et al.) with a given
  # level of experimental noise.
  #
  # Args:
  #   Q:  The Q-values of the measured datapoints.
  #   F.exp:  The y-values of the measured datapoints.
  #   cc:  Initial guess for the spline knot y-values ('c' in Fischer et al.).
  #   xi:  The spline knot Q-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve
  #   J:  (numeric matrix) Matrix for the quadratic contribution of cc to the
  #      fitness function.
  #   mu:  (numeric vector) Vector for the linear contribution of cc to the
  #      fitness function.
  #
  # Returns:
  #   The optimized values of cc.

#  opt <- tryCatch(
#    expr=trust(objfun=SofQ.fitness, parinit=cc, rinit=1e-2, rmax=1e6,
#      blather=FALSE,
#      Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
#      eta=eta, Phi=Phi, J=J, mu=mu),
#    error=function(e) {
#      return (list(argument=NA))  # Just forget it for now; none of these alternate methods are working!
#      warning(e)
#      cat("Trying optim instead!\n")
#      opt <- optim(par=cc, fn=SofQ.optim.fn, gr=SofQ.optim.gr, method="CG",
#        Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
#        eta=eta, Phi=Phi, J=J, mu=mu)
#      opt$argument <- opt$par
#      return (opt)
#    })

  opt <- trust(objfun=SofQ.fitness, parinit=cc, rinit=1e-2, rmax=1e6,
    blather=FALSE, Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr, sigma=sigma,
    lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu)
  return (opt)
}

Plot.S.and.G <- function(Q, F.exp, xi, cc, Phi, sigma.current, sigma) {
  Layout <- grid.layout(nrow=1, ncol=2,
    widths=unit(c(2, 1), rep("null", 2)))
  LayoutNewGridPage(Layout=Layout)

  # Plot S(Q) in the left-hand region
  data.F <- data.frame(Q=Q, F.exp=F.exp, F.plus=F.exp + sigma.current,
    B=Phi %*% cc)
  melt.F <- melt(data.F, id.vars="Q")
  sig.current <- min(sigma.current)
  sig.true <- min(sigma)
  plot.F <- (ggplot(data=melt.F, aes(x=Q, y=value, colour=variable))
    + geom_line()
    + opts(title=sprintf("S(Q) for sigma = %f (cmp. to %f)", sig.current, sig.true))
    + scale_colour_manual(values=c(F.exp="black", F.plus="blue", B="red"))
    )
  print(plot.F, vp=Subplot(1, 1))

  # Plot G(r) in the right-hand region
  dr <- 0.01    # Angstrom
  r.max <- 0.9  # Angstrom
  r <- seq(from=dr, to=r.max, by=dr)
  rho.0 <- 0.0926573  # Calculated by hand for Laves paper
  K.G <- NoiseCovarianceMatrixGFromS(r=r, Q=Q, sigma=sigma)
  G.err <- sqrt(diag(K.G))
  plot.bounds <- 3 * G.err
  G.true <- -4 * pi * r * rho.0
  y.range <- range(c(G.true + plot.bounds, G.true - plot.bounds))
  M <- FourierSineTransformMatrixQR(Q=Q, r=r)
  data.G <- data.frame(r=r, G.true=G.true, G.plus=(G.true+G.err),
    G.minus=(G.true - G.err), G.FT=M %*% (F.exp - Phi %*% cc))
  melt.G <- melt(data.G, id.vars="r")
  plot.G <- (ggplot(data=melt.G, aes(x=r, y=value, colour=variable))
    + geom_line()
    + opts(title=sprintf("G(r) for sigma = %f (cmp. to %f)", sig.current, sig.true))
    + scale_colour_manual(values=c(G.true="black", G.FT="red", G.plus="blue", G.minus="blue"))
    + ylim(y.range)
    )
  print(plot.G, vp=Subplot(1, 2))
}

SofQ.OptimizeCCShrinkingSigma <- function(Q, F.exp, xi, cc=(0 * xi), p.bgr,
  sigma, lambda, eta, Phi=spline.matrix(x=Q, spline.x=xi), J, mu) {
  # Optimizes cc (i.e., the spline values; 'c' in Fischer et al.) by gradually
  # decreasing the experimental noise, from lambda down to its true value of
  # sigma.
  #
  # Args:
  #   Q:  The Q-values of the measured datapoints.
  #   F.exp:  The y-values of the measured datapoints.
  #   xi:  The spline knot x-values.
  #   cc:  Initial guess for the spline knot y-values ('c' in Fischer et al.).
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   Phi:  The spline matrix taking cc into the background curve
  #   J:  (numeric matrix) Matrix for the quadratic contribution of cc to the
  #      fitness function.
  #   mu:  (numeric vector) Vector for the linear contribution of cc to the
  #      fitness function.
  #
  # Returns:
  #   A numeric vector containing optimal values for cc.
  sigma.current <- pmax(sigma, max(lambda, sigma))
  i <- 0
  require("Cairo")
  lower.sigma.by <- 1
  repeat {
    sigma.prev <- sigma.current
    opt <- SofQ.OptimizeCC(Q=Q, F.exp=F.exp, xi=xi, cc=cc, p.bgr=p.bgr,
      sigma=sigma.current, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu)
    cc <- opt$argument
    if (length(which(sigma.current > sigma)) <= 0) {
      break
    }
    #CairoPNG(filename=sprintf("S-and-G_optimizing-c_%03d.png", i),
    #  width=1000, height=600)
    #Plot.S.and.G(Q=Q, F.exp=F.exp, xi=xi, cc=cc, Phi=Phi,
    #  sigma.current=sigma.current, sigma=sigma)
    sigma.current <- pmax(sigma.prev / (1 + lower.sigma.by), sigma)
    #dev.off()
    i <- i + 1
  }
  return (cc)
}

SofQ.OptimizeXi <- function(Q, F.exp, xi, min.dx, p.bgr, sigma, lambda, eta,
  K.G.inv, M, z, temp.init=500, reduce.new.best=0.95, reduce=0.995,
  c.at.1=1000) {
  # Optimizes xi (i.e., the spline locations) by simulated annealing Markov
  # Chain Monte Carlo (MCMC).
  #
  # Args:
  #   Q:  The Q-values of the measured datapoints.
  #   F.exp:  The y-values of the measured datapoints.
  #   xi:  The starting spline knot x-values.
  #   p.bgr:  The probability that a single pixel contains "only" background.
  #   sigma:  The Gaussian noise level for background-only points
  #   lambda:  Expected value (over background) of signal-containing pixels
  #   eta:  A limit to the straightness sensitivity
  #   K.G.inv:  The inverse of the noise covariance matrix for low-r G(r)
  #   M:  Fourier transform matrix taking F(Q) into low-r G(r)
  #   z:  Expected Fourier transform of the true B(Q) at low-r
  #   temp.init:  Initial "temperature" for simulated annealing.
  #   reduce.new.best:  The factor to reduce the temperature when we find a new
  #      best value of psi
  #   reduce:  The factor to reduce the temperature when we accept a move which
  #      is *not* the best-so-far.
  #   c.at.1:  The (c)oncentration at temperature=1, i.e., a number which
  #      describes how closely the new values should stick to the present ones
  #      (higher means tighter; any positive number will do).
  #   J:  (numeric matrix) Matrix for the quadratic contribution of cc to the
  #      fitness function.
  #   mu:  (numeric vector) Vector for the linear contribution of cc to the
  #      fitness function.
  #
  # Returns:
  #   The optimal xi values.
  ComputeJAndMu <- function() {
    assign("Phi", spline.matrix(x=Q, spline.x=xi), envir=parent.frame())
    M.Phi <- M %*% Phi
    K.G.inv.M.Phi <- K.G.inv %*% M.Phi
    assign("mu", t(K.G.inv.M.Phi) %*% z, envir=parent.frame())
    assign("J", t(K.G.inv.M.Phi) %*% M.Phi, envir=parent.frame())
  }
  temp <- temp.init
  ComputeJAndMu()
  time.start <- proc.time()["elapsed"]
  cc <- SofQ.OptimizeCCShrinkingSigma(Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu)
  fit.best <- SofQ.fitness(params=cc, Q=Q, F.exp=F.exp, xi=xi, p.bgr=p.bgr,
    sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu,
    calc.grad=FALSE, calc.hess=FALSE)
  xi.traces <- data.frame(
    elapsed=(proc.time()["elapsed"] - time.start),
    xi=matrix(xi, nrow=1), cc=matrix(cc, nrow=1),
    log.p=0, imp=0, proposal.prob=0,
    temp=temp, psi=fit.best$value)
  fit.new <- fit.best
  while (temp > 1) {
    # Draw new xi-values (and calculate proposal probability ratio!)
    time.start <- proc.time()["elapsed"]
    xi.new <- DrawNewXi(xi=xi, min.dx=min.dx, concentration=(c.at.1 / temp))
    ComputeJAndMu()
    # Optimize c-values for the new xi
    cc <- SofQ.OptimizeCCShrinkingSigma(Q=Q, F.exp=F.exp, xi=xi.new$xi,
      p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu)
    fit.prev <- fit.new
    fit.new <- SofQ.fitness(params=cc, Q=Q, F.exp=F.exp, xi=xi.new$xi,
      p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu,
      calc.grad=FALSE, calc.hess=FALSE)
    # Calculate acceptance probability
    improvement <- (fit.new$value - fit.prev$value)
    log.prob <- (improvement / temp) + xi.new$d.log.p
    accepted <- (exp(log.prob) > runif(1))
    # Set new value of xi and adjust the temperature
    temp.prev <- temp
    if (accepted) {
      xi <- xi.new$xi
      if (fit.new$value < fit.best$value) {
        fit.best <- fit.new
        temp <- temp * reduce.new.best
      } else {
        temp <- temp * reduce
      }
    }
    xi.traces <- rbind(deparse.level=0, xi.traces, data.frame(
        elapsed=(proc.time()["elapsed"] - time.start),
        xi=matrix(xi.new$xi, nrow=1), cc=matrix(cc, nrow=1),
        log.p=log.prob, imp=improvement, proposal.prob=xi.new$d.log.p,
        temp=temp.prev, psi=fit.new$value))
    save(file="xi_traces.ROBJ", xi.traces, Q, F.exp)
  }
  return (xi.traces)
}

SofQ.QuenchedXiOptimization <- function(Q, F.exp, xi, min.dx, p.bgr, sigma,
  lambda, eta, z, K.G.inv, M=M, temp.init=500, reduce.new.best=0.95,
  reduce=0.995, c.at.1=1000)
{
  xi.traces <- SofQ.OptimizeXi(Q=Q, F.exp=F.exp, xi=xi, min.dx=min.dx,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, K.G.inv=K.G.inv, M=M,
    z=z, temp.init=temp.init, reduce.new.best=reduce.new.best, reduce=reduce,
    c.at.1=c.at.1)
  iter <- 1
  all.xi.traces <- data.frame()
  repeat {
    xi.traces$iter <- iter
    all.xi.traces <- rbind(deparse.level=0, all.xi.traces, xi.traces)
    save(file="all_xi_traces.ROBJ", all.xi.traces, Q, F.exp)
    i.best <- max(which(all.xi.traces$psi == min(all.xi.traces$psi)))
    temp.quench <- all.xi.traces$temp[i.best] * 0.5
    if (temp.quench < 1) {
      break;
    }
    iter <- iter + 1
    i.xi <- grep("^xi\\.\\d+$", names(all.xi.traces))
    xi <- unname(unlist(all.xi.traces[i.best, i.xi]))
    xi.traces <- SofQ.OptimizeXi(Q=Q, F.exp=F.exp, xi=xi, min.dx=min.dx,
      p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, K.G.inv=K.G.inv, M=M,
      z=z, temp.init=temp.quench, reduce.new.best=reduce.new.best,
      reduce=reduce, c.at.1=c.at.1)
  }
  return (all.xi.traces)
}

SofQ.TestOptimizingCC <- function(laves.data, rho.0) {
  sigma <- 0.0136
  lambda <- 1
  p.bgr <- 0.5
  dx.min <- 3
  Q <- laves.data$Q
  N <- length(Q)
  xi <- c(Q[1], 4.5, 7.5, 11, 14.5, 21, Q[N])
  E <- length(xi)
  cc <- rnorm(n=E)
  S <- laves.data$S.Q.raw
  i.high <- which(Q > 17 & Q < 23)
  S.base <- mean(S[i.high])
  sigma <- 0.0136 * sqrt(S / S.base)  # Approximating "scaled Poisson"
  eta <- 0.0001542393  # precomputed
  F.exp <- Q * (S - 1)
  dr <- 0.01    # Angstrom
  r.max <- 0.9  # Angstrom
  r <- seq(from=dr, to=r.max, by=dr)
  Phi <- spline.matrix(x=Q, spline.x=xi)
  M.Phi <- FourierSineTransformMatrixQR(Q=Q, r=r) %*% Phi
  K.G <- NoiseCovarianceMatrixGFromS(r=r, Q=Q, sigma=sigma)
  diag(K.G) <- diag(K.G) + abs(min(eigen(K.G)$values)) * 1e4  # avoid singularity
  K.G.inv.M.Phi <- solve(K.G) %*% M.Phi
  z <- FourierSineTransformQR(f.Q=F.exp, Q=Q, r=r) + 4 * pi * rho.0 * r
  mu <- t(K.G.inv.M.Phi) %*% z
  J <- t(K.G.inv.M.Phi) %*% M.Phi
  #mu <- 0 * mu; J <- 0 * J
  cc.best <- SofQ.OptimizeCCShrinkingSigma(Q=Q, F.exp=F.exp, xi=xi, cc=cc,
    p.bgr=p.bgr, sigma=sigma * Q, lambda=lambda, eta=eta, Phi=Phi, J=J, mu=mu)
  return(cc.best)
}

SofQ.TestOptimizingXi <- function(laves.data, rho.0=0.0926573) {
  lambda <- 1
  p.bgr <- 0.5
  dx.min <- 3
  Q <- laves.data$Q
  N <- length(Q)
  xi <- c(Q[1], 4.5, 7.6, 11, 14.5, 21, Q[N])
  E <- length(xi)
  cc <- rnorm(n=E)
  S <- laves.data$S.Q.raw
  i.high <- which(Q > 17 & Q < 23)
  S.base <- mean(S[i.high])
  sigma <- 0.0136 * sqrt(S / S.base)  # Approximating "scaled Poisson"
  #eta <- Eta(sigma=sigma, dx.min=dx.min, dx.total=diff(range(xi)))
  eta <- 0.0001542393  # precomputed from above
  F.exp <- Q * (S - 1)
  dr <- 0.01    # Angstrom
  r.max <- 0.9  # Angstrom
  r <- seq(from=dr, to=r.max, by=dr)
  z <- FourierSineTransformQR(f.Q=F.exp, Q=Q, r=r) + 4 * pi * rho.0 * r
  M <- FourierSineTransformMatrixQR(Q=Q, r=r)
  K.G <- NoiseCovarianceMatrixGFromS(r=r, Q=Q, sigma=sigma)
  diag(K.G) <- diag(K.G) + abs(min(eigen(K.G)$values)) * 1e4  # avoid singularity
  K.G.inv <- solve(K.G)
  xi.traces <- SofQ.QuenchedXiOptimization(Q=Q, F.exp=F.exp, xi=xi, min.dx=dx.min, p.bgr=p.bgr,
    sigma=sigma * Q, lambda=lambda, eta=eta, z=z, K.G.inv=K.G.inv, M=M)
  return(xi.traces)
}

##############################################################################
# MAIN LOGIC FUNCTIONS                                                       #
##############################################################################

NoisyTestFunction <- function(x, min.dx=1.0, w.range=c(0.1, 0.3), n.knots=6,
  n.peaks=5, max.peak.height=10, sigma.noise=0.5) {
  # Generates a random function with smooth background, returning separate
  # curves for the background, the peaks, and the noise.
  #
  # Args:
  #   x:  numeric vector of x-values
  #   min.dx:  minimum spacing between consecutive spline knots
  #   w.range:  range in peak widths
  #   n.knots:  number of spline knots for background (less means smoother)
  #   n.peaks:  number of peaks to add
  #   max.peak.height:  Maximum height of any peak
  #   sigma.noise:  noise level
  #
  # Returns:
  #   list() with the following elements:
  #     $curves:  full-res curves (data.frame() with following columns):
  #       $x:  x-values for this function (same as this function received)
  #       $total:  total curve (bgr + peaks + noise)
  #       $bgr:  background contribution
  #       $peaks:  peaks that get added to bgr
  #       $noise:  IID Gaussian noise
  #     $knots:  spline knots (data.frame() with following columns):
  #       $x:  x-values of knots
  #       $y:  y-values of knots

  # Calculate background signal
  spline.wiggle.room <- WiggleRoom(E=n.knots, min.dx=min.dx,
    x.range=diff(range(x)))
  space.fractions <- rgamma(n=(n.knots - 1), shape=1)  # arbitrary
  space.fractions <- space.fractions / sum(space.fractions)
  cumulative.extra.dist <- c(0, cumsum(space.fractions) * spline.wiggle.room)
  knots.x <- x[1] + ((1:n.knots) - 1) * min.dx + cumulative.extra.dist
  knots.y <- rnorm(n=n.knots)
  bgr <- spline(x=knots.x, y=knots.y, xout=x)$y
  # Calculate random peaks and noise
  amplitude <- runif(n=n.peaks, min=0, max=max.peak.height)
  width <- runif(n=n.peaks, min=w.range[1], max=w.range[2])
  buf <- 2 * w.range[2]
  centre <- runif(n=n.peaks, min=min(x) + buf, max=max(x) - buf)
  separate.peaks <- mapply(
    FUN=function(amp, wid, cen) amp * exp (-0.5 * ((x - cen) / wid) ^ 2),
    amp=amplitude, wid=width, cen=centre)
  peaks <- rowSums(separate.peaks)
  noise <- rnorm(sd=sigma.noise, n=length(x))
  total <- bgr + peaks + noise
  # Bundle these curves into data structures and return them:
  curves <- data.frame(x=x, total=total, bgr=bgr, peaks=peaks, noise=noise)
  knots <- data.frame(x=knots.x, y=knots.y)
  return (list(curves=curves, knots=knots))
}

PlotNoisyTestFunction <- function(x, min.dx=1.0, w.range=c(0.1, 0.3),
  n.knots=6, n.peaks=5, max.peak.height=10, sigma.noise=0.5) {
  NTF <- NoisyTestFunction(x=x, min.dx=min.dx, w.range=w.range,
    n.knots=n.knots, n.peaks=n.peaks, max.peak.height=max.peak.height,
    sigma.noise=sigma.noise)
  melted.curves <- melt(NTF$curves, id.vars="x")
  p <- (ggplot(data=melted.curves, aes(x=x))
    + geom_line(aes(y=value, colour=variable, size=variable))
    + geom_point(data=NTF$knots, aes(y=y))
    + scale_colour_manual(values=c(
        total="black", bgr="red", noise="grey", peaks="yellow"))
    + scale_size_manual(values=c(
        total=1.5, bgr=1.5, noise=0.5, peaks=0.5))
    )
  print(p)
}

TestGradient <- function() {
  sigma <- 1
  lambda <- 10
  NTF <- NoisyTestFunction(x=seq(from=0, to=10, length.out=200), sigma.noise=sigma, max.peak.height=lambda)
  xi <- NTF$knots$x
  E <- length(xi)
  cc <- rnorm(n=E)
  y <- NTF$curves$total
  p.bgr <- 0.5
  psi.base <- Psi(x=NTF$curves$x, y=y, cc=cc, p.bgr=p.bgr, xi=xi, sigma=sigma,
    lambda=lambda, calc.grad=TRUE, calc.hess=TRUE)
  d.cc <- rnorm(n=E, sd=1)
  for (i in 1:10) {
    psi.new <- Psi(x=NTF$curves$x, y=y, cc=(cc + d.cc), p.bgr=p.bgr, xi=xi,
      sigma=sigma, lambda=lambda, calc.grad=FALSE, calc.hess=FALSE)
    d.psi.actual <- psi.new$funct - psi.base$funct
    d.psi.expect <- psi.base$grad %*% d.cc
    if (!identical(psi.base$hess, NA)) {
      d.psi.expect <- d.psi.expect + 0.5 * t(d.cc) %*% psi.base$hess %*% d.cc
    }
    cat(sprintf("Ratio: %16.13f (Actual %8.6g / Expected %8.6g)\n",
        d.psi.actual / d.psi.expect, d.psi.actual, d.psi.expect))
    d.cc <- d.cc * 0.1
  }
}

ppp <- function(x) {
  return (
    ggplot(data=x[nrow(x):1, ], aes(x=i, y=ratio, colour=as.factor(iter))) 
    + geom_point()
    + scale_colour_brewer())
}

TestOptimizingCC <- function() {
  sigma <- 0.3
  lambda <- 10
  p.bgr <- 0.5
  dx.min <- 1
  NTF <- NoisyTestFunction(x=seq(from=0, to=10, length.out=200),
    sigma.noise=sigma, max.peak.height=lambda, min.dx=dx.min)
  xi <- NTF$knots$x
  E <- length(xi)
  cc <- rnorm(n=E)
  eta <- Eta(sigma=sigma, dx.min=dx.min, dx.total=diff(range(xi)))
  y <- NTF$curves$total
  cc.best <- OptimizeCCShrinkingSigma(curves=NTF$curves, cc=cc, xi=xi,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)
  psi.guessed <- PsiWrapped(params=cc.best, x=NTF$curves$x, y=NTF$curves$total, xi=xi,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)$value
  psi.best <- PsiWrapped(params=NTF$knots$y, x=NTF$curves$x, y=NTF$curves$total, xi=xi,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)$value
  cat("Fit value (psi): guessed =", psi.guessed, "and best is", psi.best, "\n")
  return(cc.best)
}

TestOptimizingXi <- function() {
  sigma <- 0.3
  lambda <- 10
  p.bgr <- 0.5
  dx.min <- 1
  NTF <- NoisyTestFunction(x=seq(from=0, to=10, length.out=200),
    sigma.noise=sigma, max.peak.height=lambda, min.dx=dx.min)
  xi <- NTF$knots$x
  xi <- seq(from=min(NTF$knots$x), to=max(NTF$knots$x), length.out=length(xi))
  E <- length(xi)
  eta <- Eta(sigma=sigma, dx.min=dx.min, dx.total=diff(range(xi)))
  y <- NTF$curves$total
  xi.traces <- QuenchedXiOptimization(curves=NTF$curves, xi=xi, min.dx=dx.min,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)
  i.best <- which(xi.traces$psi == min(xi.traces$psi))[1]
  xi.best <- unname(unlist(xi.traces[i.best, 1:E]))
  cc.best <- OptimizeCCShrinkingSigma(curves=NTF$curves, xi=xi.best,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)
  cc.true <- OptimizeCCShrinkingSigma(curves=NTF$curves, xi=NTF$knots$x,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)

  psi.best.guess <- PsiWrapped(params=cc.best, x=NTF$curves$x,
    y=NTF$curves$total, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
    eta=eta)$value
  psi.true.xi <- PsiWrapped(params=cc.true, x=NTF$curves$x, y=NTF$curves$total,
    xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta)$value
  psi.true.all <- PsiWrapped(params=NTF$knots$y, x=NTF$curves$x,
    y=NTF$curves$total, xi=xi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
    eta=eta)$value

  cat(sprintf("psi: guessed xi = %10g; true xi = %10g; true (xi, cc) = %10g\n",
      psi.best.guess, psi.true.xi, psi.true.all))
  return(xi.traces)
}

TestMultiIndexGradAndHess <- function(function.name, rand.seed=1,
  sigma=0.3, lambda=10, p.bgr=0.5) {
  # Tests that the gradient and Hessian are working correctly for one of the
  # "multi-index" (i.e., per-datapoint) quantities (i.e., f(c) or h(c))
  #
  # Args:
  #   function.name:  The function to test (either PsiFBgr or PsiHSig).
  #
  # Returns:
  #   Nothing; used for its side-effect (prints to stdout).
  set.seed(rand.seed)
  x <- seq(from=0, to=10, length.out=200)
  NTF <- NoisyTestFunction(x=x, sigma.noise=sigma, max.peak.height=lambda)
  xi <- NTF$knots$x
  E <- length(xi)
  cc <- NTF$knots$y
  cc <- rnorm(n=E)
  y <- NTF$curves$total
  Phi <- spline.matrix(x=x, spline.x=xi)
  funct <- get(function.name)
  f.base <- funct(cc=cc, d=y, Phi=Phi, p.bgr=p.bgr, sigma=sigma, lambda=lambda,
    calc.grad=TRUE, calc.hess=TRUE)
  d.cc <- rnorm(n=length(cc), sd=1)
  score.data <- data.frame()
  for (i in 1:10) {
    f.new <- funct(cc=(cc + d.cc), d=y, Phi=Phi, p.bgr=p.bgr, sigma=sigma,
      lambda=lambda, calc.grad=FALSE, calc.hess=FALSE)
    df.actual <- f.new$funct - f.base$funct
    hess.vector <- 0.5 * ReduceFrom3D(H=f.base$hess, v=d.cc, premultiply=TRUE)
    df.expect <- (f.base$grad + hess.vector) %*% d.cc
    cat(sprintf("Ratio: %16.14f (Actual %8.6g / Expected %8.6g)\n",
        (df.actual / df.expect)[75], df.actual[75], df.expect[75]))
    score <- log(df.actual / df.expect - 1) / log(10)
    score.data <- rbind(score.data, data.frame(score=score, iter=i))
    d.cc <- d.cc * 0.1
  }
  stat_sum_single <- function(fun, geom="point", ...) { 
    # From http://had.co.nz/ggplot2/stat_summary.html
    stat_summary(fun.y=fun, colour="red", geom=geom, size = 1, ...) 
  }
  p <- (ggplot(data=score.data, aes(x=iter, y=score))
    + stat_boxplot(aes(group=iter))
    + stat_sum_single(mean, geom="line")
    + opts(title=sprintf("Testing Hessian for function '%s'", function.name))
    + scale_x_continuous("Iteration number (div. dc by 10 each time)")
    + scale_y_continuous("Log10 (actual / expected - 1)")
    #+ geom_point()
    )
  print(p)
  return (p)
}

EvaluateEtaUniversal <- function(w) {
  eta.vals <- NA * w
  for (i in 1:length(w)) {
    eta.vals[i] <- EtaUniversal(width=w[i])
    cat(sprintf("eta(%8.1f) = %12.3f\n", w[i], eta.vals[i]))
  }
  return (eta.vals)
}

TryFittingCalculatedLaves <- function(laves.data, sigma=0.001, lambda=0.327,
  dx.min=5, p.bgr=0.5) {
  xi <- c(min(laves.data$Q), 6, 12, 17.1, max(laves.data$Q))
  E <- length(xi)
  curves <- data.frame(x=laves.data$Q, total=laves.data$S.Q.calc)
  eta <- Eta(sigma=sigma, dx.min=dx.min, dx.total=diff(range(curves$x)))
  opt <- QuenchedXiOptimization(curves=curves, xi=xi, min.dx=dx.min,
    p.bgr=p.bgr, sigma=sigma, lambda=lambda, eta=eta, c.at.1=1e4)
  return (opt)
}

laves.data <- read.table("ZrTiNi_C14-Laves.txt", sep="\t", header=TRUE)
