#!/usr/bin/R

# Create the datafile for F(Q) = Q[S(Q) - 1],
# for simulated bulk gold with sigma = 0.01 nm

# Charles R. Hogg III
# 2012-07-20

# Functions to calculate F(Q)
F.per.nm <- function(Q, G, r) {
  if (!require("gppois")) {
    stop("We need the gppois package for the Widths function")
  }
  dr <- Widths(r)
  return (sum(G * sin(Q * r) * dr))
}
F.Q <- function(Q, G, r) {
  return (sapply(X=Q, FUN=F.per.nm, G=G, r=r))
}

###############################################################################
# MAIN LOGIC

# For which Q-values should we compute F(Q)?  (I'm using the same ones from a
# RMC output file which Igor gave me, except I'm converting to units of 1/nm
# instead of 1/A.)
Q <- seq(from=5.1, to=299.9, by=0.1)

# Read in the original G(r) data
g.data <- read.table(file="Au_g-G_double_rmax-78nm.dat", sep="\t", header=TRUE)

# If we simply Fourier-transform this, we'll get "termination ripples", so
# first convolve with a Gaussian envelope
sigma.r <- 25  # nm
g.data$G.damped <- g.data$G.of.r * exp(-(g.data$r.nm / sigma.r) ^ 2)

# Now we can do the FT.  (I wonder how adaptive quadrature compares to the
# "dumb" method, in terms of speed and accuracy...?)
Au.F.Q <- F.Q(Q=Q, G=g.data$G.damped, r=g.data$r.nm)
bulkAuSig0p01nm <- data.frame(Q=Q, F.Q=Au.F.Q, S.Q=(1 + Au.F.Q / Q))
save(file="../../data/bulkAuSig0p01nm.rda", bulkAuSig0p01nm)
