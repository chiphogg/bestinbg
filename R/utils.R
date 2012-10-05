#' Inverts the result of the 'order' function.
#'
#' @param i  Numeric vector with a permutation of the integers from
#'     \code{1:length(i)}.
#'
#' @return Numeric vector \code{r} such that \code{r[i] == 1:length(i)}.
invert.order <- function(i) {
  r <- i + NA
  for (j in 1:length(i)) {
    r[i[j]] <- j
  }
  return (r)
}

#' Compute the width for each x (specifically, its Voronoi cell size) to aid
#' in numerical integration.
#'
#' @param x  A sorted numeric vector of x-values
Dx <- function(x) {
  n <- length(x)
  i <- order(x)
  x.sort <- x[i]
  dx <- diff(c(x.sort[1], 0.5 * (x.sort[-1] + x.sort[-n]), x.sort[n]))
  return (dx[invert.order(i)])
}

#' Setup a new grid page with the specified Layout.  
#' (Adapted from http://bit.ly/9m4zyD)
#'
#' @param Layout  Result of calling grid.layout function with the desired
#'      arrangement.
LayoutNewGridPage <- function(Layout, ...) {
  grid.newpage()
  pushViewport(viewport(layout = Layout))
}

#' Terse wrapper to return a viewport into the specified row (x) and 
#' column (y) of a viewport.
#' (Adapted from http://bit.ly/9m4zyD)
#'
#' @param x  Numeric; the row(s) of the layout.
#' @param y  Numeric; the column(s) of the layout.
Subplot <- function(x, y) {
  viewport(layout.pos.row=x, layout.pos.col=y)
}

#' Plots noisy data, the true background curve, and an estimate for the
#' background curve.
#'
#' @param curves  Something like the output of NoisyTestFunction()$curves: a
#'      data.frame having columns $x, $total, and $bgr (as well as $peaks and
#'      $noise).
#' @param xi  The spline knot x-values.
#' @param cc  The spline knot y-values ('c' in Fischer et al.).
#' @param label  The title for the plot.
#'
#' @return a \code{ggplot()} object with the desired plot.
PlotBackgroundGuess <- function(curves, xi, cc, label="Background guess") {
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

