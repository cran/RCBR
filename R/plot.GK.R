#' Plot a GK object
#'
#' Given a fitted model by the Guatier-Kitamura procedure plot the estimated density contours
#'
#' @param  x is the fitted GK object
#' @param  ... other arguments to pass to \code{contour}, notably e.g. \code{add = TRUE}
#' @return nothing (invisibly)
#' @importFrom graphics contour
#' @export
plot.GK <- function(x, ...){
    contour(x$u, x$v, matrix(x$w,length(x$u), length(x$v)), ...)
}

