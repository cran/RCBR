#' Plot a KW2 object
#'
#' Given a fitted model by the rcbr NPMLE procedure plot the estimated mass points
#'
#' @param  x is the fitted NPMLE object
#' @param  smooth is a parameter to control bandwidth of the smoothing if a contour plot
#'	of the estimated density is desired, default is no smoothing and only
#'	the mass points of the discrete estimate are plotted.
#' @param  pal a color palette 
#' @param  inches as used in \code{symbols} to control size of mass points
#' @param  N scaling of the color palette
#' @param  tol tolerance for size of mass points
#' @param  ... other arguments to pass to \code{symbols}, notably e.g. \code{add = TRUE}
#' @return nothing (invisibly)
#' @importFrom grDevices heat.colors
#' @importFrom graphics contour
#' @export
plot.KW2 <- function(x, smooth = 0, pal = NULL, inches = 1/6, N = 25, tol = 0.001, ...){
    W <- x$W
    s <- (W > tol)
    D <- data.frame(w = W[s], x = x$uv[s,1], y = x$uv[s,2])
    if(smooth > 0){
	uv <- x$uv
        mu <- mv <- 40
        ru <- range(uv[s, 1])
        rv <- range(uv[s, 2])
        u <- seq(ru[1] - 4 * smooth, ru[2] + 4 * smooth, length = mu)
        v <- seq(rv[1] - 4 * smooth, rv[2] + 4 * smooth, length = mv)
        fs <- matrix(0, mu, mv)
        for (i in 1:mu) {
            for (j in 1:mv) {
                Z <- cbind(u[i] - uv[, 1], v[j] - uv[, 2])
                fs[i, j] <- sum(mvtnorm::dmvnorm(Z, sigma = diag(2) * smooth) * W)
            }
        }
        contour(u, v, matrix(fs,mu,mv), ...)
    }
    else{
	if(!length(pal))
	    pal <- heat.colors(N, alpha = .5)
	with(D, {
	 bg <- round(1 + 10 * (1 - w/max(w)))
	 symbols(x = x, y = y, circles = sqrt(w), inches = inches,  
		    fg = "grey30", bg = bg, ...)
	 palette(pal)
	})
    }
}

