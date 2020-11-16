#'
#' NPMLE fitting for random coefficient binary response model
#'
#' Exact NPMLE fitting requires that the \code{uv} argument contain a matrix
#' whose rows represent points in the interior of the locally maximal polytopes
#' determined by the hyperplane arrangement of the observations.  If it is not
#' provided it will be computed afresh here; since this can be somewhat time
#' consuming, \code{uv} is included in the returned object so that it can be
#' reused if desired.  Approximate NPMLE fitting can be achieved by specifying
#' an equally spaced grid of points at which the NPMLE can assign mass using
#' the arguments \code{u} and \code{v}.  If the design matrix \code{X} contains
#' only 2 columns, so we have the Cosslett, aka current status, model then the
#' polygons in the prior description collapse to intervals and the default method
#' computes the locally maximal count intervals and passes their interior points
#' to the optimizer of the log likelihood.  Alternatively, as in the bivariate
#' case one can specify a grid to obtain an approximate solution.
#'
#' @param  x the design matrix expected to have an intercept column of
#' ones as the first column, the last column is presumed to contain values of
#' the covariate that is designated to have coefficient one.
#' @param  y the binary response.
#' @param  control is a list of parameters for the fitting, see
#' \code{KW.control} for further details.
#' @return a list with components:
#' \itemize{
#'   \item  uv evaluation points for the fitted distribution
#'   \item  W estimated mass associated with the \code{uv} points
#'   \item  logLik the loglikelihood value of the fit
#'   \item  status mosek solution status
#' }
#' @author  Jiaying Gu and Roger Koenker 
#' @references
#' Gu, J. and R. Koenker (2018)  Nonparametric maximum likelihood estimation 
#' of the random coefficients binary choice model, preprint.
#' @keywords  nonparametrics  
#' @export
rcbr.fit.KW2 <- function (x, y, control){
    uv <- control$uv
    u <- control$u
    v <- control$v
    initial <- control$initial
    epsbound <- control$epsbound
    epstol <- control$epstol
    presolve <- control$presolve
    verb <- control$verb
    if(length(u) * length(v)){
	du <- diff(u)[1]
	dv <- diff(v)[1]
	if(!all(abs(diff(u) - du) < epstol)) stop("u-grid must be equally spaced")
	if(!all(abs(diff(v) - dv) < epstol)) stop("v-grid must be equally spaced")
	uv <- expand.grid(u, v)
	cells <- NULL
	locmax <- NULL
    }
    else if(!length(uv)){ # Find locally maximal polygons
	cells <- NICER(x[,1:2], -x[,3], initial, verb, epsbound, epstol) 
	sv <- cells$SignVector
	sv[which(y==0),] <- -sv[which(y==0),]
	locmax <- neighbours(sv)
	uv <- t(cells$w[,locmax])
    }
    n <- NROW(x)
    p <- NCOL(x)
    w <- rep(1, n)/n
    B <- as.matrix(cbind(uv, 1))
    d <- rep(1, length(nrow(uv)))
    A <- x %*% t(B)
    A <- 1 * (A >= 0)
    A <- (y == 1) * A + (y == 0) * (1 - A)
    f <- KWDual(A, d, w, verb = verb)  # FIXME:  pass other control parameters!!!
    W <- f$f
    x1 <- x[y==1,]
    x0 <- x[y==0,]
    pos <- apply(x1 %*% t(B), 1, function(a) sum(W[a > 0]))
    neg <- 1 - apply(x0 %*% t(B), 1, function(a) sum(W[a > 0]))
    logLik <- sum(log(c(pos, neg)))
    z <- list(uv = uv, W = W,  logLik = logLik, cells = cells, locmax = locmax, x = x, status = f$status)
    class(z) <- c("KW2", "density")
    return(z)
}
