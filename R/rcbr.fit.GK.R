#' Gautier and Kitamura (2013) bivariate random coefficient binary response
#'
#' This is an implementation based on the matlab version of Gautier and
#' Kitamura's deconvolution method for the bivariate random coefficient 
#' binary response model.  Methods based on the fitted object are provided
#' for \code{predict}, \code{logLik} and \code{plot}.requires orthopolynom
#' package for Gegenbauer polynomials
#'
#' @param  X the design matrix expected to have an intercept column of
#'	ones as the first column.
#' @param  y the binary response.
#' @param  control is a list of tuning parameters for the fitting,see
#'	\code{GK.control} for further details.
#' @return a list with components:
#' \describe{
#'   \item{u}{grid values}
#'   \item{v}{grid values}
#'   \item{w}{estimated function values on 2d u x v grid}
#'   \item{X}{design matrix}
#'   \item{y}{response vector} 
#' }
#' @author Gautier and Kitamura for original matlab version, Jiaying Gu
#' and Roger Koenker for the R translation.
#'
#' @references
#' Gautier, E. and Y. Kitamura (2013)  Nonparametric estimation in random coefficients
#' binary choice models, \emph{Ecoonmetrica}, 81, 581-607.
#'
#' @keywords  nonparametrics  
#' @export
rcbr.fit.GK <- function(X, y, control){
    u <- control$u; v <- control$v
    T <- control$T; TX <- control$TX
    Mn <- control$Mn; tol <- control$tol
    du <- diff(u)[1]
    dv <- diff(v)[1]
    if(!all(abs(diff(u) - du) < tol)) stop("u-grid must be equally spaced")
    if(!all(abs(diff(v) - dv) < tol)) stop("v-grid must be equally spaced")
    stopifnot(requireNamespace("orthopolynom"))
    V <- function(d) (2 * pi^(d/2))/gamma(d/2) # Surface area of the (d-1)-sphere
    H <- function(n,d) ((2*n+d-2)*choose(n+d-2,d-2)/(n+d-2))
    lambda <- function(n, d) { #See Proposition A.4 in GK.
	p <- (n-1)/2 
	if(n == 1) V(d-1)/(d-1)
	else (-1)^p * V(d-1) * prod(seq(1,2*p-1,2))/prod(seq(d-1, d+2*p-1,2))
    }
    Gegen <- function(T,d) # yields T+1 functions
       lapply(orthopolynom::gegenbauer.polynomials(T,(d-2)/2), "as.function") 
    Chi <- function(n, d, T, l = 3, s = 3) # See KG Prop A.3, and Section 6
       (1 - (Zeta(n,d)/(Zeta(T,d) + 1))^(s/2))^l 
    Zeta <- function(n,d) n * (n+d-2) # See GK Lemma A.1
    fX <- function(x, X, TX) {
	n <- nrow(X)
	d <- ncol(X)
	G <- Gegen(TX, d)
	S <- sapply(1:(TX+1), function(i) mean(G[[i]](X %*% x)))
	s <- sapply(1:(TX+1), function(i) G[[i]](1))
	R <- Chi(1:(TX+1),d,TX+1) * H(0:TX,d)/s
	max(sum(R * S)/V(d),0)
    }
    X <- as.matrix(X/sqrt(apply(X^2,1,sum)))
    n <- nrow(X)
    d <- ncol(X)
    B <- cbind(expand.grid(u,v), 1)
    Bsd <- sqrt(apply(B^2,1,sum))
    B <- as.matrix(B/Bsd)
    G <- Gegen(2*T + 2, d)
    R <- sapply(0:(T-1), function(p) Chi(2*p+1,d,2*T+1) * H(2*p+1,d)/
       (lambda(2*p+1,d)*G[[2*p+2]](1)))
    fx <- pmax(sapply(1:n, function(i) fX(X[i,],X,TX)), Mn)
    w <- rep(0, nrow(B))
    for(j in 1:nrow(B)){
	S <- sapply(0:(T-1), function(p) mean((2 * y - 1) * G[[2*p + 2]](X %*% B[j,])/fx))
	w[j] <- max(2 * sum(R * S)/V(d), 0)
	}
    w <- du * dv * w/Bsd^3
    z <- list(u = u, v = v, w = w, X = X, y = y)
    class(z) <- "GK"
    return(z)
}
