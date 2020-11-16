#' Prediction of Marginal Effects
#'
#' Given a fitted model by the rcbr NPMLE procedure predictions are made
#' at new design points given by the \code{newdata} argument.
#'
#' @param  object is the fitted NPMLE object
#' @param  \dots  is expected to contain an argument \code{newdata}
#' @return a vector pf predicted probabilities 
#' @seealso \code{bound.KW2} for a prediction function with bounds
#' @export
predict.KW2 <- function(object, ...){
    uv <- object$uv
    W <- object$W
    mu <- 40
    mv <- 40
    eps <- 0.0001
    dots <- list(...)
    if(!length(dots$newdata))
	stop("No newdata to predict at.")
    nd <- as.matrix(dots$newdata)
    if(length(dots$smooth) && (dots$smooth > 0)){
       ru <- range(uv[W > eps,1])
       rv <- range(uv[W > eps,2])
       s <- dots$smooth
       u <- seq(ru[1] - 3 * s, ru[2] + 3 * s, length = mu)
       v <- seq(rv[1] - 3 * s, rv[2] + 3 * s, length = mv)
       fs <- matrix(0, mu, mv)
       for(i in 1:mu){
          for(j in 1:mv){
             Z <- cbind(u[i]-uv[,1],v[j]-uv[,2])
             fs[i,j] <- sum(mvtnorm::dmvnorm(Z, sigma = diag(2) * s) * W)
             }
        }
        uv <- expand.grid(u,v)
        W <- c(fs)
        W <- W/sum(W)
    }
    X <- cbind(1, nd)
    B <- as.matrix(cbind(uv, 1))
    B <- B/sqrt(apply(B^2, 1, sum)) # is this necessary?
    A <- X %*% t(B)
    apply(A, 1, function(a) sum(W[a > 0]))
}

