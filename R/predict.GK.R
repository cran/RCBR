#' Prediction of Marginal Effects
#'
#' Given a fitted model by the Gautier Kitamura procedure predictions are made
#' at new design points given by the \code{newdata} argument.
#'
#' @param  object is the fitted object of class "GK"
#' @param  \dots  is expected to contain an argument \code{newdata}
#' @return a vector pf predicted probabilities 
#' @export
predict.GK <- function(object, ...){
    # Add smoothing option after further testing
    dots <- list(...)
    if(!length(dots$newdata))
	stop("No newdata to predict at.")
    X <- cbind(1, as.matrix(dots$newdata, ncol = 2))
    B <- as.matrix(cbind(expand.grid(object$u,object$v), 1))
    B <- B/sqrt(apply(B^2, 1, sum)) # is this necessary?
    A <- X %*% t(B)
    apply(A, 1, function(a) sum(object$w[a > 0]))
}
