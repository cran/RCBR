#' log likelihood for Gautier Kitamura procedure
#' @param object a fitted object of class "GK"
#' @param ... other parameters for logLik
#' @return a scalar log likelihood
#' @export
#'
logLik.GK <- function(object, ...) {
    w <- object$w/sum(object$w)
    uv <- expand.grid(object$u, object$v)
    B <- as.matrix(cbind(uv, 1))
    Xp <- object$X[object$y==1,]
    Xn <- object$X[object$y==0,]
    pos <- apply(Xp %*% t(B), 1, function(a) sum(w[a > 0]))
    neg <- 1 - apply(Xn %*% t(B), 1, function(a) sum(w[a > 0]))
    sum(log(c(pos, neg)))
}

