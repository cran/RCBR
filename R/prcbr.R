#' Profiling estimation methods for RCBR models
#' 
#' Profile likelihood and (GEE) score methods for estimation of random coefficient binary
#' response models.  This function is a wrapper for \code{rcbr} that uses the offset
#' argument to implement estimation of additional fixed parameters.  It may be useful to
#' restrict the domain of the optimization over the profiled parameters, this can be
#' accomplished, at least for box constraints by setting \code{omethod = "L-BFGS-B"}
#' and specifying the \code{lo} and \code{up} accordingly.
#'
#' @param formula is of the extended form enabled by the \pkg{Formula} package.
#' In the Cosslett, or current status, model the formula takes the form 
#' \code{y ~ v | z} where \code{v} is the covariate designated to have coefficient
#' one, and \code{z} is another covariate or group of covariates that are assumed
#' fixed coefficients that are to be estimated.   
#' @param b0 is either an initial value of the parameter for the Z covariates
#' or a matrix of such values, in which case optimization occurs over this discrete
#' set, when there is only one covariate then b0 is either scalar, or a vector.
#' @param data data frame for formula variables
#' @param logL if logL is TRUE the log likelihood is optimized, otherwise
#' a GEE score criterion is minimized.
#' @param omethod optimization method for \code{optim}, default "BFGS".
#' @param lo lower bound(s) for the parameter domain
#' @param up upper bound(s) for the parameter domain
#' @param \dots other arguments to be passed to \code{rcbr.fit} to control fitting.
#' @return a list comprising the components:
#' \describe{
#'	\item{bopt}{output of the optimizer for the profiled parameters beta}
#'	\item{fopt}{output of the optimizer for the random coefficients eta}
#' }
#' @importFrom stats optim
#' @importFrom stats stepfun
#' @importFrom stats logLik
#' @importFrom Formula Formula
#' @export
prcbr <- function (formula, b0, data, logL = TRUE, omethod = "BFGS", lo = -Inf, 
    up = Inf, ...) 
{
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0)
    mf <- mf[c(1, m)]
    F <- Formula::Formula(formula)
    mf[[1]] <- as.name("model.frame")
    mf$formula <- F
    mf <- eval(mf, parent.frame())
    y <- model.response(mf)
    X <- model.matrix(F, data = mf, rhs = 1)
    Z <- model.matrix(F, data = mf, rhs = 2)[, -1, drop = FALSE]
    p <- NCOL(Z)
    b0 <- matrix(b0, p)
    M <- ncol(b0)
    cntl <- KW.control(...)
    eps <- .Machine$double.eps^(2/3)
    plogL <- function(b, X, y, Z, cntl) {
        -logLik(rcbr.fit(X, y, offset = Z %*% b, mode = "KW", 
            cntl))
    }
    pglogL <- function(b, X, y, Z, cntl) {
        f <- rcbr.fit(X, y, offset = Z %*% b, mode = "KW", cntl)
        if (class(f)[[1]] == "KW1") {
            Fhat <- stepfun(f$x, cumsum(c(0, f$y)))
            FZb <- Fhat(X[, 2] + Z %*% b)
        }
        else if (class(f)[[1]] == "KW2") {
            Fhatu <- stepfun(f$uv[, 1], cumsum(c(0, f$W)))
            Fhatv <- stepfun(f$uv[, 2], cumsum(c(0, f$W)))
            Fhat <- function(u, v) pmin(Fhatu(u), Fhatv(v))
            FZb <- Fhat(X[, 3], X[, 2] + Z %*% b)
        }
        FZb <- FZb * ((eps < FZb) && (FZb < (1 - eps)))
        c(crossprod(crossprod(Z, y - FZb)/length(y)))
    }
    if (logL) {
        if (M == 1) {
            g <- optim(b0, plogL, method = omethod, X = X, y = y, Z = Z, 
		       lower = lo, upper = up, cntl = cntl)
            bhat <- g$par
        }
        else {
            g <- apply(b0, 2, FUN = function(b) plogL(b, X, y, 
                Z, cntl = cntl))
            bhat <- b0[, which(g == max(g))]
        }
        f <- rcbr.fit(X, y, offset = Z %*% bhat, mode = "KW", 
            cntl)
    }
    else {
        if (M == 1) {
            g <- optim(b0, pglogL, method = omethod, X = X, y = y, Z = Z, 
		       lower = lo, upper = up, cntl = cntl)
            bhat <- g$par
        }
        else {
            g <- apply(b0, 2, FUN = function(b) pglogL(b, X, 
                y, Z, cntl = cntl))
            bhat <- b0[, which(g == min(g))]
        }
        f <- rcbr.fit(X, y, offset = Z %*% bhat, mode = "KW", 
            cntl)
    }
    list(bopt = g, fopt = f)
}
