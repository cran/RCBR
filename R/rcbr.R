#' Estimation of Random Coefficient Binary Response Models
#'
#' Two methods are implemented for estimating binary response models  
#' with random coefficients:  A nonparametric maximum likelihood method
#' proposed by Cosslett (1986) and extended by Ichimura and Thompson (1998),
#' and a (hemispherical) deconvolution method proposed by Gautier and
#' and Kitamura (2013).  The former is closely related to the NPMLE
#' for mixture models of Kiefer and Wolfowitz (1956).  The latter is an
#' R translation of the matlab implementation of Gautier and Kitamura.
#'
#' The \code{predict} method produces estimates of the probability of a "success"
#' (y = 1) for a particular vector, \code{(z,v)},  when aggregated over the estimated
#' distribution of random coefficients.
#'
#' The \code{logLik} produces an evaluation of the log likelihood value
#' associated with a fitted model.
#'
#' @param formula an expression of the generic form \code{y ~ z + v} where
#' \code{y} is the observed binary response, \code{z} is an observed covariate
#' with a random coefficient, and \code{v} is an observed covariate with
#' coefficient normalize to be one.  If \code{z} is not present then the
#' model has only a random "intercept" coefficient  and thus corresponds
#' to the basic model of Cosslett (1983);  this model is also referred to
#' as the current status model in the biostatistics literature, see Groeneboom
#' and Hendrikx (2016).  When \code{z} is present there are random coefficients
#' associated with both the intercept and \code{z}. 
#' @param data is a \code{data.frame} containing the data referenced in the
#' formula.
#' @param  subset specifies a subsample of the data used for fitting the model
#' @param  offset specifies a fixed shift in \code{v} representing the
#' potential effect of other covariates having fixed coefficients that may be
#' useful for profile likelihood computations.  (Should be vector of the same
#' length as \code{v}.
#' @param mode controls whether the Gautier and Kitamura, "GK", or Kiefer and
#' Wolfowitz, "KW" methods are used.
#' @param ...  miscellaneous other arguments to control fitting.  
#' See \code{GK.control} and \code{KW.control} for further details.
#'
#' @return of object of class \code{GK}, \code{KW1}, with components described in 
#' further detail in the respective fitting functions.
#'
#' @author  Jiaying Gu and Roger Koenker
#'
#' @references Kiefer, J. and J. Wolfowitz (1956) Consistency of the Maximum
#' Likelihood Estimator in the Presence of Infinitely Many Incidental
#' Parameters, \emph{Ann. Math. Statist}, 27, 887-906.
#'
#' Cosslett, S. (1983) Distribution Free Maximum Likelihood Estimator of the 
#' Binary Choice Model, \emph{Econometrica}, 51, 765-782.
#'
#' Gautier, E. and Y. Kitamura (2013)  Nonparametric estimation in random coefficients
#' binary choice models, \emph{Ecoonmetrica}, 81, 581-607.
#'
#' Gu, J. and R. Koenker (2020)  Nonparametric Maximum Likelihood Methods for 
#' Binary Response Models with Random Coefficients, \emph{J. Am. Stat Assoc} 
#'
#' Groeneboom, P. and K. Hendrickx (2016) Current Status Linear Regression,
#' preprint available from \url{https://arxiv.org/abs/1601.00202}.
#'
#' Ichimura, H. and T. S. Thompson, (1998) Maximum likelihood estimation of a binary 
#' choice model with random coefficients of unknown distribution," 
#' \emph{Journal of Econometrics},  86, 269-295.
#'
#' @examples{
#' if(packageVersion("Rmosek") > "8.0.0"){
#'     # Simple Test Problem for rcbr
#'     n <- 100
#'     B0 = rbind(c(0.7,-0.7,1),c(-0.7,0.7,1))
#'     z <- rnorm(n)
#'     v <- rnorm(n)
#'     s <- sample(0:1, n, replace = TRUE)
#'     XB0 <- cbind(1,z,v) %*% t(B0)
#'     u <- s * XB0[,1] + (1-s) * XB0[,2]
#'     y <- (u > 0) - 0
#'     D <- data.frame(z = z, v = v, y = y)
#'     f <- rcbr(y ~ z + v, mode = "KW", data = D)
#'     plot(f)
#'     # Simple Test Problem for rcbr
#'     set.seed(15)
#'     n <- 100
#'     B0 = rbind(c(0.7,-0.7,1),c(-0.7,0.7,1))
#'     z <- rnorm(n)
#'     v <- rnorm(n)
#'     s <- sample(0:1, n, replace = TRUE)
#'     XB0 <- cbind(1,z,v) %*% t(B0)
#'     u <- s * XB0[,1] + (1-s) * XB0[,2]
#'     y <- (u > 0) - 0
#'     D <- data.frame(z = z, v = v, y = y)
#'     f <- rcbr(y ~ z + v, mode = "GK", data = D)
#'     contour(f$u, f$v, matrix(f$w, length(f$u)))
#'     points(x = 0.7, y = -0.7, col = 2)
#'     points(x = -0.7, y = 0.7, col = 2)
#'     f <- rcbr(y ~ z + v, mode = "GK", data = D, T = 7)
#'     contour(f$u, f$v, matrix(f$w, length(f$u)))
#'     points(x = 0.7, y = -0.7, col = 2)
#'     points(x = -0.7, y = 0.7, col = 2)
#'     }
#' }
#' @keywords nonparametric
#' @importFrom stats model.matrix 
#' @importFrom stats model.response 
#' @importFrom stats terms
#' @export

rcbr <- function (formula, data, subset, offset, mode = "GK", ...)
{
    cl <- match.call()
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action", "weights",
        "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- terms(formula, data = data)
    Y <- model.response(mf, "any")
    X <- model.matrix(mt, mf)
    n <- length(Y)
    control <- switch(mode, GK = GK.control(n, ...), KW = KW.control(...))
    rval <- rcbr.fit(X,Y, mode = mode, control = control)
    return(rval)
}
