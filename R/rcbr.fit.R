#' Fitting of Random Coefficient Binary Response Models
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
#' @param  x design matrix
#' @param  y binary response vector
#' @param  offset specifies a fixed shift in \code{v} representing the
#' potential effect of other covariates having fixed coefficients that may be
#' useful for profile likelihood computations.  (Should be vector of the same
#' length as \code{v}.
#' @param mode controls whether the Gautier and Kitamura, "GK", or Kiefer and
#' Wolfowitz, "KW" methods are used.
#' @param control  control parameters for fitting methods
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
#' Gautier, E. and Y. Kitamura (2013)  Nonparametric estimation in random coefficients
#' binary choice models, \emph{Ecoonmetrica}, 81, 581-607.
#'
#' Groeneboom, P. and K. Hendrickx (2016) Current Status Linear Regression,
#' preprint available from \url{https://arxiv.org/abs/1601.00202}.
#'
#' Ichimuma, H. and T. S. Thompson, (1998) Maximum likelihood estimation of a binary 
#' choice model with random coefficients of unknown distribution," 
#' \emph{Journal of Econometrics},  86, 269-295.
#'
#' @keywords nonparametric
#' @export
rcbr.fit <- function(x, y, offset = NULL, mode = "KW", control) {
    if(length(offset)) x[,ncol(x)] <- x[,ncol(x)] + offset
    if(mode == "KW") mode <- paste("KW", ncol(x) - 1, sep = "")
    switch(mode,
	   GK = rcbr.fit.GK(x, y, control),
	   KW1 = rcbr.fit.KW1(x, y, control), # Cosslett case
	   KW2 = rcbr.fit.KW2(x, y, control)
	   )
}
