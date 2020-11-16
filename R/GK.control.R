#' Control parameters for Gautier-Kitamura bivariate random coefficient binary response
#'
#' These parameters can be passed via the \code{...} argument of the \code{rcbr} function.
#' defaults as suggested in Gautier and Kitamura matlab code
#'
#' @param n  the sample size 
#' @param u  grid values for intercept coordinate
#' @param v  grid values for slope coordinate
#' @param T  Truncation parameter for numerator must grow "sufficiently slowly with n"
#' @param TX  Truncation parameter for denomerator must grow "sufficiently slowly with n"
#' @param Mn  Trimming parameter "chosen to go to 0 slowly with n"
#'
#' @return updated list
#' @export
GK.control <- function(n, u = -20:20/10, v = -20:20/10, T = 3, TX = 10, Mn = 1/log(n)^2) { 
    list(u = u, v = v, T = T, TX = TX, Mn = Mn, tol = .Machine$double.eps^0.5)
}
