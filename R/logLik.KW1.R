#' log likelihood for KW1 procedure
#' @param object a fitted object of class "KW1"
#' @param ... other parameters for logLik
#' @return a scalar log likelihood
#' @export
#'
logLik.KW1 <- function(object, ...) object$logLik
