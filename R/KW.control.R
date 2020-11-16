#' Control parameters for NPMLE of bivariate random coefficient binary response
#'
#' These parameters can be passed via the \code{...} argument of the \code{rcbr} function.
#' The first three arguments are only relevant if full cell enumeration is employed for
#' bivariate version of the NPMLE.
#'
#' @param uv matrix of evaluation points for potential mass points
#' @param u grid of evaluation points for potential mass points
#' @param v grid of evaluation points for potential mass points
#' @param initial initial point for cell enumeration algorithm
#' @param epsbound controls how close witness points can be to vertices of a cell
#' @param epstol zero tolerance for witness solutions
#' @param presolve controls whether Mosek does a presolve of the LP
#' @param verb controls verbosity of Mosek solver 0 implies it is quiet
#'
#' @return updated list
#' @export
KW.control <- function(uv = NULL, u = NULL, v = NULL, initial = c(0,0), 
		       epsbound = 1, epstol = 1e-07, presolve = 1, verb = 0)
    list(uv = uv, u = u, v = v, initial = initial, epsbound = epsbound, 
	 epstol = epstol, presolve = presolve, verb = verb)
