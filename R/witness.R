#' Find witness point
#'
#' Find (if possible) an interior point of a polytope solving a linear program
#'
#' Solves LP:  max over {w,eps} {eps | SAw - eps >= Sb, 0 < eps <= epsbound}
#' S is diag(s),  if at the solution eps > 0, then w is a valid interior point
#' otherwise the LP fails to find an interior point, another s must be tried.
#' Constructs a problem formulation that can be passed to Rmosek for solution.  
#'
#' @param  A Is a n by d matrix of hyperplane slope coefficients.
#' @param  b Is an n vector of hyperplane intercept coefficients.
#' @param  s Is an n vector of signs. 
#' @param  epsbound Is a scalar tolerance controlling how close the witness 
#'	point can be to an edge of the polytope.
#' @param  epstol Is a scalar tolerance for the LP convergence.
#' @param  presolve Controls whether Mosek should presolve the LP.
#' @param  verb Controls verbosity of Mosek solution.
#' @return List with components:
#'  \itemize{
#'	\item  w proposed interior point at solution 
#'	\item  fail indicator of whether w is a valid interior point
#'  }
#' @export
witness <- function(A,b, s,epsbound = 1,epstol = 1e-07, presolve = 1, verb = 0){
  if(!is.list(s)) s <- list(s)
  m = ncol(A)
  n = nrow(A)
  if (is.null(m)){  
  	m = length(A)
  	n = 1
  }
  # Interior defines a Linear Programming Problem 	
  Interior <- list()
  Interior$sense <- "max"
  Interior$c <- c(rep(0, m), 1)
  blx <- c(rep(-Inf,m),0)    # lower bound of eps set as 0 
  bux <- c(rep(Inf,m),epsbound)
  Interior$bx <- rbind(blx,bux)
  lapply(s, FUN=function(item){
    S <- if (length(item) > 1) diag(item) else item
    A = cbind(S %*% A, matrix(-1, nrow = n,1))
    Interior$A <- as(A, 'dgCMatrix')
    buc <- c(rep(Inf, n))
    blc <- c(S %*% b)
    Interior$bc <- rbind(blc,buc)
    Interior$iparam$MSK_IPAR_PRESOLVE_USE <- presolve  
    r <- Rmosek::mosek(Interior,opts = list(verbose = verb))
    w = r$sol$bas$xx[1:m]
    eps <- r$sol$bas$xx[m+1]  
    fail <- eps <= epstol
    list(w = w, fail = fail, s = item, eps = eps)  
  })
}
