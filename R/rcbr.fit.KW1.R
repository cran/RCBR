#' NPMLE fitting for the Cosslett random coefficient binary response model
#'
#' This is the original one dimensional version of the Cosslett model, also
#' known as the current status model:
#' \deqn{P(y = 1 | v) = \int I (\eta > v)dF(\eta).}
#' invoked with the formula \code{y ~ v}.  By default the algorithm computes a vector
#' of potential locations for the mass points of \eqn{\hat F} by finding interior
#' points of the intervals between the ordered \code{v}, and then solving a convex
#' optimization problem to determine these masses.  Alternatively, a vector of
#' predetermined locations can be passed via the control argument.  Additional
#' covariate effects can be accommodated by either specifying a fixed offset in
#' the call to \code{rcbr} or by using the profile likelihood function \code{prcbr}.
#' 
#' @param  X the design matrix expected to have an intercept column of
#' ones as the first column, the last column is presumed to contain values of
#' the covariate that is designated to have coefficient one.
#' @param  y the binary response.
#' @param  control is a list of parameters for the fitting, see
#' \code{KW.control} for further details.
#' @return a list with components:
#' \itemize{
#'   \item  x evaluation points for the fitted distribution
#'   \item  y estimated mass associated with the \code{v} points
#'   \item  logLik the loglikelihood value of the fit
#'   \item  status mosek solution status
#' }
#' @author  Jiaying Gu and Roger Koenker 
#' @references
#' Gu, J. and R. Koenker (2018)  Nonparametric maximum likelihood estimation 
#' of the random coefficients binary choice model, preprint.
#' @keywords  nonparametrics  
#' @export
rcbr.fit.KW1 <- function (X, y, control) {
    x <- X[,2]
    v <- control$v
    n <- length(y)
    if (length(v) == 1) 
	v <- seq(min(x) - 1e-4, max(x) + 1e-4, length = v)
    shift <- 0
    if(!length(v)){
    	if (length(unique(x))==n) u <- sort(x)
	else{ 
	    # when there are ties in x and their y disagree, break the tie by 
	    # shifting those x with y =0 to the right by a magnitude of shift
	    shift <- min(abs(diff(unique(sort(x)))))/2  
	    nn <- names(table(x))[which(tapply(y,x,mean) != 0 & tapply(y,x,mean)!=1)] 
	    x[x %in% nn & y==0] <- x[x %in% nn & y==0]+shift
	    u <- sort(unique(x))
	    }
	m <- length(u)
	c <- rep(0,m+1)
	for(k in 1:m)
	    c[k] = sum(x[y==0]< u[k]) + sum(x[y==1]>= u[k])
	c[m+1] = sum(x[y==0]<= u[m])
	if (sum(diff(c)==0)>0){ 
	    # counts has flat regions, reduce to unique counts, and recover locations of local maxima
	    if (sum(diff(c)==0)==(length(c)-1))
        	localmax = 1:length(c)
	    else{
		loc <- which(diff(c)!=0) 
		id <- rep(1:(length(loc)+1), c(loc[1],diff(loc), length(c)-loc[length(loc)]))
		creduce = c[c(loc,max(loc)+1)]
		s <- sign(diff(creduce))
		ns = length(s)
		lmax <- which(s[-ns] - s[-1] == 2)
		if(s[1] == -1) lmax = c(0,lmax)
		if(s[ns] == 1 ) lmax = c(lmax,ns)
		localmax = c(1:length(c))[id%in%(lmax+1)]
	    }
	}
	else{
	    s <- sign(diff(c))
	    ns = length(s)
	    lmax <- which(s[-ns] - s[-1] == 2)
	    if(s[1] == -1) lmax = c(0,lmax)
	    if(s[ns] == 1 ) lmax = c(lmax,ns)
	    localmax = lmax+1
	    }
	# use right end point of the partition as support point. 
	# For the last open interval use max(u) + shift
	v <- c(u, max(u) + shift)[localmax]  
    }
    w <- rep(1, n)/n
    d <- rep(1, length(v))
    A <- outer(x, v, ">=")
    A <- (y == 1) * A + (y == 0) * (1 - A)
    if (sum(rowSums(A)==0)>0){
    	warning("zero probability occurs, switched to full data support")
    	v = sort(unique(x))
    	}
    d <- rep(1, length(v))
    A <- outer(x, v, ">=")
    A <- (y == 1) * A + (y == 0) * (1 - A)
    f <- KWDual(A, d, w)
    logLik <- sum(log(A %*% f$f))
    z <- list(x = v, y = f$f, logLik = logLik, status = f$status)
    class(z) <- c("KW1", "density")
    return(z)
}
