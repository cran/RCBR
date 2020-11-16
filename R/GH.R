#' Current Status Linear Regression
#'
#' Groeneboom and Hendrickx semiparametric binary response estimator (scalar case)
#' score estimator based on NPMLE avoids any smoothing
#' proposed by Groneboom and Hendrickx (2018).  
#' 
#' @param b  parameter vector (fix last entry as a known number, usually 1 or -1, for normalization) 
#' @param X  design matrix
#' @param y  binary response vector
#' @param eps trimming tolerance parameter
#' @return  A list with components:
#'	\itemize{
#'	\item  evaluation of a score function at parameter value
#'	\item  estimated standard error 
#'	\item  sindex single index linear predictor
#'	}
#' @references
#'	Groeneboom, P. and K. Hendrickx (2018) Current Status Linear Regression, 
#'  Annals of Statistics, 46, 1415-1444,
#' @importFrom REBayes Cosslett
#' @export
GH <- function(b, X, y, eps = 0.001){
	# input b is a vector with the last entry fixed (normalization)
	sindex = as.numeric(X%*%b)
	o = order(sindex)
	fstar <- Cosslett(sindex, y, v = sort(sindex))
	Fhat = cumsum(fstar$y)
	crossprod(X[o,1], (y[o] - Fhat)*(1*(Fhat >= eps & Fhat <= 1-eps)))/nrow(X)
}
#' Current Status Linear Regression Standard Errors
#'
#' Groeneboom and Hendrickx semiparametric binary response estimator (scalar case)
#' score estimator based on NPMLE avoids any smoothing
#' proposed by Groneboom and Hendrickx (2018).  
#' 
#' @param bstar  parameter vector (fix last entry as a known number, usually 1 or -1, for normalization) 
#' @param X  design matrix
#' @param y  binary response vector
#' @param hc  kernel bandwidth (used for the standard error estimation)
#' @param eps trimming tolerance parameter
#' @return  A list with components:
#'	\itemize{
#'	\item  evaluation of a score function at parameter value
#'	\item  estimated standard error 
#'	\item  sindex single index linear predictor
#'	}
#' @references
#'	Groeneboom, P. and K. Hendrickx (2018) Current Status Linear Regression, 
#'  Annals of Statistics, 46, 1415-1444,
#' @importFrom stats sd 
#' @importFrom stats var 
#' @importFrom REBayes Cosslett
#' @export
GH.se <- function(bstar, X, y, eps = 0.001, hc = 2){
	# covariance involves a density function, approximated by kernel smoothed estimates. 
	# input b is a vector with the last entry fixed (normalization)
	p = ncol(X) 
	if (p > 2) stop("only scalar case allowed")
	dX = X[,1]-mean(X[,1])
	sindex = as.numeric(X%*%bstar)
	o = order(sindex)
	n = length(y)
	h = hc * sd(sindex) * n^(-1/7)
	K <- function(u){
	(1*(abs(u)<=1))*35*(1-u^2)^3/32
	}
	fstar <- Cosslett(sindex, y, v = sort(sindex))
	Fhat <- cumsum(fstar$y)
	f <- rep(NA, n)
	for (i in 1:n){
		Kz <- K((sindex[i] - sort(sindex))/h)/h
		f[i] <- crossprod(Kz, fstar$y)
		}
	fhat = f[o]	
	A = mean(fhat*(dX[o]^2)*(1*(Fhat >= eps & Fhat <= 1-eps)))
	#B = sum(Fhat[-trim] * (1-Fhat[-trim]) * (X[o,2][-trim]-mean(X[-trim,2]))^2)/nrow(X)
	B = var(X[o,1] * (y[o] - Fhat)*(1*(Fhat >= eps & Fhat <= 1-eps)))
	sqrt(B/(A^2)/nrow(X))	
	}		
	
