#' Prediction of Bounds on Marginal Effects
#'
#' Given a fitted model by the exact NPMLE procedure prediction is made
#' at a new design point with lower and upper bounds for the prediction due
#' to ambiguity of the assignment of mass within the cell enumerated polygons.
#'
#' @param  object is the fitted NPMLE object
#' @param  \dots  is expected to contain an argument \code{newdata}
#' @return a list consisting of the following components:
#'	\describe{
#'		\item{phat}{Point prediction}
#'		\item{lower}{lower bound prediction}
#'		\item{upper}{upper bound prediction}
#'		\item{xpoly}{indices of crossed polygons}
#'	}
#' @seealso \code{predict.KW2} for a simpler prediction function without bounds
#' @author Jiaying Gu
#' @export
bounds.KW2 <- function(object, ...){
    x <- object$x
    p <- ncol(x)
    dots <- list(...)
    newdata <- dots$newdata
    if(!length(newdata))
	stop("No newdata to predict at.")
    if (length(newdata) > p) 
	stop("Only one new data point at a time!")
    A <- x[,1:(p-1)]
    b <- -x[,p]
    SignVector = object$cells$SignVector
    w <- object$cells$w
    locmax <- object$locmax
    duprec <- which(A[,2]/A[,1] == newdata[2]/newdata[1] & b/A[,1]== -newdata[3]/newdata[1])
    if (length(duprec) > 0){
	A = A[-duprec,,drop=FALSE]
	b = b[-duprec]
	SignVector = SignVector[-duprec,]
	}
    Anew = rbind(A,newdata[1:(p-1)])
    bnew = c(b, -newdata[p])
    nobs <- nrow(Anew)
    tempdim = sum(!is.na(SignVector[1,]))
    SignVector <- rbind(SignVector,NA)
    SignVector[nrow(SignVector),1:tempdim] <- 
	sign(as.vector(t(w[,1:tempdim]) %*% Anew[nobs,] - bnew[nobs]))
    if(any(SignVector[nobs,1:tempdim]==0)) 
	stop("Error: newdata falls on existing interior points: try perturbing newdata")
    colset <- polyzone(SignVector, w, Anew, bnew)
    locmaxset <- c(1:ncol(SignVector))[locmax]
    crossedpoly  <- sort(colset[which(colset%in%locmaxset)])
    crossedpolyf <- object$W[which(locmaxset%in%colset)]
    currentsign <- 0.5*(SignVector[nobs, locmaxset]+1)
    currentpred <- currentsign %*% object$W
    crossedpolysign <- SignVector[nobs, crossedpoly]
    upperpred <- sum(crossedpolyf[which(crossedpolysign== -1)]) + currentpred
    lowerpred <- currentpred - sum(crossedpolyf[which(crossedpolysign == 1)])
    list(phat = currentpred, upper = upperpred, lower = lowerpred, xpoly = crossedpoly)
}

