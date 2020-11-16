#' New Incremental Cell Enumeration (in) R
#'
#' Find interior points and cell counts of the polygons (cells) formed by a
#' line arrangement.  
#'
#' Modified version of the algorithm of Rada and Cerny (2018).
#' The main modifications include preprocessing as hyperplanes are added
#' to determine which new cells are created, thereby reducing the number of
#' calls to the witness function to solve LPs, and treatment of degenerate
#' configurations as well as those in "general position."  When the hyperplanes
#' are in general position the number of polytopes (cells) is determined by the
#' elegant formula of Zazlavsky (1975) \deqn{m = {n \choose d} + n + 1}.  In
#' degenerate cases, i.e. when hyperplanes are not in general position, the
#' number of cells is more complicated as considered by Alexanderson and Wetzel (1981).
#' The function \code{polycount} is provided to check agreement with their results
#' in an effort to aid in the selection of tolerances for the \code{witness} function.
#' Current version is intended for use with \eqn{d = 2}, but the algorithm is adaptable to
#' \eqn{d > 2}, and there is an experimental version called \code{NICERd} in the package.
#'
#' @param  A is a n by 2 matrix of slope coefficients
#' @param  b is an n vector of intercept coefficients
#' @param  initial origin for the interior point vectors \code{w}
#' @param  epsbound is a scalar tolerance controlling how close the witness point
#'	can be to an edge of the polytope
#' @param  epstol is a scalar tolerance for the LP convergence
#' @param  verb controls verbosity of Mosek solution 
#'
#' @return  A list with components:
#' \itemize{
#'	\item  SignVector a n by m matrix of signs determining position of cell relative 
#'		to each hyperplane.
#'	\item  w a d by m matrix of interior points for the m cells
#'	}
#' @examples{
#' if(packageVersion("Rmosek") > "8.0.0"){
#'     A = cbind(c(1,-1,1,-2,2,1,3), c(1,1,1,1,1,-1,-2))
#'     B = matrix(c(3,1,7,-2,7,-1,1), ncol = 1)
#'     plot(NULL,xlim = c(-10,10),ylim = c(-10,10))
#'     for (i in 1:nrow(A))
#' 	  abline(a = B[i,1]/A[i,2], b = -A[i,1]/A[i,2],col = i)
#'     f = NICER(A, B)
#'     for (j in 1:ncol(f$SignVector))
#'     	  points(f$w[1,j], f$w[2,j], cex = 0.5)
#'     }
#' }
#' @references
#' Alexanderson, G.L and J.E. Wetzel, (1981) Arrangements of planes in space, 
#'	Discrete Math, 34, 219--240.
#' Gu, J. and R. Koenker (2020)  Nonparametric Maximum Likelihood Methods for
#' 	Binary Response Models with Random Coefficients, \emph{J. Am. Stat Assoc}
#' Rada, M. and M. Cerny (2018) A new algorithm for the enumeration of cells 
#'	of hyperplane arrangements and a comparison with Avis and Fukada's reverse 
#'	search, SIAM J. of Discrete Math, 32, 455-473.
#' Zaslavsky, T. (1975) Facing up to arrangements:  Face-Count Formulas for 
#'	Partitions of Space by Hyperplanes, Memoirs of the AMS, Number 154.
#' @export
#'
NICER <- function(A, b, initial = c(0,0), verb = TRUE, epsbound = 1, epstol = 1e-07){
  if (ncol(A) > 2) stop("Error: NICER only works for dim = 2, use NICERd")
  if (sum(abs(A %*% initial - b)) < epstol) 
      stop("Error: initial point falls on some of the hyperplanes")
  n <- nrow(A)
  maxdim <- choose(n,2) + n + 1   # maximum column dimension of SignVector for arrangement in R2
  SignVector <- matrix(NA, n, maxdim)  
  w <- matrix(0, 2, maxdim)  
  eps <- matrix(0,1,maxdim)
  mode(SignVector) <- "integer" 
  i <- 0 
  if (t(as.vector(A[(1:(i+1)),]))%*%as.vector(initial)>b[(1:(i+1))])
    {SignVector[(i+1),1] <- 1L
  }else{
    SignVector[(i+1),1] <- -1L}
  w[,1] <- initial
  eps[,1] <- epsbound
  s <- SignVector[1:(i+1),1] * c(rep(1, i),-1)
  test <- witness(A[(1:(i+1)),], b[(1:(i+1))], s = s, epsbound = epsbound, epstol = epstol)
  if (!test[[1]]$fail){
    SignVector[(i+1),2] <- test[[1]]$s
    w[,2] <- c(test[[1]]$w)
    eps[,2] <- c(test[[1]]$eps)
  }
  # for each i columns of SignVector and w are filled with values
  # test if the corresponding w is on the upper side of the hyperplane 
  for (i in 1:(n-1)){ # adding hyperplane 2 to n
    tempdim = sum(!is.na(SignVector[1,]) )  
    SignVector[i+1,1:tempdim] <- sign(as.vector(t(w[,1:tempdim])%*%A[i+1,] - b[i+1]))  
    # the following lines avoid the situation that the newly added line crosses 
    # any of the existing interior points. If so, find a new interior point 
    # N.B. seems using the same epsbound does not help, so I've tentatively 
    # set it at 0.1, this needs further testing)
    if(any(SignVector[i+1,1:tempdim]==0)){
      for(k in which(SignVector[i+1,1:tempdim]==0)){
        testk <- witness(A[1:(i+1),], b[1:(i+1)], s= c(SignVector[1:i, k],1), epsbound=epsbound, epstol=epstol)
        if(!testk[[1]]$fail){
          w[,k] <- c(testk[[1]]$w)
          eps[,k] <- c(testk[[1]]$eps)
	  signs <- sign(t(w[,k])%*%A[i+1,] - b[i+1])
	  mode(signs) <- "integer"
          SignVector[i+1,k] <- signs
        }
      }
    }
    # if the newly added line does not coincide with previous existing lines
    # pre-process to see which polygons creates new cells with the addition the hyperplane
    # otherwise, no new polygons should be added and move on. 
    dup <- duplicated(cbind(A[1:(i+1),]/A[1:(i+1),which.max(abs(A[i+1,])>0)], b[1:(i+1)]/A[1:(i+1),which.max(abs(A[i+1,])>0)]), MARGIN=1)
    if (dup[i+1] == 0){
    	colset <- polyzone(SignVector, w, A[1:(i+1),], b[1:(i+1)])
    # Find new interior point by flipping the sign of the last entry in relevant SignVector    
      s = lapply(colset, FUN=function(x) SignVector[1:(i+1),x] * c(rep(1,i),-1))
      testlist <- witness(A[(1:(i+1)),], b[(1:(i+1))], s = s, epsbound = epsbound, epstol = epstol)
      testlist.keep = which(!unlist(lapply(testlist,"[[",2))) # keep only those that has fail = FALSE
      if (length(testlist.keep)>0){
        for(l in 1:length(testlist.keep)){
          SignVector[1:(i+1),tempdim+l] <- testlist[[testlist.keep[l]]]$s
          w[,tempdim+l] <- c(testlist[[testlist.keep[l]]]$w)
          eps[,tempdim+l] <- c(testlist[[testlist.keep[l]]]$eps)
        }
      }
    }
    if ((verb != 0) && round((i+1)/10)==(i+1)/10) 
	print(paste0(i+1, " :: Polygons: ", sum(!is.na(SignVector[1,]) )), quote=FALSE)				
  }
  w = w[,colSums(is.na(SignVector))<nrow(SignVector)]
  SignVector = SignVector[,colSums(is.na(SignVector))<nrow(SignVector)]
  eps = eps[,1:ncol(w)]
  if (dim(SignVector)[2] != polycount(A,b)) 
      print("warning: polycount does not match, consider changing epsbound.")
  list(w = w, SignVector = SignVector, eps = eps, A = A, b = b)
}
