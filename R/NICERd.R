#' New (Accelerated) Incremental Cell Enumeration (in) R
#'
#' Find interior points and cell counts of the polygons (polytopes) formed by a
#' hyperplane arrangement.
#'
#' Modified version of the algorithm of Rada and Cerny (2018).
#' The main modifications include preprocessing as hyperplanes are added
#' to determine which new cells are created, thereby reducing the number of
#' calls to the witness function to solve LPs, and treatment of degenerate
#' configurations as well as those in "general position." (for \eqn{d=2} for now).  When the hyperplanes
#' are in general position the number of cells (polytopes) is determined by the
#' elegant formula of Zaslavsky (1975) \deqn{m = {n \choose d} + n + 1}.  In
#' degenerate cases, i.e. when hyperplanes are not in general position, the
#' number of cells is more complicated as considered by Alexanderson and Wetzel (1981).
#' The function \code{polycount} is provided to check agreement with their results
#' in an effort to aid in the selection of tolerances for the \code{witness} function for arrangement in \eqn{d=2}.
#' The current version is intended mainly for use with \eqn{d = 2}, but the algorithm is adapted to
#' the general position setting with \eqn{d > 2}, although it requires hyperplanes in general position and may require some patience when both
#' the sample size is large.
#' if hyperplanes not general position (i.e. all cross at origin), turn off accelerate
#'
#' @param  A is a n by d matrix of hyperplane slope coefficients
#' @param  b is an n vector of hyperplane intercept coefficients
#' @param  initial origin for the interior point vectors \code{w}
#' @param  epsbound is a scalar tolerance controlling how close the witness point
#'	can be to an edge of the polytope
#' @param  epstol is a scalar tolerance for the LP convergence
#' @param  verb controls verbosity of Mosek solution 
#' @param accelerate allows the option to turn off acceleration step (turned off by default)
#' @return  A list with components:
#' \itemize{
#'	\item  SignVector a n by m matrix of signs determining position of cell relative 
#'		to each hyperplane.
#'	\item  w a d by m matrix of interior points for the m cells
#'	}
#' @references
#' Alexanderson, G.L and J.E. Wetzel, (1981) Arrangements of planes in space, 
#'	Discrete Math, 34, 219--240.
#' Rada, M. and M. Cerny (2018) A new algorithm for the enumeration of cells 
#'	of hyperplane arrangements and a comparison with Avis and Fukada's reverse 
#'	search, SIAM J. of Discrete Math, 32, 455-473.
#' Zaslavsky, T. (1975) Facing up to arrangements:  Face-Count Formulas for 
#'	Partitions of Space by Hyperplanes, Memoirs of the AMS, Number 154.
#' @export
#'
NICERd <- function(A, b, initial = rep(0, ncol(A)), verb = TRUE, accelerate = FALSE, epsbound = 1, epstol = 1e-07){
  if (sum(abs(A %*% initial - b)) < epstol) 
      stop("Error: initial point falls on some of the hyperplanes")
  m <- ncol(A)  # if arrangement in 3d then m = 3 
  n <- nrow(A)
  maxdim <- sum(choose(n,0:m))   # maximum column dimension of SignVector 
  SignVector <- matrix(NA, n, maxdim)  
  w <- matrix(0, ncol(A), maxdim)  
  eps <- matrix(0,1,maxdim)
  mode(SignVector) <- "integer" 
  i <- 0 
  if (t(as.vector(A[(1:(i+1)),]))%*%as.vector(initial)>b[(1:(i+1))]){
    SignVector[(i+1),1] <- 1L
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
    # pre-process to see which polygons creates new cells with the addition the hyperplane (normalize each plane by the first non-zero entry of the row)
    # otherwise, no new polygons should be added and move on. 
    dup <- duplicated(cbind(A[1:(i+1),]/A[1:(i+1),which.max(abs(A[i+1,])>0)], b[1:(i+1)]/A[1:(i+1),which.max(abs(A[i+1,])>0)]), MARGIN=1)
    if (dup[i+1] == 0){
      if (i > m-1){
      	if(accelerate){
        combs = combn((i+1), m)
		combs = combs[,which(combs==(i+1),arr.ind=TRUE)[,2]]
		vertex <- matrix(NA, nrow = m, ncol = ncol(combs))
  
        for (k in 1:ncol(combs)){ 
            Atemp = A[combs[,k],]
			btemp = b[combs[,k]]
			vertexsolve = try(solve(Atemp)%*%btemp, silent = TRUE)
			if(inherits(vertexsolve, "try-error")) stop("Hyperplanes not in general position, turn off accelerate")
			vertex[,k] = vertexsolve
        }
        if (sum(!is.na(vertex))==0){
          colset = c(1:tempdim)
        }else{ # feas records sign constraints of each existing hyperplane wrt to each vertex
          feas <- matrix(NA_integer_, nrow = i, ncol = ncol(combs))  
          Atemp = A[1:i,]  
          btemp = b[1:i] 
          for (j in 1:ncol(combs)){
            if (!is.na(vertex[1,j])){
              feas[(1:i)[!(1:i)%in%combs[,j]],j] = sign(round(Atemp[!((1:i)%in%combs[,j]),,drop=FALSE]%*%vertex[,j]-btemp[!((1:i)%in%combs[,j])],digits = 12)) 
	    }else{
              feas[,j] <- NA_integer_
          }
          }
          rownames(feas) <- 1:i
          feas <- feas[which(rowSums(is.na(feas))!=i),,drop=FALSE]
          feas[which(feas==0)] <- NA_integer_
          feas <- unique(feas,MARGIN=1)
          tempSignVector <- SignVector[1:i, 1:tempdim]
          # if all rows have only (m-1) NA (i.e. general position case)
          # enumerate the NA by 1 and -1 exhausts all possible polygons that may have been crossed 
          # otherwise throw up error  
          # colset collects all polygons that are crossed by the newly added hyperplane
          if (sum(apply(feas,2,function(x) sum(is.na(x))==(m-1)))==ncol(feas)){
          	feaset = NULL
            li <- expand.grid(rep(list(c(-1,1)), m-1))
            for (j in 1:ncol(feas)){
            	feastemp = t(matrix(rep(feas[,j], nrow(li)), nrow = i))
            	feaset = c(feaset, lapply(1:nrow(li), FUN = function(x) replace(feastemp[x,], is.na(feastemp[x,]), as.numeric(li[x,]))))
            	}
            feaset = matrix(unlist(feaset), nrow = i)
            feaset = unique(feaset,MARGIN = 2)
            colset <-which(duplicated(cbind(tempSignVector,feaset),fromLast=TRUE,MARGIN=2))
          }
	  else{
	  	stop("Hyperplanes not in general position, turn off acceleration")
	  	}
        }
     }else{
        colset = c(1:tempdim)
        }
        }
        else{
        colset = c(1:tempdim)
        }
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
  if (dim(SignVector)[2] != sum(choose(n, 0:m))) 
      print("warning: number of cells does not reach upper bound, either hyperplanes not in general position or consider changing epsbound.")
  list(w = w, SignVector = SignVector, eps = eps, A = A, b = b)
}
