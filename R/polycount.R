#' Check Cell Count  for degenerate hyperplane arrangements
#'
#' When the hyperplane arrangement is degenerate, i.e. not in general position,
#' the number of distinct cells can be checked against the formula of
#' Alexanderson and Wetzel (1981).
#'
#' @param  A is a n by m matrix of hyperplane slope coefficients
#' @param  b is an n vector of hyperplane intercept coefficients
#' @param  maxints  is maximum number of lines allowed to cross at the same vertex
#'
#' @return  number of distinct cells
#' @references
#' Alexanderson, G.L and J.E. Wetzel, (1981) Arrangements of planes in space, 
#'	Discrete Math, 34, 219--240.
#' @importFrom utils combn
#' @export
#'
polycount <- function(A, b, maxints=10){
  if (ncol(A) >2) stop("only implemented for dim = 2")
  Ab <- cbind(A,b)
  dup <- which(duplicated(Ab, MARGIN=1))
  if (length(dup)>0){
  	A = A[-dup,,drop=FALSE]
  	b = b[-dup]
  	}
  Ab = cbind(A,b)
  Ab2 = cbind(-A,-b)
  dup2 <- which(duplicated(rbind(Ab,Ab2),MARGIN=1))-nrow(A)
  if (length(dup2)>0){
    duploc <- (1:(length(dup2)/2))*2
    Ab = Ab[-dup2[duploc],,drop=FALSE]
  }
  A <- Ab[,1:2,drop=FALSE]
  b <- Ab[,3]
  if (nrow(A) < 2) stop("number of effective planes < 2")
  combs <- combn(nrow(A),2) # A here must have redundant rows removed to work (i.e. won't work if nrow(A) ==2)
  vertices <- vapply(1:ncol(combs), FUN=function(j){
    if(A[combs[1,j],1]/A[combs[1,j],2]!= A[combs[2,j],1]/A[combs[2,j],2])
      solve(A[combs[,j],],b[combs[,j]]) 
    else 
      rep(NA,2)}, FUN.VALUE=numeric(2))
  vertices.orig <- vertices
  vertices <- vertices[,!is.na(colSums(vertices))]
  vertices <- round(vertices,digits=12)
  if (is.null(dim(vertices))){
  	output = 4  # has only 2 lines
  	}else
  	{
  if (dim(vertices)[1]>0 &dim(vertices)[2]==0){
  	output = nrow(A)+1   # all lines parallel to each other
  	}else
  	{  
  	vertices <- vertices[,order(vertices[1,],vertices[2,])]
  whichdup <- which(duplicated(vertices,MARGIN=2))
  diffwhich <- diff(whichdup); diffwhich[diffwhich!=1] <- 0
  reps <- vapply(strsplit(paste0(diffwhich,collapse=""),"0")[[1]], FUN=nchar, FUN.VALUE=numeric(1))+2
  if(!all(as.numeric(names(table(reps))) %in% choose(3:maxints,2))) 
      stop("Error: inappropriate number of repetitions")
  counts <- c(NA, NA, sapply(3:maxints, FUN=function(i) sum(reps==choose(i,2))))
  counts[2] <- if(maxints==2) ncol(vertices) else ncol(vertices) - sum(choose(3:maxints,2)*counts[3:maxints])
  output = 1 + nrow(A) + sum((1:(maxints-1))*counts[2:maxints])
  }
  }
  output
}	
