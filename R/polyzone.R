#' Identify crossed polygons from existing cells when adding a new line (works only for dim = 2)
#' 
#' Given an existing cell configuration represented by the Signvector and
#' associated interior points w, identify the polygons crossed by the next
#' new line.
#'
#' @param SignVector current SignvVctor matrix
#' @param w associated interior points
#' @param A design matrix for full problem aka [1,z]
#' @param b associated final column of design matrix aka [v]
#' @return vector of indices of crossed polygons
#' @author Jiaying Gu
#' @export
polyzone <- function(SignVector, w, A, b){
	if (ncol(A) > 2) stop("polyzone only works for dim = 2")
    tempdim = sum(!is.na(SignVector[1,]))
    N = nrow(A)
    if (N >2){
	vertex <- matrix(NA, nrow = 2, ncol = (N-1))  
	for (k in 1:(N-1)){
    	   if (A[k,1]/A[k,2] != A[N,1]/A[N,2])  # new line not parallel to existing lines 
               vertex[,k] <- solve(rbind(A[k,], A[N,]), rbind(b[k], b[N]))  
        }
    if (sum(!is.na(vertex))==0)
      colset = c(1:tempdim)
    else{ # feas records sign constraints of each existing hyperplanes wrt to each vertex
      feas <- matrix(NA_integer_, nrow = (N-1), ncol = (N-1))  
      Atemp = A[1:(N-1),]  
      btemp = b[1:(N-1)] 
      for (j in 1:(N-1)){
        if (!is.na(vertex[1,j]))
          feas[j, -j] <- sign(Atemp[-j,]%*%vertex[,j]-btemp[-j])  
	else
          feas[j,] <- NA_integer_
      }
      rownames(feas) <- 1:(N-1)
      feas <- feas[which(rowSums(is.na(feas))!=(N-1)),,drop=FALSE]
      feas[which(feas==0)] <- NA_integer_
      feas <- unique(feas,MARGIN=1)
      tempSignVector <- SignVector[1:(N-1), 1:tempdim]
      # if all rows have only one NA (i.e. are in general position)
      # enumerate the NA by 1 and -1 exhausts all possible polygons that may have been crossed 
      # otherwise, distinguish rows that have one NA (set1na) and more than one NA (set2na)
      # set2na only occurs when the new added hyperplane crosses an existing vertex, 
      # colset collects all polygons that are crossed by the new added hyperplane
      if (sum(is.na(feas))==nrow(feas)){
        fea1 = feas
        fea1[is.na(fea1)] = 1L
        fea2 = feas
        fea2[is.na(fea2)] = -1L
        feaset = rbind(fea1,fea2)
        feaset = unique(t(feaset),MARGIN = 2)
        colset <-which(duplicated(cbind(tempSignVector,feaset),fromLast=TRUE,MARGIN=2))
      }
      else{
        set1na <- rowSums(is.na(feas))==1
        set2na <- rowSums(is.na(feas))>1
        fea1na1 <- feas[set1na,]
        fea1na1[is.na(fea1na1)] = 1L
        fea1na2 <- feas[set1na,]
        fea1na2[is.na(fea1na2)] = -1L
        fea1naset = rbind(fea1na1,fea1na2)
        fea1naset = unique(t(fea1naset),MARGIN=2)
        colset1na <- which(duplicated(cbind(tempSignVector,fea1naset), fromLast = TRUE, MARGIN = 2))
        fea2na <- feas[set2na,,drop=FALSE]
        constraintmat <- abs(fea2na)==1  
        constraintmat <- constraintmat & !is.na(constraintmat) # picking out location of constraints in feas
        step <- 0.1*c(1,-1)/A[N,] 
        colsettemp <- lapply(1:nrow(constraintmat), FUN=function(j){
          constraint <- constraintmat[j,]
          realj <- as.integer(rownames(constraintmat)[j])
          clines <- which(!constraint)
          clines <- clines[clines!=realj]
          linecols <- which(colSums(abs(tempSignVector[constraint,,drop=FALSE]-fea2na[j,constraint]))==0)
          remcols <- logical(length(linecols))
          for(line in clines){
            vert <- solve(A[c(N,line),],b[c(N,line),drop=FALSE])
            vert_m <- vert - step
            vert_p <- vert + step
            goodsigns <- sign(A[c(realj,line),]%*%cbind(vert_m,vert_p) - b[c(realj,line)])
     	    csum1 <- colSums(abs(tempSignVector[c(realj,line),linecols,drop = FALSE]-goodsigns[,1])) != 0 
	    csum2 <- colSums(abs(tempSignVector[c(realj,line),linecols,drop = FALSE]-goodsigns[,2])) != 0
            remcols <- remcols | csum1 & csum2
          }
          linecols <- linecols[!remcols]
        })
        colset2na <- unique(unlist(colsettemp))
        colset <- unique(c(colset1na, colset2na))
        }
      }
   }
   else
     colset = c(1:2)
   colset
}
