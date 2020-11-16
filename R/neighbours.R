#' Check Neighbouring Cell Counts
#'
#' Compare cell counts for each cell with its neighbours and return indices
#' of the locally maximal cells.
#'
#' @param  SignVector  n by m matrix of signs produced by NICER
#'
#' @return Column indices of the cells that are locally maximal, i.e. those
#'	whose neighbours have strictly fewer cell counts.  The corresponding
#'	interior points of these cells can be used as potential mass points
#'	for the NPMLE function \code{rcbr.fit.KW}.
#' @export
neighbours <- function(SignVector){
  SVcolsums <- colSums(SignVector==1L)
  SVorder <- order(SVcolsums,decreasing=TRUE)
  SVcolsums_od <- SVcolsums[SVorder]
  SVordered <- SignVector[,SVorder]
  elim <- logical(nc <- ncol(SignVector))
  for(i in 1:nc){
    if(!elim[i]){
      count.i <- SVcolsums_od[i]
      tem_big <- SVordered[,i] - SVordered[,SVcolsums_od-count.i==1,drop=FALSE]
      tem_small <- SVordered[,i] - SVordered[,SVcolsums_od-count.i==-1,drop=FALSE]
      nbr_big <- which(colSums(abs(tem_big))==2L)
      if(!is.na(ind_b <- match(count.i+1,SVcolsums_od))) nbr_big <- nbr_big + ind_b - 1
      nbr_small <- which(colSums(abs(tem_small))==2L) + match(count.i-1,SVcolsums_od) - 1
      if(i > 1 & length(nbr_big)) elim[i] <- count.i < max(SVcolsums_od[nbr_big])
      if(i < nc & length(nbr_small)) elim[nbr_small[SVcolsums_od[nbr_small] < count.i]] <- TRUE
    }
  }
  elim[SVorder] <- elim
  !elim
}
