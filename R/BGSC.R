#' @title normalize ExpData
#' @description .AAA
#' @author Claus Weinholdt
#' @usage normalizeExoData()
#' @param pval is the dectecion p value ot the Illumina BeatChip
#' @return a \code{limma object} 
#' @export
normalizeExpData <- function(pval=0.05){
  data(ExpData, envir = environment()) 
  pe2 <- propexpr(ExpData);
  dim(pe2) <- c(2,3);
  dimnames(pe2) <- list(CellType=c("Control","EGF"),Donor=c(1,2,3));
  print(pe2) ;
  print(mean(pe2[1:2,]))
  print(ExpData)

  print(dim(ExpData))
  y2 <- neqc(ExpData,negctrl="NEGATIVE") ; NORM<-"Limma_BG_QN"     #,robust=TRUE)
  expressed <- rowSums(y2$other$Detection <= pval) == 6 #>= 3
  #expressed <- rowSums(y2$other$Detection <= pval) >= 3
  y2 <- y2[expressed,]
  print(dim(y2))
  return(y2)
  
}

main <- function(){
  
  normData <- normalizeExpData()
}