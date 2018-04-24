#' @title normalize ExpData
#' @description .AAA
#' @author Claus Weinholdt
#' @usage normalizeExoData(pval)
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

#' @title set classes a,b,c and d
#' @description set classes a,b,c and d
#' @author Claus Weinholdt
#' @usage get.Lset ()
#' @return a \code{list} with classes 
#' @export
get.Lset <- function(){
  Lsets <- list()
  Lsets[["a"]] <- list("s0"=c(1,2,3,4,5,6),"s1"=c())
  Lsets[["b"]] <- list("s0"=c(1,3,5)      ,"s1"=c(2,4,6))
  Lsets[["c"]] <- list("s0"=c(1,3,5,4)    ,"s1"=c(2,6))
  Lsets[["d"]] <- list("s0"=c(1,3,5,4,6)  ,"s1"=c(2))
  
  return(Lsets)
}

.Lset.get.var <- function(d,Lset,LsetMean,LsetN){
  #### unbiased sample variance using N-1 ! 
  
  if(LsetN['s1']>0){
    tmpvar <- (  sum((d[ Lset[["s0"]]  ] - LsetMean['s0'])^2) 
                 + sum((d[ Lset[["s1"]]  ] - LsetMean['s1'])^2)) / ((LsetN['s0']-1) + (LsetN['s1']-1))
  }else{
    tmpvar <- (sum((d[ Lset[["s0"]]  ] - LsetMean['s0'])^2))  / (LsetN['s0']-1)
  }
  return( unname(tmpvar) )
}

.Lset.get.logL <- function(d,Lset,LsetMean,LsetN,LsetVar){
  
  .NV <- function(d,mu,tau){
    # tau is precision
    # f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))  
    return(-0.5*tau*(d-mu)^2)
  }
  
  TAU <- 1/LsetVar
  
  if(LsetN['s1']>0){
    TmplogL <- -1*(sum(LsetN)/2) * log(LsetVar*2*pi) + sum( .NV(d[ Lset[["s0"]] ], LsetMean['s0'] ,TAU) , .NV(d[ Lset[["s1"]] ],LsetMean['s1'],TAU) )      
  }else{
    TmplogL <- -1*(sum(LsetN)/2) * log(LsetVar*2*pi) + sum( .NV(d[ Lset[["s0"]] ], LsetMean['s0'],TAU ) )    
  }
  
  return(TmplogL)  
}

#' @title calucate the the log likelihood for each class for each gene 
#' @description calucate the the log likelihood for each class for each gene 
#' @author Claus Weinholdt
#' @usage logLikelihoodOfnormData(normData$E)
#' @param x is maxtrix with normalized expression
#' @return a \code{matrix} with log likelihoods 
#' @export
logLikelihoodOfnormData <- function(x){  # with var-estimator

  Lsets <- get.Lset()
  
  ### init global varibales 
  ALL.MUs  <<- c()
  ALL.VARs <<- c()
  
  logL <- t(apply(x,1,function(d){
    res <- lapply(Lsets , function(Lset){ 
      LsetN <- sapply(Lset,length )
      LsetMean <- sapply(Lset,function(x) mean(d[x],na.rm = T)  )  #;if(is.nan(LsetMean['s1'])){ print("A")}
      LsetVar <- .Lset.get.var(d,Lset,LsetMean,LsetN)
      LsetlogL <- .Lset.get.logL(d,Lset,LsetMean,LsetN,LsetVar)
      return(list("N"=LsetN , 'Mean' = LsetMean, "Var" = LsetVar , "logL" = LsetlogL))
    })
    
    ALL.MUs <<- rbind( ALL.MUs, as.vector(sapply(res,function(x) x$Mean)) )
    ALL.VARs <<- rbind( ALL.VARs, as.vector(sapply(res,function(x) x$Var)) )
    
    return( sapply(res,function(x) x$logL))
  }))
  
  dimnames( ALL.MUs ) <- list(rownames(x),paste0(rep(names(Lsets),each=2),c('0','1')))
  dimnames( ALL.VARs ) <- list(rownames(x),names(Lsets))
  
  return(logL)
}

#' @title maxIndex of a vector
#' @description get the maximal Index
#' @author Claus Weinholdt
#' @usage maxIndex(x)
#' @param x is a vector
#' @return a \code{value}  
#' @export
maxIndex<-function(x){
  return(which(x == max(x)))
}

#' @title minIndex of a vector
#' @description get the minimal Index
#' @author Claus Weinholdt
#' @usage minIndex(x)
#' @param x is a vector
#' @return a \code{value}  
#' @export
minIndex<-function(x){
  return(which(x == min(x)))
}

#' @title calucate BIC and AIC
#' @description calucate BIC and AIC
#' @author Claus Weinholdt
#' @usage lInformationCriterion(logL,npar,k , IC='BIC')
#' @param logL is log liklelihood
#' @param npar represents the number of parameters in the fitted model
#' @param k is number of ovservations 
#' @param IC is BIC or AIC 
#' @return a \code{value}  
#' @export
InformationCriterion <- function(logL,npar,k , IC='BIC'){
  
  if(length(logL) != length(npar) || length(logL) != length(k)){
    print('npar not for each class') 
  }else{
    if(IC=='AIC'){
      #AIC== -2*loglikelihood + k*npar, npar represents the number of parameters in the fitted model, and k = 2 
      aic <- (-2)*logL + k*npar
      return(aic)
    }else if(IC=='BIC'){
      #BIC == AIC(object, ..., k = log(nobs(object)))
      
      k=log(k) # k -> number of ovservations (6)
      #npar represents the number of parameters
      
      bic <- (-2)*logL + npar *k 
      #alternative
      #bic<- (-2)*logL + npar *(k+log(2*pi))
      return(bic)
    }else{
      print(c(k,"<2 -- ERROR"))
    }
    
  }  
}

#' @title get InformationCriterion
#' @description get InformationCriterion (AIC or BIC) of LogLik matrix 
#' @author Claus Weinholdt
#' @usage get.IC(L,npar,k,IC)
#' @param L is log liklelihood matrix
#' @param npar represents the number of parameters in the fitted model
#' @param k is number of ovservations 
#' @param IC is BIC or AIC 
#' @return a \code{list}  
#' @export
get.IC <- function(L,npar,k,IC = 'BIC'){
  IC <- t(apply(L,1,function(x) InformationCriterion(x,npar,k,IC = IC) ))
  return(IC)
} 

#' @title get Posterior from BIC
#' @description get Posterior from BIC
#' @author Claus Weinholdt
#' @usage get.Posterior(B,Pis)
#' @param B is BIC matrix
#' @param Pis class prior
#' @return a \code{matrix} of  Posterior 
#' @export
get.Posterior <- function(B, Pis= c(0.7,0.1,0.1,0.1) ) {
  .qB<-function(B){ return(exp(-B/2)) }
  P_x_m <- t(log(apply(B,MARGIN=1,FUN = .qB)))
  
  .w <- function(B,pis){
    t <- c()
    for(i in 1:4){
      t[i] <- B[i]*pis[i]
    }
    return(t/sum(t))
  }
  pis=Pis  #c(0.90,0.01,0.08,0.01)
  P_m_x <- t(apply(exp(P_x_m),MARGIN=1,FUN=.w,Pis   ))
  
  colnames(P_m_x) <- colnames(B)
  
  return(P_m_x)
}

#' @title get Results of Posterior data
#' @description get the Genes assigned to the best class
#' @author Claus Weinholdt
#' @usage get.gene.classes(data, indexing="max",filter=0.75, DoPlot=FALSE)
#' @param data is Posterior matrix
#' @param indexing is the indexing method "max" or "min"
#' @param filter for Posterior 
#' @param DoPlot if TRUE plot histogram of Posterior
#' @return a \code{list} of assigned genes
#' @export
get.gene.classes <- function(data, indexing="max", filter=0.75, DoPlot=FALSE){
  
  Ind <- switch(indexing, 
                "max" =  apply(data,MARGIN=1,FUN=maxIndex),
                "min" =  apply(data,MARGIN=1,FUN=minIndex)
                )
  for(i in min(Ind):max(Ind) ) Ind[Ind==i]  <- colnames(data)[i]
  
  res <- lapply(colnames(data),function(x) { data[names(Ind)[Ind==x],]  })
  names(res) <- colnames(data)

  resFilter <- lapply(names(res) , function(x) res[[x]][ res[[x]][,x]>=filter ,  ])  
  names(resFilter) <- colnames(data)
  
  if(DoPlot){
    par(mfrow=c(2,2))
    for(i in names(res)) {hist(res[[i]][,i],50,main = i,xlab = "") ; abline(v=filter,col=2)}
    par(mfrow=c(1,1))
  }
  print(rbind("ALL"=sapply(res, nrow),"Filter"=sapply(resFilter, nrow) ))
  
  return(list("res"=res,"resFilter"=resFilter))
}

#' @title make data for enrichment analysis
#' @description make data for enrichment analysis
#' @author Claus Weinholdt
#' @param data is Posterior matrix
#' @param LsetForFC is the Lset of the data. needed for the logFC calculation 
#' @param normData is the normData object 
#' @return a \code{list} with logFC and Symbols
#' @export
make.enrichment.data <- function(data,LsetForFC,normData){

  data(Illumina_to_REFSEQ_MRNA, envir = environment()) 
  
  tmp <- normData$genes
  tmp <- tmp[rownames(data),]

  if(is.null(LsetForFC[['s1']])){
    logFC <- NA
  }else if( length(LsetForFC[['s0']]) == 1 ){
    logFC <- rowMeans(normData$E[rownames(data),LsetForFC[['s1']] ]) - (normData$E[rownames(data),LsetForFC[['s0']] ])
  }else if( length(LsetForFC[['s1']]) == 1 ){
    logFC <- (normData$E[rownames(data),LsetForFC[['s1']] ]) - rowMeans(normData$E[rownames(data),LsetForFC[['s0']] ])
  }else{
    logFC <- rowMeans(normData$E[rownames(data),LsetForFC[['s1']] ]) - rowMeans(normData$E[rownames(data),LsetForFC[['s0']] ])
  }
  tmp$logFC <- logFC
  tmp <- tmp[,c('SYMBOL','logFC')]
  
  tmp2 <- data.frame("REFSEQ_MRNA"=Illumina_to_REFSEQ_MRNA[rownames(tmp),]$REFSEQ_MRNA , "logFC" = tmp$logFC)
  rownames(tmp2) <- rownames(tmp)
  
  return( list( "SYMlogFC"=tmp,"SYM"= unique(as.character(tmp$SYMBOL)),
                "REFSEQlogFC"=tmp2,"REFSEQ"= unique(as.character(tmp2$REFSEQ_MRNA))))
}

#' @title Convert REFSEQ_MRNA string from DAVID to Illumina ids
#' @description Convert REFSEQ_MRNA string from DAVID to Illumina ids
#' @author Claus Weinholdt
#' @param string with REFSEQ_MRNA ids example:  "NM_032169, NM_000903, NM_015913"
#' @param asString if TRUE return a comma seperagted string with Illumina Ids
#' @return a \code{vector} with Illumina Ids
#' @export
REFSEQ_MRNA_to_Illm <- function(string, asString=FALSE){
  # example : string <-  "NM_032169, NM_000903, NM_015913"
  data(Illumina_to_REFSEQ_MRNA, envir = environment()) 
  setkey(Illumina_to_REFSEQ_MRNA,'REFSEQ_MRNA')
  REFSEQ <- stringr::str_replace( stringr::str_split(string,',')[[1]] , pattern = ' ', replacement = "")
  ILM <- Illumina_to_REFSEQ_MRNA[REFSEQ ,][['illumina_humanht_12_v4']]
  
  if(asString){
    ILM <- paste0(ILM,collapse = ',')
  }
  return(ILM )
}

PISStudy <- function(){
  ### 
  PIS <- list("P91"=c(0.91 ,0.03 ,0.03  ,0.03),
              "P85"=c(0.85 ,0.05 ,0.05  ,0.05),
              "P70"=c(0.7  ,0.1  ,0.1   ,0.1),
              "P55"=c(0.55 ,0.15 ,0.15  ,0.15),
              "P40"=c(0.4  ,0.2  ,0.2   ,0.2),
              "P25"=c(0.25 ,0.25 ,0.25  ,0.25)
  )
  PostClassPis  <- lapply(PIS, function(Pis){  
    tmp <- get.Posterior( normDataBIC ,Pis) 
    print(Pis)
    get.gene.classes(data=tmp,indexing="max",filter=0.75, DoPlot=TRUE) 
  })
  
}

main <- function(){
  
  # normalize Data
  normData <- normalizeExpData()

  # calulate logLik for class  Data
  Lsets <- get.Lset()
  normDataLogLik <- logLikelihoodOfnormData(normData$E)

  # calulate BIC from logLik 
  npar <- sapply(Lsets, function(x) sum(!sapply(x,is.null ) )) + 1  ## number parapeters for LogLilk -> mean + var 
  k <-  sapply(Lsets, function(x) sum( sapply(x,length) ))
  normDataBIC <- get.IC(normDataLogLik , npar, k , IC = 'BIC')
  BICminInd <- apply( normDataBIC,MARGIN = 1, FUN = minIndex)
  print(table(BICminInd))
  
  # calulate Posterior from BIC
  normDataPosterior <- get.Posterior( normDataBIC ,Pis = c(0.7,0.1,0.1,0.1))
  POSTmaxInd <- apply(normDataPosterior, MARGIN = 1 ,FUN = maxIndex)
  print(table(POSTmaxInd))
  
  PostClass <- get.gene.classes(data=normDataPosterior,indexing = "max",filter = 0.75, DoPlot = TRUE)
  
  
  ### enrichment analysis
  enrichmentFolder  <- '/Users/weinhol/GitHub/BGSC/inst/extdata/'
  
  ### do not run
  enrich <- list()
  for(i in names(PostClass[['resFilter']])){
    enrich[[i]] <-  make.enrichment.data(data = PostClass[['resFilter']][[i]] , LsetForFC = Lsets[[i]] , normData = normData )
    write.csv2(enrich[[i]]$REFSEQlogFC, paste0(enrichmentFolder,'Filtering_logFC_',i,'.csv'))
    write.table(enrich[[i]]$SYM, paste0(enrichmentFolder,'Filtering_Symbole_',i,'.txt'),quote = FALSE,row.names = FALSE,col.names = FALSE)
    write.table(enrich[[i]]$REFSEQ, paste0(enrichmentFolder,'Filtering_REFSEQ_MRNA_',i,'.txt'),quote = FALSE,row.names = FALSE,col.names = FALSE)
  }

  ChartReport <- list()
  ChartReportCategoryList  <- list()
  ChartReportCategorySigList <- list()
  ChartReportCategoryTermsList <- list()
  for(i in names(PostClass[['resFilter']])){
    ChartReport[[i]] <- fread( system.file("extdata",paste0("Filtering_REFSEQ_MRNA_",i,"_chartReport.txt"), package = "BGSC",mustWork = TRUE) )
    setkey(ChartReport[[i]],'Category')
  
    ChartReportCategorys <- unique(as.character(ChartReport[[i]]$Category))
    for(ChartReportCategory in ChartReportCategorys){
      
      ChartReportCategoryList
      ChartReportCategoryList[[ChartReportCategory]][[i]]$Genes <- sapply(ChartReportCategoryList[[ChartReportCategory]][[i]]$Genes, function(x) REFSEQ_MRNA_to_Illm(x,asString = TRUE))
      
      ChartReportCategoryList[[ChartReportCategory]][[i]] <- ChartReport[[i]][ChartReportCategory,]
      setkey(ChartReportCategoryList[[ChartReportCategory]][[i]],'Term')

      ChartReportCategorySigList[[ChartReportCategory]][[i]] <-ChartReportCategoryList[[ChartReportCategory]][[i]][which(Benjamini < 0.1),]
      setkey(ChartReportCategorySigList[[ChartReportCategory]][[i]],'Term')
      
      ChartReportCategoryTermsList[[ChartReportCategory]] <- unique(c(ChartReportCategoryTermsList[[ChartReportCategory]], ChartReportCategorySigList[[ChartReportCategory]][[i]]$Term) )
      
    }
  }
  
  ChartReportCategorySigMatList <- list()
  # ChartReportCategory <- 'KEGG_PATHWAY' #ChartReportCategory <- 'GOTERM_MF_FAT'
  for(ChartReportCategory in names(ChartReportCategorySigList)){
    print(ChartReportCategory)
    ChartReportCategorySigMat <- matrix(0, length(names(ChartReportCategorySigList[[ChartReportCategory]])), length(ChartReportCategoryTermsList[[ChartReportCategory]]),
                                        dimnames = list(names(ChartReportCategorySigList[[ChartReportCategory]]),ChartReportCategoryTermsList[[ChartReportCategory]]) )
    print(dim(ChartReportCategorySigMat))
    if(dim(ChartReportCategorySigMat)[2] > 0){
      for(i in names(ChartReportCategorySigList[[ChartReportCategory]]) ){
        tmp <- ChartReportCategorySigList[[ChartReportCategory]][[i]]$Term  
        if(length(tmp) > 0){
          ChartReportCategorySigMat[i,tmp] <- 1
        }
      }
      
      pheatmap::pheatmap(ChartReportCategorySigMat,cluster_cols = FALSE,cluster_rows = FALSE,legend_breaks = c(0,1),color = c('gray',2),main = ChartReportCategory,
                         gaps_row = c(1:nrow(ChartReportCategorySigMat) ),gaps_col = c(1:ncol(ChartReportCategorySigMat) ),border_color = 'black'
                         ,filename = paste0(enrichmentFolder,'/ChartReportSig_',ChartReportCategory,'.pdf'),width = 25,height = 10
                         )
      
      
      
      
    }  
  }
  
  

  library("org.Hs.eg.db")
  ENSEMBL2EG <- as.list(org.Hs.egENSEMBL2EG)
  ## Bimap interface:
  x <- org.Hs.egGO
  # Get the entrez gene identifiers that are mapped to a GO ID
  mapped_genes <- mappedkeys(x)
  # Convert to a list
  xx <- as.list(x[mapped_genes])
  if(length(xx) > 0) {
    # Try the first one
    got <- xx[[1]]           
    got[[1]][["GOID"]]
    got[[1]][["Ontology"]]
    got[[1]][["Evidence"]]
  }  
  EGFR.GOs <- names(xx[[ENSEMBL2EG[['ENSG00000146648']]]])
  
  EGFR_overlapping_GO <- list()
  for(i in names(PostClass[['resFilter']])){
    ChartReportCategoryList$GOTERM_MF_FAT[[i]]$GOTerm <- sapply(stringr::str_split(ChartReportCategoryList$GOTERM_MF_FAT[[i]]$Term,'~'),function(x) x[1])
    tmp <- ChartReportCategoryList$GOTERM_MF_FAT[[i]][ChartReportCategoryList$GOTERM_MF_FAT[[i]]$GOTerm %in% EGFR.GOs,]
    
    EGFR_overlapping_GO[[i]] <- rbind(EGFR_overlapping_GO[[i]],tmp )
    
  }  
  RBPs::f.input4(EGFR_overlapping_GO[['a']]$Term,
  EGFR_overlapping_GO[['b']]$Term,
  EGFR_overlapping_GO[['c']]$Term,
  EGFR_overlapping_GO[['d']]$Term,vennOut = T)
 
  tmp <- REFSEQ_MRNA_to_Illm(EGFR_overlapping_GO[['c']][EGFR_overlapping_GO[['c']]$Term == "GO:0005524~ATP binding",]$Genes)
  
  normData$E[intersect(rownames(normData$E),tmp),]
  
  table(POSTmaxInd[intersect(rownames(normData$E),tmp)])
  
    sapply(ChartReportCategoryList[[ChartReportCategory]][[i]]$Genes, function(x) REFSEQ_MRNA_to_Illm(x,asString = TRUE))
   
}
