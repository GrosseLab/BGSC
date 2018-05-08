# ----------------------------------------------------------------------
#' @title getQPCR
#' @description get qPCR
#' @author Claus Weinholdt
#' @return a \code{data.table} 
#' @export
getQPCR <- function(){   
  ### load data ---------------------------
    data(RawQPCRsf, envir = environment()) 
    qgenes <- c("CKAP2L", "ROCK1", "TPR","ALDH4A1", "CLCA2","GALNS") 
  
    # RawQPCRsf <- RawQPCR[RawQPCR$Probe == 'SF',]
    RawQPCRsf.rep <- list()
    RawQPCRsf.rep[['rep1']] <- RawQPCRsf[1:6,]
    RawQPCRsf.rep[['rep2']] <- RawQPCRsf[7:12,]
    RawQPCRsf.rep[['rep3']] <- RawQPCRsf[13:18,]

  ### parameter ---------------------------
    c0 <- get.Lset()$c$s0 ; c1 <- get.Lset()$c$s1
    # with M sample size of c0 and N sample size of c1
    M <- length(c0) ;N <- length(c1)
    Nrep <- length(RawQPCRsf.rep)
    doPlot <- FALSE
    
  ### delta Ct ---------------------------  
    delta.ct <-  lapply(qgenes,function(qg)  do.call(cbind, lapply(RawQPCRsf.rep , function(rep) rep[, colnames(rep)[grepl(qg,toupper(colnames(rep)))][1]  ]   ) ) )
    names(delta.ct) <- qgenes
    relative.Werte <-  lapply(qgenes,function(qg)  do.call(cbind, lapply(RawQPCRsf.rep , function(rep) rep[, colnames(rep)[grepl(qg,toupper(colnames(rep)))][2]  ]   ) ) )
    names(relative.Werte) <- qgenes
    
    delta.ct.Mean <- t(do.call(cbind,lapply(delta.ct, function(mat) rowMeans((mat)))))
    colnames(delta.ct.Mean) <- c(1:6)
    delta.ct.Mean.logFC <- (rowMeans(delta.ct.Mean[,c1]) - rowMeans(delta.ct.Mean[,c0]) ) * -1 
  
  ### relative.Werte ---------------------------    
    relative.Werte.geoMean <- t(do.call(cbind,lapply(relative.Werte, function(mat) rowMeans(log2(mat)))))
    colnames(relative.Werte.geoMean) <- c(1:6)
    relative.Werte.geoMean.C1 <- rowMeans(relative.Werte.geoMean[,c1])
    relative.Werte.geoMean.C0 <- rowMeans(relative.Werte.geoMean[,c0])
    relative.Werte.geoMean.C1C0.logFC <- relative.Werte.geoMean.C1 - relative.Werte.geoMean.C0
  
    if(doPlot){  
      barplot(t(cbind(relative.Werte.geoMean.C1C0.logFC,delta.ct.Mean.logFC)) ,beside=T,las=2,legend.text = c("rel",'deltaCt'),ylab = 'log2 fold change')   # qGenes$'GLIPR2' wahrscheinlich mit Fehler bei GLIPR2relative.Werte 
    }
    
  ### pooeld.var ---------------------------
    relative.Werte.var <- t(do.call(cbind,lapply(relative.Werte, function(mat) apply( log2(mat) ,1, function(vec) var(vec)  ) )))
    colnames(relative.Werte.var) <- c(1:6)
    relative.Werte.unvar <- relative.Werte.var * (Nrep - 1 ) # to get  sum_i ( x_i - mean(x) )^2  == `repUnVAR`
    relative.Werte.pooeld.var.C0 <- rowSums( relative.Werte.unvar[,c0] ) / ( sum( rep( (Nrep - 1), M) ) ) # ( sum_m  `repUnVAR`_m )  \ ( sum_m M-1 )
    relative.Werte.pooeld.var.C1 <- rowSums( relative.Werte.unvar[,c1] ) / ( sum( rep( (Nrep - 1), N) ) )
    relative.Werte.pooeld.var.C0C1_SatterthwaiteApproximation <- sqrt( (relative.Werte.pooeld.var.C0 / M) + (relative.Werte.pooeld.var.C1 / N) )
  
  ### var of all values ---------------------------
    relative.Werte.var.C0 <- sapply(relative.Werte, function(mat){ rownames(mat) <- c(1:6);  var(as.vector(log2(mat)[c0,]))} )
    relative.Werte.var.C1 <- sapply(relative.Werte, function(mat){ rownames(mat) <- c(1:6);  var(as.vector(log2(mat)[c1,]))} )
    relative.Werte.var.C0C1_SatterthwaiteApproximation <- sqrt( (relative.Werte.var.C0 / M) + (relative.Werte.var.C1 / N) )

  ### std.err of means --> same method as for Illumina  ---------------------------
    s1stderr <- apply(relative.Werte.geoMean[,c1],MARGIN = 1, function(row) plotrix::std.error(row,na.rm = T) )
    s0stderr <- apply(relative.Werte.geoMean[,c0],MARGIN = 1, function(row) plotrix::std.error(row,na.rm = T) )
    s0s1stderr = sqrt( s1stderr^2 +  s0stderr^2)
  
    if(doPlot){
      barplot(t(cbind(relative.Werte.var.C0 ,relative.Werte.pooeld.var.C0,s0stderr)) ,beside=T,las=2,legend.text = c("var",'pooled.var','stdErr'),ylab = 'var',main='C0')
      barplot(t(cbind(relative.Werte.var.C1 ,relative.Werte.pooeld.var.C1,s0stderr)) ,beside=T,las=2,legend.text = c("var",'pooled.var','stdErr'),ylab = 'var',main='C1')
      barplot(t(cbind(relative.Werte.var.C0C1_SatterthwaiteApproximation , relative.Werte.pooeld.var.C0C1_SatterthwaiteApproximation,s0s1stderr)) ,beside=T,las=2,legend.text = c("var",'pooled.var','stdErr'),ylab = 'var',main='C0C1_SatterthwaiteApproximation')
    }
  ### return---------------------------
    qCPRdataC <- data.table("rn"=names(relative.Werte) ,relative.Werte.geoMean.C1,relative.Werte.geoMean.C0,relative.Werte.geoMean.C1C0.logFC,
                             relative.Werte.pooeld.var.C0,relative.Werte.pooeld.var.C1,relative.Werte.pooeld.var.C0C1_SatterthwaiteApproximation,
                             relative.Werte.var.C0,relative.Werte.var.C1,relative.Werte.var.C0C1_SatterthwaiteApproximation,
                            s1stderr,s0stderr,s0s1stderr)
    setkey(qCPRdataC)              
    return(qCPRdataC)
  
  # vec <- c( 2,33,45,12,44,21,442)
  # sqrt( sum( (vec - mean(vec))^2 ) / (length(vec)-1) )
  # sd(vec)
  # sum( (vec - mean(vec))^2 ) 
  # var(vec) * (length(vec)-1)
}  

# ---------------------------------------------------------------------- 
#' @title normalize ExpData
#' @description .AAA
#' @author Claus Weinholdt
#' @usage normalizeExpData(DetectionPval = 0.05 , DetectionPvalNumber = "ALL")
#' @aliases normalizeExpData
#' @param DetectionPval is the dectecion p value ot the Illumina BeatChip
#' @param DetectionPvalNumber miniml number of samples with DetectionPval. If DetectionPvalNumber is "ALL" then DetectionPvalNumber equal to the number of samples
#' @return a \code{limma object} 
#' @import limma stringr
#' @export
normalizeExpData <- function(DetectionPval = 0.05 , DetectionPvalNumber = "ALL"){
  data(ExpData, envir = environment()) 
  pe2 <- limma::propexpr(ExpData);
  dim(pe2) <- c(2,3);
  dimnames(pe2) <- list(CellType=c("Control","EGF"),Donor=c(1,2,3));
  # print(pe2) ;
  # print(mean(pe2[1:2,]))
  # print(ExpData)

  qqq <- normexp.fit.detection.p(ExpData, detection.p="Detection")
  
  y2 <- limma::neqc(ExpData,negctrl = "NEGATIVE",detection.p = ExpData$other$Detection ,robust=TRUE,regular = "regular", offset = 16) ; NORM <- "Limma_BG_QN"     #,robust=TRUE)
  head(y2$E)
  
  if (DetectionPvalNumber == "ALL") {
    
    DetectionPvalNumber <- ncol(y2$other$Detection) 
    
    expressed <- rowSums(y2$other$Detection <= DetectionPval) == DetectionPvalNumber #>= 3
  
  } else {
    if( DetectionPvalNumber >  ncol(y2$other$Detection) ) DetectionPvalNumber <- ncol(y2$other$Detection) 
    
    expressed <- rowSums(y2$other$Detection <= DetectionPval) >= DetectionPvalNumber
   
  }
  
 
  y2 <- y2[expressed,]
  
  
  message(
    paste0('genes with detection pval <= ',DetectionPval,' in ',
           DetectionPvalNumber ,' of ',ncol(y2$other$Detection),' Samples --> ',
    sum(stringr::str_count(rownames(y2$E),'ILMN_')), ' of ',
    sum(stringr::str_count(rownames(ExpData$E),'ILMN_'))
  ))
  
  
  return(y2)
  

  
}

# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
.Lset.get.var <- function(d,Lset,LsetMean,LsetN){
  #### unbiased sample variance using N-1 ! 
  
  if (LsetN['s1'] > 0) {
    tmpvar <- ( sum((d[ Lset[["s0"]]  ] - LsetMean['s0'])^2) 
              + sum((d[ Lset[["s1"]]  ] - LsetMean['s1'])^2)) / ((LsetN['s0']-1) + (LsetN['s1']-1))
  }else{
    tmpvar <- (sum((d[ Lset[["s0"]]  ] - LsetMean['s0'])^2))  / (LsetN['s0']-1)
  }
  return( unname(tmpvar) )
}

# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
#' @title calucate the the log likelihood for each class for each gene 
#' @description calucate the the log likelihood for each class for each gene 
#' @author Claus Weinholdt
#' @usage logLikelihoodOfnormData(x)
#' @param x is maxtrix with normalized expression e.g. normData$E
#' @return a \code{list} with log likelihoods , all means and all variances 
#' @export
logLikelihoodOfnormData <- function(x){  # with var-estimator

  Lsets <- get.Lset()
  
  tmpRes <- apply(x,1,function(d){
    res <- lapply(Lsets , function(Lset){ 
      LsetN <- sapply(Lset,length )
      LsetMean <- sapply(Lset,function(x) mean(d[x],na.rm = T)  )  #;if(is.nan(LsetMean['s1'])){ print("A")}
      LsetVar <- .Lset.get.var(d,Lset,LsetMean,LsetN)
      LsetlogL <- .Lset.get.logL(d,Lset,LsetMean,LsetN,LsetVar)
      return(list("N" = LsetN , 'Mean' = LsetMean, "Var" = LsetVar , "logL" = LsetlogL))
    })
    
    ALL.MUs <- as.vector(sapply(res,function(x) x$Mean)) 
    ALL.VARs <- as.vector(sapply(res,function(x) x$Var)) 
    logL <-  as.vector(sapply(res,function(x) x$logL ))
           
    return( list("logL" = logL,"ALL.MUs" = ALL.MUs, "ALL.VARs" = ALL.VARs) )
  })
  
  logL <- do.call(rbind,lapply(tmpRes,function(x) x[["logL"]]) )
  ALL.MUs <- do.call(rbind,lapply(tmpRes,function(x) x[["ALL.MUs"]]) )
  ALL.VARs <- do.call(rbind,lapply(tmpRes,function(x) x[["ALL.VARs"]]) )
  
  colnames( logL ) <- names(Lsets)
  colnames( ALL.MUs ) <-  paste0(rep(names(Lsets),each = 2),c('0','1'))
  colnames( ALL.VARs ) <- names(Lsets)
  
  # print(data.table::data.table(logL,keep.rownames = T))
  # print(data.table::data.table(ALL.MUs,keep.rownames = T))
  # print(data.table::data.table(ALL.VARs,keep.rownames = T))
   
  return(list("logL" = logL,"ALL.MUs" = ALL.MUs, "ALL.VARs" = ALL.VARs) )
}

# ----------------------------------------------------------------------
#' @title maxIndex of a vector
#' @description get the maximal Index
#' @author Claus Weinholdt
#' @usage maxIndex(x)
#' @param x is a vector
#' @return a \code{value}  
#' @export
maxIndex <- function(x){
  return(which(x == max(x)))
}

# ----------------------------------------------------------------------
#' @title minIndex of a vector
#' @description get the minimal Index
#' @author Claus Weinholdt
#' @usage minIndex(x)
#' @param x is a vector
#' @return a \code{value}  
#' @export
minIndex <- function(x){
  return(which(x == min(x)))
}

# ----------------------------------------------------------------------
#' @title calucate BIC and AIC
#' @description calucate BIC and AIC
#' @author Claus Weinholdt
#' @param logL is log liklelihood
#' @param npar represents the number of parameters in the fitted model
#' @param k is number of ovservations 
#' @param IC is BIC or AIC 
#' @return a \code{value}  
#' @export
InformationCriterion <- function(logL,npar,k , IC = 'BIC'){
  
  # x = the observed data;
  # k = the number of data points in x, the number of observations, or equivalently, the sample size;
  # npar = the number of free parameters to be estimated. 
  # p(x|npar) = the probability of the observed data given the number of parameters; or, the likelihood of the parameters
  # given the dataset;
  # logL = log of the maximized value of the likelihood function for the estimated model.
  
  if (length(logL) != length(npar) || length(logL) != length(k)) {
    print('npar not for each class') 
  } else {
    if (IC == 'AIC') {
      #AIC== -2*loglikelihood + k*npar, npar represents the number of parameters in the fitted model, and k = 2 
      aic <- (-2)*logL + k*npar
      return(aic)
    } else if (IC == 'BIC') {
      #BIC == AIC(object, ..., k = log(nobs(object)))
      
      k <- log(k) # k -> number of ovservations (6)
      #npar represents the number of parameters
      
      bic <- (-2)*logL + npar * k 
      #alternative
      #bic<- (-2)*logL + npar *(k+log(2*pi))
      return(bic)
    } else {
      print(c(k,"<2 -- ERROR"))
    }
    
  }  
}

# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
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

# ----------------------------------------------------------------------
#' @title get Results of Posterior data
#' @description get the Genes assigned to the best class
#' @author Claus Weinholdt
#' @usage get.gene.group(data, indexing="maximal",filter=0.75, DoPlot=FALSE)
#' @param data is Posterior matrix
#' @param indexing is the indexing method "maximal" or "minimal"
#' @param filter for Posterior 
#' @param DoPlot if TRUE plot histogram of Posterior
#' @import graphics
#' @return a \code{list} of assigned genes
#' @export
get.gene.group <- function(data, indexing="maximal", filter=0.75, DoPlot=FALSE){
  
  Ind <- switch(indexing, 
                "maximal" =  apply(data,MARGIN=1,FUN=maxIndex),
                "minimal" =  apply(data,MARGIN=1,FUN=minIndex)
                )
  for(i in min(Ind):max(Ind) ) Ind[Ind == i]  <- colnames(data)[i]
  
  res <- lapply(colnames(data),function(x) { data[names(Ind)[Ind == x],]  })
  names(res) <- colnames(data)

  resFilter <- lapply(names(res) , function(x) res[[x]][ res[[x]][,x] >= filter ,  ])  
  names(resFilter) <- colnames(data)
  
  if(DoPlot){
    print("Histograms of approximate posterior")
    par(mfrow=c(2,2))
    for(i in names(res)) { hist(res[[i]][,i],50,main = paste0('group ',i),xlab = "approximate posterior") ; abline(v = filter,col = 2)}
    par(mfrow=c(1,1))
  }
  tmp <- rbind("#genes assigned to group" = sapply(res, nrow),"Filter" = sapply(resFilter, nrow) )
  rownames(tmp)[2] <- paste0("#genes assigned to group with Filter",filter)  
  print(paste0('Number of genes assigned to group with the ',indexing,' approximate posterior'))
  print(tmp)
  
  return(list("res" = res,"resFilter" = resFilter))
}

# ----------------------------------------------------------------------
#' @title normalize ExpData
#' @description .AAA
#' @author Claus Weinholdt
#' @usage get.log2Mean.and.log2FC(normData)
#' @param normData is the normData data object
#' @return a \code{list} with mean , std.err , log2FC and std.err of log2FC 
#' @note  \code{std.error.xy = sqrt( std.error(x)^2 +  std.error(y)^2)}
#' @import purrr plotrix
#' @export
get.log2Mean.and.log2FC <- function(normData) {
  
  Lset <- get.Lset()
  resOut <- purrr::map2(Lset,names(Lset), 
              function(set, setN){
                
                tmp.dt <- data.table::data.table()
                s0 <- set$s0
                s1 <- set$s1
                
                if ( !is.null(s1) ) {

                  if ( (length(s0) > 1 & length(s1) > 1) ) {
                    
                    s1M <- rowMeans( normData$E[,s1]) 
                    s0M <- rowMeans( normData$E[,s0]) 
                    s1s0FC <- s1M - s0M
                    
                    s1stderr <- apply(normData$E[,s1],MARGIN = 1, function(row) plotrix::std.error(row,na.rm = T) )
                    s0stderr <- apply(normData$E[,s0],MARGIN = 1, function(row) plotrix::std.error(row,na.rm = T) )
                    
                    # https://www.researchgate.net/post/Can_anyone_help_with_calculating_error_in_RT-qPCRs_fold-change_data
                    s0s1stderr = sqrt( s1stderr^2 +  s0stderr^2)
                    
                    tmp.dt <- data.table::data.table("rn" = names(s1M), s1M , s0M , s1s0FC , s1stderr, s0stderr, s0s1stderr)                      
                    data.table::setkey(tmp.dt,'rn')
                    
                  } else if ( (length(s0) == 1 & length(s1) > 1) ) {
                    
                    s1M <- rowMeans( normData$E[,s1]) 
                    s0M <- normData$E[,s0]
                    s1s0FC <- s1M - s0M
                    
                    tmp.dt <- data.table::data.table("rn" = names(s1M), s1M , s0M , s1s0FC)                  
                    data.table::setkey(tmp.dt,'rn')
                    
                  } else if ( (length(s1) == 1 & length(s0) > 1) ) { 
                    
                    s1M <- normData$E[,s1] 
                    s0M <- rowMeans( normData$E[,s0]) 
                    s1s0FC <- s1M - s0M
                    
                    tmp.dt <- data.table::data.table("rn" = names(s1M), s1M , s0M , s1s0FC)                    
                    data.table::setkey(tmp.dt,'rn')
                    
                  }
                  
                } else {
                  s0M <- rowMeans( normData$E[,s0]) 
                  
                  tmp.dt <- data.table::data.table("rn" = names(s0M) , s0M )                    
                  data.table::setkey(tmp.dt,'rn')
                }
                
                return(tmp.dt)
                
              } )
  
  ### SatterthwaiteApproximation -> sqrt( (var(x)/M ) + (var(y)/N) ) ; with M sample size of and N sample size of y
  ### equal to 
  ### std.error.xy = sqrt( std.error(x)^2 +  std.error(y)^2 )
  FYI <- FALSE
  if (FYI) {
    c0 <- c(1,3,4,5);      c1 <- c(2,6); 
    xkm <- normData$E[,c0] ; M <- length(c0)
    ykn <- normData$E[,c1] ; N <- length(c1)
    #http://www.statisticshowto.com/satterthwaite-approximation/
    #https://wolfweb.unr.edu/~ldyer/classes/396/PSE.pdf
    SatterthwaiteApproximation <- sqrt( (apply(xkm,1,var) / M) + (apply(ykn,1,var) / N) )
    hist( SatterthwaiteApproximation - resOut$c[names(SatterthwaiteApproximation),][['s0s1stderr']],100)
  }
  
  return(resOut)
  
}

# ----------------------------------------------------------------------
#' @title Density plot for NV fit 
#' @description Density plot for NV fit 
#' @author Claus Weinholdt
#' @param id is a Illuminan id
#' @param normData is the normData data object
#' @param ALL.MUs matrix with all means from logLikelihoodOfnormData()
#' @param ALL.VARs matrix with all variances from logLikelihoodOfnormData()
#' @param useGroup if set to a group (eg. "a","b","c" or "d") the row is colored
#' @param DOplot if TRUE plot is printed
#' @param basesize size of font
#' @param GrAblack if TRUE coloring group 'a' in black because we do not use a Indicator variable
#' @param onlySYMBOL if TRUE title is only gene SYMBOL
#' @return a \code{list} with mean , std.err , log2FC and std.err of log2FC 
#' @import ggplot2 gridExtra grid stats
#' @export
Density.NV.fit.plot <- function(id, normData, ALL.MUs, ALL.VARs, useGroup = NA, DOplot = FALSE, basesize = 14, GrAblack = TRUE , onlySYMBOL = FALSE ){
  Lset <- get.Lset()
  indicatorTMP <- lapply(Lset,function(x){
    vec <- rep(0,6);names(vec) <- c(1:6)
    if (!is.null(x$s1)) { vec[x$s1] <- 1 }
    return(vec)
  })
  
  IDs.dt <- data.table::data.table(normData$genes,keep.rownames = T,key = 'rn')
  
  pg.M <- ALL.MUs[id,]; names(pg.M) <- colnames(ALL.MUs)
  pg.s <- sqrt(ALL.VARs[id,]); names(pg.s) <- colnames(ALL.VARs)
  
  logE = normData$E[id,]
  tmp = data.frame(logE=rep(logE,4),
                   density = rep( 0, length(rep(logE,4))),
                   group    = c( rep('a',6),rep('b',6),rep('c',6),rep('d',6) ),
                   # indicator = factor(c( rep(0,6) , c(0,1,0,1,0,1), c(0,1,0,0,0,1) ,c(0,1,0,0,0,0) ))
                   indicator = factor( c(indicatorTMP$a , indicatorTMP$b , indicatorTMP$c ,indicatorTMP$d ))
  )
  
  x <- seq(4, 20, length=1000)
  dM <- round(max(
    max(dnorm(x,mean =  pg.M['a0'], sd = pg.s['a'])),
    max(dnorm(x,mean =  pg.M['b1'], sd = pg.s['b'])),
    max(dnorm(x,mean =  pg.M['c1'], sd = pg.s['c'])),
    max(dnorm(x,mean =  pg.M['d1'], sd = pg.s['d'])) )+0.5) 
  
  col2= RColorBrewer::brewer.pal(8,'Paired')[c(6,2)] 
  
  .thememap <- function(base_size = 12, legend_key_size = 0.4, base_family = "", col = "grey70") {
    ggplot2::theme_gray(base_size = base_size, base_family = base_family) %+replace% 
      ggplot2::theme(title = ggplot2::element_text(face="bold", colour=1,angle=0           ,vjust= 0.0,           size=base_size),
                     axis.title.x = ggplot2::element_text(face="bold", colour=1, angle=0   ,vjust= 0.0,           size=base_size),
                     # axis.text.x  = ggplot2::element_text(face="bold", colour=1, angle = -30 , vjust = 1, hjust = 0, size=base_size),
                     axis.text.x  = ggplot2::element_text(face="bold", colour=1,  size=base_size),
                     strip.text.x = ggplot2::element_text(face="bold", colour=1, angle=0    ,vjust= 0.5,           size=base_size),
                     axis.title.y = ggplot2::element_text(face="bold", colour=1, angle=90   ,vjust= 1.5, hjust=.5, size=base_size),
                     axis.text.y  = ggplot2::element_text(face="bold", colour=1,                                  size=base_size),
                     legend.text  = ggplot2::element_text(face="bold" ,colour=1, angle=0  ,vjust= 0.0,             size=base_size),
                     legend.title = ggplot2::element_text(face="bold" ,colour=1, angle=0  ,vjust= 0.2,             size=base_size),
                     
                     
                     axis.ticks =  ggplot2::element_line(colour = "grey70", size = 0.5),
                     panel.grid.major =  ggplot2::element_line(colour = col, size = 0.2),
                     panel.grid.minor =  ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill="white",size = 0.2,),
                     # #panel.grid.minor.y = element_line(size=3),
                     # panel.grid.major = ggplot2::element_line(colour = "white"),
                     
                     # Force the plot into a square aspect ratio
                     # aspect.ratio = 1,
                     # aspect.ratio = 9 / 16,
                     
                     legend.key.size  = ggplot2::unit(legend_key_size, "cm"),
                     plot.title = ggplot2::element_text(hjust = 0.5 , vjust= 1)          
      )
  }
  
  ymax <- 0
  if( max(logE,na.rm = T) > 0 )   ymax = max(logE,na.rm = T) else ymax = 1  
  ymin <- -0
  if( min(logE,na.rm = T) > 0)   ymin = min(logE,na.rm = T) else ymin = 0  
  Lims <- c( floor(ymin) - 2  , ceiling(ymax) + 2 )
  
  
  if (onlySYMBOL) {
    labstitle=paste0(IDs.dt[id,][['SYMBOL']])
  } else {    
    labstitle=paste0(id,' -- ',IDs.dt[id,][['SYMBOL']])
  }  
  
  g2 = ggplot(tmp, aes(x = logE,y = density,colour = indicator)) +
    scale_y_continuous(limits = c(-0.5, dM),breaks = seq(0,dM,1)) +
    # scale_x_continuous(limits = c(4.5, 15),breaks = seq(4.5,15,1)) +
    scale_x_continuous(limits = Lims ,breaks = seq(Lims[1],Lims[2],1)) +
    geom_point(shape = 20,size = 0)  + 
    facet_wrap(~ group,nrow = 4)  +
    scale_color_manual(values = col2 , name = "Indicator variable", labels = list("g = 0", "g = 1" )) +
    labs(title = labstitle ,y = "Density", x = "Logarithmic expression levels") +
    .thememap(base_size = basesize,0.6) +
    theme(legend.background = element_rect(fill = "grey90", size=.5, linetype = "dotted")) +
    theme(legend.position = "bottom") 
  
  if (!is.na(useGroup)) {
    bestSet <- subset(tmp, group == useGroup)
    g2 <- g2 + geom_rect(data = bestSet, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), colour = NA , fill = "yellow", alpha = .01)
  }
  
  g3 <- g2 + 
    with(tmp[tmp$group == "a",],stat_function(data = tmp[tmp$group == "a" & tmp$indicator == 0,],fun = dnorm,args = list(mean =  pg.M['a0'], sd = pg.s['a']), size = 1.2) ) +
    with(tmp[tmp$group == "b",],stat_function(data = tmp[tmp$group == "b" & tmp$indicator == 0,],fun = dnorm,args = list(mean =  pg.M['b0'], sd = pg.s['b']), size = 1.2) ) +
    with(tmp[tmp$group == "b",],stat_function(data = tmp[tmp$group == "b" & tmp$indicator == 1,],fun = dnorm,args = list(mean =  pg.M['b1'], sd = pg.s['b']), size = 1.2) ) +
    with(tmp[tmp$group == "c",],stat_function(data = tmp[tmp$group == "c" & tmp$indicator == 0,],fun = dnorm,args = list(mean =  pg.M['c0'], sd = pg.s['c']), size = 1.2) ) +
    with(tmp[tmp$group == "c",],stat_function(data = tmp[tmp$group == "c" & tmp$indicator == 1,],fun = dnorm,args = list(mean =  pg.M['c1'], sd = pg.s['c']), size = 1.2) ) +
    with(tmp[tmp$group == "d",],stat_function(data = tmp[tmp$group == "d" & tmp$indicator == 0,],fun = dnorm,args = list(mean =  pg.M['d0'], sd = pg.s['d']), size = 1.2) ) +
    with(tmp[tmp$group == "d",],stat_function(data = tmp[tmp$group == "d" & tmp$indicator == 1,],fun = dnorm,args = list(mean =  pg.M['d1'], sd = pg.s['d']), size = 1.2) )
  
  g4 <- g3 + 
    geom_point(data = tmp, aes(x = logE,y = density),shape = 20,size = 3.5) +
    geom_point(data = tmp, aes(x = logE,y = density),shape = 21,size = 3,colour = "black")
  
  ### coloring group 'a' in black because we do not use a Indicator variable 
  if (GrAblack) {
    SetA <- subset(tmp, group == 'a')
    g4 <- g4 + 
      with(tmp[tmp$group == "a",],stat_function(data = tmp[tmp$group == "a" & tmp$indicator == 0,],fun = dnorm,args = list(mean =  pg.M['a0'], sd = pg.s['a']), size = 1.2 ,colour = "darkslategray") ) +
      geom_point(data = SetA, aes(x = logE,y = density), shape = 20,size = 3.5,colour = "slategray") +
      geom_point(data = SetA, aes(x = logE,y = density), shape = 21,size = 3 , colour = "black")
  }
  
  if (DOplot) { print(g4) }
  
  colors()[grep("gray",colors())]
  
  # g4 + annotate("rect", xmin=Lims[1], xmax=Lims[2], ymin=0, ymax=Inf, alpha=0.2, fill="red") 
  
  return(g4)
}

# ----------------------------------------------------------------------
#' @title make.plot.data.FC.Ill.qPCR
#' @description make.plot.data.FC.Ill.qPCR
#' @author Claus Weinholdt
#' @param qCPRdata is data.table with qPCR data
#' @param qgenesIDs Ids of qCPR genes
#' @param MeanFoldChangeClass is data.table with Illumina data
#' @param class for the plots 
#' @return a \code{data.frame} as plot input 
#' @export
make.plot.data.FC.Ill.qPCR <- function(qCPRdata, qgenesIDs ,MeanFoldChangeClass, class = "c"){
  TMP <- data.frame()
  for (tmpG in qCPRdata$rn  ) {
    expC <- MeanFoldChangeClass[[class]][qgenesIDs[[tmpG]], ]
    TMPexp <- data.frame( "Gene" = tmpG,
                          "ExpID" = expC$rn,
                          # "Set"="Illumina",
                          "Set" = "Microarray",
                          "FC" = expC[["s1s0FC"]] ,
                          "stderr" = expC[["s0s1stderr"]]
    )
    
    TMPqpc <- data.frame( "Gene" = tmpG,
                          # "Set" = "RT-qPCR",
                          "Set" = "qPCR",
                          "FC" = qCPRdata[tmpG,][["relative.Werte.geoMean.C1C0.logFC"]] ,
                          #"stderr" = qCPRdata[tmpG,][["relative.Werte.pooeld.var.C0C1_SatterthwaiteApproximation"]],
                          "stderr" = qCPRdata[tmpG,][["s0s1stderr"]]
    )
    
    for (i in 1:nrow(TMPexp)) {
      tmp <- TMPqpc 
      tmp$ExpID <-  as.character(TMPexp[i,]$ExpID)
      TMP <- rbind(TMP , rbind(TMPexp[i, ],tmp[,  colnames(TMPexp) ]))
    }
  }  
  TMP$pid <- paste0(TMP$Gene,'::',TMP$ExpID)
  TMP$pid <- factor(TMP$pid)
  return(TMP)
  
}   

# ----------------------------------------------------------------------
#' @title make.plot.data.exp.Ill.qPCR
#' @description make.plot.data.exp.Ill.qPCR
#' @author Claus Weinholdt
#' @param qCPRdata is data.table with qPCR data
#' @param MeanFoldChangeClass is data.table with Illumina data
#' @param class for the plots 
#' @return a \code{data.frame} as plot input 
#' @export
make.plot.data.exp.Ill.qPCR <- function(qCPRdata,MeanFoldChangeClass, class="c"){
  TMP <- data.frame()
  for (tmpG in qCPRdata$rn  ) {
    expC <- MeanFoldChangeClass[[class]][qgenesIDs[[tmpG]], ]
    TMPexp1 <- data.frame("Gene" = tmpG,
                          "ExpID" = expC$rn,
                          'Set' = paste0(class,'1'),
                          "Mean" = expC$s1M ,
                          "stderr" = expC$s1stderr)
    TMPexp0 <- data.frame("Gene" = tmpG,
                          "ExpID" = expC$rn,
                          'Set' = paste0(class,'0'),
                          "Mean" = expC$s0M ,
                          "stderr" = expC$s0s1stderr)
    
    TMP <-  rbind(TMP,rbind(TMPexp0,TMPexp1) )
    
  }  
  # TMP <- TMP[ TMP$Gene != 'KIF5C',]
  TMP$pid <- paste0(TMP$Gene,'::',TMP$ExpID)
  TMP$pid <- factor(TMP$pid)
  return(TMP)
  
}

# ----------------------------------------------------------------------
#' @title ggplot thememap for barplot
#' @description ggplot thememap for barplot
#' @author Claus Weinholdt
#' @param base_size size of font
#' @param legend_key_size size of key for the legend 
#' @param base_family base family
#' @param colgridmajor color of major grid 
#' @return a \code{theme} 
#' @export
thememapBarplot <- function(base_size = 12, legend_key_size = 0.4, base_family = "", colgridmajor = "grey70") {
  ggplot2::theme_gray(base_size = base_size, base_family = base_family) %+replace% 
    ggplot2::theme(title = ggplot2::element_text(face="bold", colour=1,angle=0           ,vjust= 0.0,           size=base_size),
                   axis.title.x = ggplot2::element_text(face="bold", colour=1, angle=0   ,vjust= 0.0,           size=base_size),
                   # axis.text.x  = ggplot2::element_text(face="bold", colour=1, angle = -30 , vjust = 1, hjust = 0, size=base_size),
                   axis.text.x  = ggplot2::element_text(face="bold", colour=1, angle = 270          , hjust = 0, size=base_size),
                   strip.text.x = ggplot2::element_text(face="bold", colour=1, angle=0    ,vjust= 0.5,           size=base_size),
                   axis.title.y = ggplot2::element_text(face="bold", colour=1, angle=90   ,vjust= 1.5, hjust=.5, size=base_size),
                   axis.text.y  = ggplot2::element_text(face="bold", colour=1,                                  size=base_size),
                   legend.text  = ggplot2::element_text(face="bold" ,colour=1, angle=0  ,vjust= 0.0,             size=base_size),
                   legend.title = ggplot2::element_text(face="bold" ,colour=1, angle=0  ,vjust= 0.2,             size=base_size),
                   
                   
                   axis.ticks =  ggplot2::element_line(colour = "grey70", size = 0.5),
                   panel.grid.major =  ggplot2::element_line(colour = colgridmajor, size = 0.2),
                   panel.grid.minor =  ggplot2::element_blank(),
                   panel.background = ggplot2::element_rect(fill="white",size = 0.2,),
                   # #panel.grid.minor.y = element_line(size=3),
                   # panel.grid.major = ggplot2::element_line(colour = "white"),
                   
                   # Force the plot into a square aspect ratio
                   # aspect.ratio = 1,
                   # aspect.ratio = 9 / 16,
                   
                   legend.key.size  = ggplot2::unit(legend_key_size, "cm"),
                   plot.title = ggplot2::element_text(hjust = 0.5 , vjust= 1)          
    )
}

# ----------------------------------------------------------------------
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
  } else if( length(LsetForFC[['s0']]) == 1 ) {
    logFC <- rowMeans(normData$E[rownames(data),LsetForFC[['s1']] ]) - (normData$E[rownames(data),LsetForFC[['s0']] ])
  } else if( length(LsetForFC[['s1']]) == 1 ) {
    logFC <- (normData$E[rownames(data),LsetForFC[['s1']] ]) - rowMeans(normData$E[rownames(data),LsetForFC[['s0']] ])
  } else {
    logFC <- rowMeans(normData$E[rownames(data),LsetForFC[['s1']] ]) - rowMeans(normData$E[rownames(data),LsetForFC[['s0']] ])
  }
  tmp$logFC <- logFC
  tmp <- tmp[,c('SYMBOL','logFC')]
  
  tmp2 <- data.frame("REFSEQ_MRNA" = Illumina_to_REFSEQ_MRNA[rownames(tmp),]$REFSEQ_MRNA , "logFC" = tmp$logFC)
  rownames(tmp2) <- rownames(tmp)
  
  return( list( "SYMlogFC" = tmp, "SYM" = unique(as.character(tmp$SYMBOL)),
                "REFSEQlogFC" = tmp2, "REFSEQ" = unique(as.character(tmp2$REFSEQ_MRNA)))
          )
}

# ----------------------------------------------------------------------
#' @title Convert REFSEQ_MRNA string from DAVID to Illumina ids
#' @description Convert REFSEQ_MRNA string from DAVID to Illumina ids
#' @author Claus Weinholdt
#' @param string with REFSEQ_MRNA ids example:  "NM_032169, NM_000903, NM_015913"
#' @param asString if TRUE return a comma seperagted string with Illumina Ids
#' @return a \code{vector} with Illumina Ids
#' @import stringr
#' @export
REFSEQ_MRNA_to_Illm <- function(string, asString=FALSE){
  # example : string <-  "NM_032169, NM_000903, NM_015913"
  data(Illumina_to_REFSEQ_MRNA, envir = environment()) 
  data.table::setkey(Illumina_to_REFSEQ_MRNA,'REFSEQ_MRNA')
  REFSEQ <- stringr::str_replace( stringr::str_split(string,',')[[1]] , pattern = ' ', replacement = "")
  ILM <- Illumina_to_REFSEQ_MRNA[REFSEQ ,][['illumina_humanht_12_v4']]
  
  if(asString){
    ILM <- paste0(ILM,collapse = ',')
  }
  return(ILM)
}

# study the effect of prior ----------------------------------------------------------------------
PIS_Study <- function(normDataBIC, qgenesIDs){
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
    get.gene.group(data=tmp,indexing="maximal",filter = 0.75, DoPlot=TRUE) 
    
  })
  
  print(lapply(PostClassPis,function(x) x$res$c[as.character(do.call(c,qgenesIDs)),]))
  tmp <- f.input.list(lapply(PostClassPis,function(x) rownames(x$resFilter$c) )[-1])
  tmp <- f.input.list(lapply(PostClassPis,function(x) rownames(x$resFilter$c) )[-6])
  
  return(PostClassPis)
}

# ----------------------------------------------------------------------
makeExpData <- function() {
  wd2 <- "/Volumes/ianvsITZroot/home/adsvy/Kappler/Kappler_Wichmann_Medizin"
  x <- limma::read.ilmn(files=paste0(wd2,"/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_nobkgd_SF767-1-799-6.txt")  ,sep="\t",
                        ctrlfiles=paste0(wd2,"/KW_120813_1343/Analysen_SF767-1-799-6/ControlProbeProfile_SF767-1-799-6.txt"),
                        other.columns=c("Detection","BEAD_STDERR","Avg_NBEADS")
  )
  
  ##### 2. ### reduce to our 6 samples !!!! 
  x2<-x
  targets<-c()
  targets$CellType=c("KO","KO_EGF","ALL","ALL_EGF","14","14_EGF")
  set<-c(1:6) #SF 767
  # set<-c(7:12) # x 999
  x2$E<-x2$E[,set]
  x2$other$Detection<-x2$other$Detection[,set]
  x2$other$BEAD_STDERR<-x2$other$BEAD_STDERR[,set]
  x2$other$Avg_NBEADS<-x2$other$Avg_NBEADS[,set]
  x2$targets<-targets
  
  ExpData <- x2
  setwd('/Users/weinhol/GitHub/BGSC')
  devtools::use_data(ExpData,pkg = 'BGSC')
  #data(ercc, envir = environment()) 
}

# ----------------------------------------------------------------------
#' @import utils
makeQPCR <- function(){
  RawQPCR <- read.csv(system.file("extdata", "RawQPCR.csv", package = "BGSC"),sep = ';')
  # setwd('/Users/weinhol/GitHub/BGSC')
  # devtools::use_data(RawQPCR)
  
  qgenes<-c("CKAP2L", "ROCK1", "TPR","ALDH4A1", "CLCA2","GALNS") 
  RawQPCR.paper <- RawQPCR[,c(colnames(RawQPCR)[1:2],as.vector(sapply(qgenes,function(x) colnames(RawQPCR)[grepl(x,toupper(colnames(RawQPCR)))] )))]
  RawQPCRsf <- RawQPCR.paper[ RawQPCR.paper$Probe == 'SF',]
  RawQPCR799 <- RawQPCR.paper[ RawQPCR.paper$Probe == '799',]
  setwd('/Users/weinhol/GitHub/BGSC')
  devtools::use_data(RawQPCRsf)
  
}