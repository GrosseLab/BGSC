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
#' @usage logLikelihoodOfnormData(normData$E)
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
  BICminInd <- apply( normDataBIC,MARGIN=1,FUN=minIndex)
  print(table(BICminInd))
  
  # calulate Posterior from BIC
  normDataPosterior <- get.Posterior( normDataBIC ,Pis = c(0.7,0.1,0.1,0.1))
  POSTmaxInd <- apply(normDataPosterior,MARGIN=1,FUN=maxIndex)
  print(table(POSTmaxInd))
    
  ###
  t<-normDataPosterior[POSTmaxInd==1,] ;normDataPosterior_a<-t[order(t[,1],decreasing=T),] 
  t<-normDataPosterior[POSTmaxInd==2,] ;normDataPosterior_b<-t[order(t[,2],decreasing=T),]
  t<-normDataPosterior[POSTmaxInd==3,] ;normDataPosterior_c<-t[order(t[,3],decreasing=T),]
  t<-normDataPosterior[POSTmaxInd==4,] ;normDataPosterior_d<-t[order(t[,4],decreasing=T),]
  normDataPosterior_a.75<-normDataPosterior_a[normDataPosterior_a[,1]>0.75,]
  normDataPosterior_b.75<-normDataPosterior_b[normDataPosterior_b[,2]>0.75,]
  normDataPosterior_c.75<-normDataPosterior_c[normDataPosterior_c[,3]>0.75,]
  normDataPosterior_d.75<-normDataPosterior_d[normDataPosterior_d[,4]>0.75,]
  rbind('ALL'=c(
    'A'=dim(normDataPosterior_a)[1] ,
    'B'=dim(normDataPosterior_b)[1] ,
    'C'=dim(normDataPosterior_c)[1],
    'D'=dim(normDataPosterior_d)[1]) ,
    'x75'=c(
      dim(normDataPosterior_a.75)[1] ,
      dim(normDataPosterior_b.75)[1] ,
      dim(normDataPosterior_c.75)[1] ,
      dim(normDataPosterior_d.75)[1] )
  )
  
  
  ###
  PPP <- list("P91"=c(0.91 ,0.03 ,0.03  ,0.03),
              "P85"=c(0.85 ,0.05 ,0.05  ,0.05),
              "P70"=c(0.7  ,0.1  ,0.1   ,0.1),
              "P55"=c(0.55 ,0.15 ,0.15  ,0.15),
              "P40"=c(0.4  ,0.2  ,0.2   ,0.2),
              "P25"=c(0.25 ,0.25 ,0.25  ,0.25)
              )
  PPPtab <- lapply(PPP  , function(Pis){  tmp <- get.Posterior( normDataBIC ,Pis)
         tmpmaxInd <- apply(tmp,MARGIN=1,FUN=maxIndex)
         print(table(tmpmaxInd)) 
         
         t<-tmp[tmpmaxInd==1,] ;tmp_a<-t[order(t[,1],decreasing=T),] 
         t<-tmp[tmpmaxInd==2,] ;tmp_b<-t[order(t[,2],decreasing=T),]
         t<-tmp[tmpmaxInd==3,] ;tmp_c<-t[order(t[,3],decreasing=T),]
         t<-tmp[tmpmaxInd==4,] ;tmp_d<-t[order(t[,4],decreasing=T),]
         tmp_a.75<-tmp_a[tmp_a[,1]>0.75,]
         tmp_b.75<-tmp_b[tmp_b[,2]>0.75,]
         tmp_c.75<-tmp_c[tmp_c[,3]>0.75,]
         tmp_d.75<-tmp_d[tmp_d[,4]>0.75,]
         
         rbind('ALL'=c(
           'A'=dim(tmp_a)[1] ,
           'B'=dim(tmp_b)[1] ,
           'C'=dim(tmp_c)[1],
           'D'=dim(tmp_d)[1]) ,
           'x75'=c(
             dim(tmp_a.75)[1] ,
             dim(tmp_b.75)[1] ,
             dim(tmp_c.75)[1] ,
             dim(tmp_d.75)[1] )
         )
         # rownames(tmp_c.75)
         
  })
  # print( PPPtab)
  # RBPs::f.input4(PPPtab$P85,PPPtab$P70,PPPtab$P55,PPPtab$P25)
  
}
