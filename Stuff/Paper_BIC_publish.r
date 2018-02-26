# on='/Volumes/sol.informatik.uni-halle.de'
# source(paste(sep="",on,"/Promotion/Kappler_Wichmann_Medizin/Auswertung/KW_functions.r"))
# source(paste(sep="",on,"/Promotion/Kappler_Wichmann_Medizin/Auswertung/MY_DE_Bayes.functions.r"))
# 

wd='/home/adsvy/Kappler/Kappler_Wichmann_Medizin/Auswertung'
# wd = '/Volumes/work512/home/weinhol/Promotion/Kappler_Wichmann_Medizin/Auswertung'
# wd='/home/weinhol/Promotion/Kappler_Wichmann_Medizin/Auswertung'

wd='/Volumes/ianvsITZ/Kappler/Kappler_Wichmann_Medizin/Auswertung'
setwd(wd)
source(paste0("KW_functions.r"))
source(paste0("MY_DE_Bayes.functions.r"))

library('grDevices')#,lib.loc ="/Library/Frameworks/R.framework/Versions/3.3/Resources/library" )

## save.image("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2016_9-14.RData")
##load("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2016_9-14.RData")

## save.image("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2016_9-22.RData")
load("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2016_9-22.RData")

# wd='/Users/weinhol/Promotion/ServerTMP/Kappler_Wichmann_Medizin/Auswertung'

a=c(42,55,24,13,44,30 )
boxplot(a)
par(mar=c(10,4,4,2) + 0.1,mfrow=c(1,1))

plot(density(rnorm(2000,mean(a),sd(a))) )
abline(v=mean(a))
dnorm(1000,mean(a),sd(a))

read.Data<-function(){
  # on<-"~";
  # on='/Volumes/sol.informatik.uni-halle.de'
  x <- read.ilmn(files=paste0(wd,"/../KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_nobkgd_SF767-1-799-6.txt")  ,sep="\t",
                     ctrlfiles=paste0(wd,"/../KW_120813_1343/Analysen_SF767-1-799-6/ControlProbeProfile_SF767-1-799-6.txt"),
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
  
  
  pe2 <- propexpr(x2);dim(pe2) <- c(2,3);dimnames(pe2) <- list(CellType=c("Control","EGF"),Donor=c(1,2,3));pe2 ;mean(pe2[1:2,])
  
  return(x2)
}

normalize.and.fiter.Data<-function(x2,pval=0.05){
  y2 <- neqc(x2,negctrl="NEGATIVE") ; NORM<-"Limma_BG_QN"     #,robust=TRUE)
  expressed <- rowSums(y2$other$Detection <= pval) == 6 #>= 3
  #expressed <- rowSums(y2$other$Detection <= pval) >= 3
  
  y2 <- y2[expressed,]
  print(dim(y2))
  return(y2)
}

ALL.MUs=c()
ALL.SDs=c()
get.ML.logLs.var.NEW<-function(x){  # with var-estimator

  NV<-function(d,mu,tau){
    #tau is precision
    # f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))  
    return(-0.5*tau*(d-mu)^2)
  }
  
  get.var<-function(d,sub,MUs,N,i0,i1){
    return(  ( sum((d[sub[[i0]]]-MUs[i0])^2) + sum((d[sub[[i1]]]-MUs[i1])^2)) / (  (N[i0]-1) + (N[i1]-1) ) ) 
  }
  
  calc.MUs<-function(d,sub){  
    MUs<-c()
    for(i in 1:length(sub)){
      MUs[i]<-mean(d[sub[[i]]])
    }
    return(MUs)
  }
  
  logLs<-function(d,sub,N){
    
    MUs<-calc.MUs(d,sub)
    
    ALL.MUs <<- rbind( ALL.MUs,MUs)
    
    VARs<-c()
    VARs[1]<- 1/(N[1]-1)*sum( (d[a0]-MUs[1])^2)
    VARs[2]<-get.var(d,sub,MUs,N,2,3)
    VARs[3]<-get.var(d,sub,MUs,N,4,5)
    VARs[4]<-get.var(d,sub,MUs,N,6,7)
    SDs<-sqrt(VARs)
    
    ALL.SDs <<- rbind( ALL.SDs,SDs)
    
    TAUs<-1/VARs
    
    logL<-c()   
    logL[1]<-( -(N[1]/2) * log(VARs[1]*2*pi) +(sum( (NV(d[sub[[1]]],MUs[1],TAUs[1]))  ) ) )    
    logL[2]<-( -(N[1]/2) * log(VARs[2]*2*pi) +(sum( (NV(d[sub[[2]]],MUs[2],TAUs[2])) , (NV(d[sub[[3]]],MUs[3],TAUs[2])) ) ) )    
    logL[3]<-( -(N[1]/2) * log(VARs[3]*2*pi) +(sum( (NV(d[sub[[4]]],MUs[4],TAUs[3])) , (NV(d[sub[[5]]],MUs[5],TAUs[3])) ) ) )    
    logL[4]<-( -(N[1]/2) * log(VARs[4]*2*pi) +(sum( (NV(d[sub[[6]]],MUs[6],TAUs[4])) , (NV(d[sub[[7]]],MUs[7],TAUs[4])) ) ) )
    
    alt.meth<-function(){   
      normprob = function (x,m,v) { (1/sqrt(2*pi*v))*exp(-((x-m)^2)/(2*v))  }   
         
      logL2<-c()  
        logL2[1] <- log(prod(c(normprob(d[sub[[1]]][1],MUs[1],VARs[1]),
                               normprob(d[sub[[1]]][2],MUs[1],VARs[1]),
                               normprob(d[sub[[1]]][3],MUs[1],VARs[1]),
                               normprob(d[sub[[1]]][4],MUs[1],VARs[1]),
                               normprob(d[sub[[1]]][5],MUs[1],VARs[1]),
                               normprob(d[sub[[1]]][6],MUs[1],VARs[1]))))
    
       logL2[2] <- log(prod(c(normprob(d[sub[[2]]][1],MUs[2],VARs[2]),
                              normprob(d[sub[[2]]][2],MUs[2],VARs[2]),
                              normprob(d[sub[[2]]][3],MUs[2],VARs[2]),
                              normprob(d[sub[[3]]][1],MUs[3],VARs[2]),
                              normprob(d[sub[[3]]][2],MUs[3],VARs[2]),
                              normprob(d[sub[[3]]][3],MUs[3],VARs[2]))))
                
      logL2[3] <- log(prod(c(normprob(d[sub[[4]]][1],MUs[4],VARs[3]),
                             normprob(d[sub[[4]]][2],MUs[4],VARs[3]),
                             normprob(d[sub[[4]]][3],MUs[4],VARs[3]),
                             normprob(d[sub[[4]]][4],MUs[4],VARs[3]),
                             normprob(d[sub[[5]]][1],MUs[5],VARs[3]),
                             normprob(d[sub[[5]]][2],MUs[5],VARs[3]))))
      
      logL2[4] <- log(prod(c(normprob(d[sub[[6]]][1],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][2],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][3],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][4],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][5],MUs[6],VARs[4]),
                             normprob(d[sub[[7]]][1],MUs[7],VARs[4]))))
      
      # ### Test for 
      # x=D.n$E
      # g="ILMN_1730999" #TPR
      # d=x[g,]
      # logL2
      # [1] -5.534718 -3.177679  1.601310 -4.269458
      # > logL
      # [1] -5.534718 -3.177679  1.601310 -4.269458
      # 
         
    }     
         
    get.single.NormProb<-function(){
      normprob = function (x,m,v) { (1/sqrt(2*pi*v))*exp(-((x-m)^2)/(2*v))  }   
      
      normprobS<-list()  
      normprobS[[1]] <- c(normprob(d[sub[[1]]][1],MUs[1],VARs[1]),
                             normprob(d[sub[[1]]][2],MUs[1],VARs[1]),
                             normprob(d[sub[[1]]][3],MUs[1],VARs[1]),
                             normprob(d[sub[[1]]][4],MUs[1],VARs[1]),
                             normprob(d[sub[[1]]][5],MUs[1],VARs[1]),
                             normprob(d[sub[[1]]][6],MUs[1],VARs[1]))

      normprobS[[2]] <- c(normprob(d[sub[[2]]][1],MUs[2],VARs[2]),
                             normprob(d[sub[[2]]][2],MUs[2],VARs[2]),
                             normprob(d[sub[[2]]][3],MUs[2],VARs[2]),
                             normprob(d[sub[[3]]][1],MUs[3],VARs[2]),
                             normprob(d[sub[[3]]][2],MUs[3],VARs[2]),
                             normprob(d[sub[[3]]][3],MUs[3],VARs[2]))
      
      normprobS[[3]] <-  c(normprob(d[sub[[4]]][1],MUs[4],VARs[3]),
                             normprob(d[sub[[4]]][2],MUs[4],VARs[3]),
                             normprob(d[sub[[4]]][3],MUs[4],VARs[3]),
                             normprob(d[sub[[4]]][4],MUs[4],VARs[3]),
                             normprob(d[sub[[5]]][1],MUs[5],VARs[3]),
                             normprob(d[sub[[5]]][2],MUs[5],VARs[3]))
      
      normprobS[[4]] <-  c(normprob(d[sub[[6]]][1],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][2],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][3],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][4],MUs[6],VARs[4]),
                             normprob(d[sub[[6]]][5],MUs[6],VARs[4]),
                             normprob(d[sub[[7]]][1],MUs[7],VARs[4]))
     
      SetNames=names(normprobS[[2]])
      normprobS[[1]] <- normprobS[[1]][SetNames]
      normprobS[[2]] <- normprobS[[2]][SetNames]
      normprobS[[3]] <- normprobS[[3]][SetNames]
      normprobS[[4]] <- normprobS[[4]][SetNames]
      
      normprobS.tab <- do.call(rbind,normprobS);rownames(normprobS.tab)=c("a","b","c","d")
  
      write.table(normprobS.tab,file=paste0(wd,'/TPR_normprobS.tab.txt'),row.names=T,col.names=NA,quote=F,sep="\t")
      saveRDS(normprobS.tab,file=paste0(wd,'/TPR_normprobS.tab.RDS'))
      apply(normprobS.tab, 1, prod)
  
    }
         
    return(logL)
  }
  
  N<-c()
  a0<-c(1,2,3,4,5,6);   N[1]  <- length(a0);
  b0<-c(1,3,5);         N[2] <- length(b0);
  b1<-c(2,4,6);         N[3] <- length(b1);
  c0<-c(1,3,5,4);       N[4] <- length(c0);
  c1<-c(2,6);           N[5] <- length(c1);
  d0<-c(1,3,5,4,6);     N[6] <- length(d0);
  d1<-c(2);             N[7] <- length(d1);
  sub<-list(a0=a0,b0=b0,b1=b1,c0=c0,c1=c1,d0=d0,d1=d1)
  
  ALL.MUs<<-c()
  ALL.SDs<<-c()
  require(reshape)
  logL<-apply(x,1,logLs,sub,N)
  
  restructur<-function(A){
    A<-t(A)
    colnames(A)<-c("a","b","c","d")
    return(A)
  }
  
  logL<-restructur(logL)  
  
  rownames(ALL.MUs) <-rownames(x)
  colnames(ALL.MUs) <-names(sub)
  rownames(ALL.SDs) <-rownames(x)
  colnames(ALL.SDs) <-c("a",'b','c','d')
  ALL.MUs <<-ALL.MUs
  ALL.SDs <<-ALL.SDs
  
  return(logL)
}

myIC<-function(logL,npar,k=2){
  if(k==2){
    #AIC== -2*loglikelihood + k*npar, npar represents the number of parameters in the fitted model, and k = 2 
    aic<- (-2)*logL + k*npar
    return(aic)
  }else if(k>2){
    #BIC == AIC(object, ..., k = log(nobs(object)))
    
    k=log(k) # k -> number of ovservations (hier 6)
    #npar represents the number of parameters
    
    bic<- (-2)*logL + npar *k 
    #alternative
    #bic<- (-2)*logL + npar *(k+log(2*pi))
    return(bic)
  }else{
    print(c(k,"<2 -- ERROR"))
  }
}

get.posterior<-function(B, Pis=c(0.90,0.01,0.08,0.01) ) {
  qB<-function(B){
    return(exp(-B/2))
  }
  
  P_x_m<-apply(B$IC,MARGIN=1,FUN = qB)
  P_x_m<-t(log(P_x_m))

    w<-function(B,pis){
      t<-c()
      for(i in 1:4){
        t[i]<-B[i]*pis[i]
      }
      return(t/sum(t))
    }
  pis=Pis  #c(0.90,0.01,0.08,0.01)
  P_m_x<-apply(exp(P_x_m),MARGIN=1,FUN=w,pis   )
  P_m_x<-t(P_m_x)

  return(P_m_x)
}

D<-read.Data()
D.n<-normalize.and.fiter.Data(D)
D.n_L<-get.ML.logLs.var.NEW(D.n$E);range(D.n_L)
D.n_L.BIC <- make.IC(D.n_L,c(2,3),6) # BIC with 6 = number of values   
D.n_L.BIC.Post<-get.posterior(D.n_L.BIC,Pis = c(0.7,0.1,0.1,0.1))
m<-apply((D.n_L.BIC.Post),MARGIN=1,FUN=maxIndex)
print(table(m))

t<-D.n_L.BIC.Post[m==1,] ;D.n_L.BIC.Post_a<-t[order(t[,1],decreasing=T),] 
t<-D.n_L.BIC.Post[m==2,] ;D.n_L.BIC.Post_b<-t[order(t[,2],decreasing=T),]
t<-D.n_L.BIC.Post[m==3,] ;D.n_L.BIC.Post_c<-t[order(t[,3],decreasing=T),]
t<-D.n_L.BIC.Post[m==4,] ;D.n_L.BIC.Post_d<-t[order(t[,4],decreasing=T),]
D.n_L.BIC.Post_a.85<-D.n_L.BIC.Post_a[D.n_L.BIC.Post_a[,1]>0.75,]
D.n_L.BIC.Post_b.85<-D.n_L.BIC.Post_b[D.n_L.BIC.Post_b[,2]>0.75,]
D.n_L.BIC.Post_c.85<-D.n_L.BIC.Post_c[D.n_L.BIC.Post_c[,3]>0.75,]
D.n_L.BIC.Post_d.85<-D.n_L.BIC.Post_d[D.n_L.BIC.Post_d[,4]>0.75,]
rbind('ALL'=c(
'A'=dim(D.n_L.BIC.Post_a)[1] ,
'B'=dim(D.n_L.BIC.Post_b)[1] ,
'C'=dim(D.n_L.BIC.Post_c)[1],
'D'=dim(D.n_L.BIC.Post_d)[1]) ,
'x75'=c(
dim(D.n_L.BIC.Post_a.85)[1] ,
dim(D.n_L.BIC.Post_b.85)[1] ,
dim(D.n_L.BIC.Post_c.85)[1] ,
dim(D.n_L.BIC.Post_d.85)[1] )
)


ILMN_SYMBOL.dt <- data.table(cbind(rownames(D$genes),D$genes$SYMBOL) );colnames(ILMN_SYMBOL.dt) <- c('ILMN','symbol');setkey(ILMN_SYMBOL.dt,'ILMN')
ILMN_SYMBOL.dt$symbol <- toupper(ILMN_SYMBOL.dt$symbol)

x75symb = c( length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_a.85),]$symbol))-1 + sum(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_a.85),]$symbol=='')  ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_b.85),]$symbol))-1 + sum(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_b.85),]$symbol=='')  ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_c.85),]$symbol))-1 + sum(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_c.85),]$symbol=='')  ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_d.85),]$symbol))-1 + sum(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_d.85),]$symbol=='')  )


ALLsymb = c( length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_a),]$symbol))-1   ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_b),]$symbol))-1   ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_c),]$symbol))-1   ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_d),]$symbol))-1  )
x75symb = c( length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_a.85),]$symbol))-1   ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_b.85),]$symbol))-1   ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_c.85),]$symbol))-1   ,
             length(unique(ILMN_SYMBOL.dt[rownames(D.n_L.BIC.Post_d.85),]$symbol))-1  )

tttt = round(D.n_L.BIC.Post_c.85[Paper.qgenesIDs,],3)
rownames(tttt) <- names(Paper.qgenesIDs)
Paper.qgenesIDs



# IDs<-cbind(rownames(D$genes),D$genes);colnames(IDs)<-c("ILMN","SYMBOL")

GO{
  
  library("org.Hs.eg.db")
  require("AnnotationDbi")
  require("GO.db")
  require("KEGG.db")
  require("annotate")
  # loading the org.At.tair.db annotation Package -> updated biannually
  require("org.At.tair.db")
  require("GOstats")
  library("AnnotationForge")
  require("GSEABase")
  
  require("GO.db")
  go_df <- data.frame(GOID=unlist(eapply(GOTERM, function(x) x@GOID)), Term=unlist(eapply(GOTERM, function(x) x@Term)), Ont=unlist(eapply(GOTERM, function(x) x@Ontology)))
  
  xx <- as.list(org.Hs.egGO2EG)
  if(length(xx) > 0){
    # Gets the entrez gene ids for the top 2nd and 3nd GO identifiers
    goids <- xx$'GO:0019903'
    
    goids <- xx$'GO:0097159'
    # Gets the entrez gene ids for the first element of goids
    setkey(EntrezID_ILM.dt,'entrez_id')
    tmp <- EntrezID_ILM.dt[as.integer(goids),]
    
    setkey(entrez.dt,ILMN)
    entrez.dt[as.character(tmp$ILMN),]
    
    genes
    
    
    goframeData.dt <- data.table(goframeData);colnames(goframeData.dt) <-c('GO','Evidence','entrez_id')
    setkey(goframeData.dt,'GO')
    expGo <- goframeData.dt['GO:0019903',]
    setkey(expGo,'entrez_id')
    
    length(unique(EntrezID_ILM.dt[expGo,]$symbol))
    sum(duplicated(universe))
    
    # Evidence code for the mappings
    names(goids[[1]])
  }
  # For org.Hs.egGO2ALLEGS
  xx <- as.list(org.Hs.egGO2ALLEGS)
  if(length(xx) > 0){
    # Gets the Entrez Gene identifiers for the top 2nd and 3nd GO identifiers
    goids <- xx[2:3]
    # Gets all the Entrez Gene identifiers for the first element of goids
    goids[[1]]
    # Evidence code for the mappings
    names(goids[[1]])
  }
  
  
  
  frame = toTable(org.Hs.egGO)
  goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
  head(goframeData)
  goFrame=GOFrame(goframeData,organism="Homo sapiens")
  goAllFrame=GOAllFrame(goFrame)
  gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
  
  KEGGframe = toTable(org.Hs.egPATH)
  KEGGframeData = data.frame(KEGGframe$path_id, KEGGframe$gene_id)
  head(KEGGframeData)
  KEGGFrame=KEGGFrame(KEGGframeData,organism="Homo sapiens")
  KEGGgsc <- GeneSetCollection(KEGGFrame, setType = KEGGCollection())
  
  universe = Lkeys(org.Hs.egGO) #EntrezID
  
'Symx75'= c(
  length( unique(D$genes[rownames(D.n_L.BIC.Post_a.85),"SYMBOL"]) ) ,
  length( unique(D$genes[rownames(D.n_L.BIC.Post_b.85),"SYMBOL"]) ),
  length( unique(D$genes[rownames(D.n_L.BIC.Post_c.85),"SYMBOL"]) ),
  length( unique(D$genes[rownames(D.n_L.BIC.Post_d.85),"SYMBOL"]) ) )

  EntrezID=fread(paste0(wd,"/../../Kappler_HIF1/TFBS/genes_EntrezID.txt"))
  EntrezID$symbol <- toupper(EntrezID$symbol)
  
  ILMN_SYMBOL.dt <- data.table(IDs[,c(1,2)]);colnames(ILMN_SYMBOL.dt)[2] <- 'symbol';
  ILMN_SYMBOL.dt$symbol <- toupper(ILMN_SYMBOL.dt$symbol)
  setkey(ILMN_SYMBOL.dt,'symbol')
  
  EntrezID_ILM.dt <- ILMN_SYMBOL.dt[EntrezID,]

  EGFR_BP=c("GO:0000165","GO:0000186","GO:0001503","GO:0001892","GO:0001934","GO:0001942","GO:0006950","GO:0007165","GO:0007166","GO:0007173","GO:0007202","GO:0007411","GO:0007435","GO:0007611","GO:0008283","GO:0008284","GO:0008543","GO:0016337","GO:0018108","GO:0021795","GO:0030335","GO:0031659","GO:0035413","GO:0038095","GO:0042059","GO:0042177","GO:0042327","GO:0043006","GO:0043066","GO:0043406","GO:0045087","GO:0045429","GO:0045739","GO:0045740","GO:0045944","GO:0046777","GO:0048011","GO:0048015","GO:0048146","GO:0048546","GO:0050679","GO:0050730","GO:0050999","GO:0051205","GO:0051897","GO:0060571","GO:0070141","GO:0070374","GO:0071364","GO:0071392")
  EGFR_CC=c("GO:0000139","GO:0005615","GO:0005634","GO:0005737","GO:0005768","GO:0005789","GO:0005886","GO:0010008","GO:0016020","GO:0016021","GO:0016323","GO:0030122","GO:0031965","GO:0045121","GO:0048471","GO:0070435")
  EGFR_MF=c("GO:0003682","GO:0003690","GO:0004709","GO:0004713","GO:0004714","GO:0004716","GO:0004888","GO:0005006","GO:0005515","GO:0005524","GO:0019899","GO:0019903","GO:0030235","GO:0042802","GO:0046982","GO:0051015")
  EGRF_kegg=c("01521","01522","04010","04012","04014","04015","04020","04060","04066","04068","04072","04144","04151","04320","04510","04520","04540","04810","04912","04915","04921","05120","05160","05200","05205","05206","05212","05213","05214","05215","05218","05219","05223","05224","05230","05231")
  
  SysDate="2016-12-19" ; #Sys.Date()
  ResFolder=paste0(wd,'/GOStat/',SysDate,'_')
  for(set in c('Gr_a','Gr_b','Gr_c','Gr_d')){
    # set='Gr_c'
    setkey(EntrezID_ILM.dt,'ILMN')  
     entrez.dt =  switch (set,
          'Gr_a' = EntrezID_ILM.dt[rownames(D.n_L.BIC.Post_a.85),],
          'Gr_b' = EntrezID_ILM.dt[rownames(D.n_L.BIC.Post_b.85),],
          'Gr_c' = EntrezID_ILM.dt[rownames(D.n_L.BIC.Post_c.85),],
          'Gr_d' = EntrezID_ILM.dt[rownames(D.n_L.BIC.Post_d.85),]
      )
   
     ### https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4012721/ 
    PMC4012721  = f.input2(toupper(c("Rab25",'Rab17','RIN2','SPRY1','DNM3','MUC1','SGK1','ECH1')),entrez.dt$symbol)$inter
    print(PMC4012721)
    entrez.dt[symbol %in% PMC4012721 ,]
    entrez.dt[symbol %in% EGFR_genes,]
    
    genes = entrez.dt$entrez_id
    genes = unique(as.character(genes[!is.na(genes)]))
    
    write.table( entrez.dt ,paste0(ResFolder,set,'_input.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
    write.table( entrez.dt[symbol %in% EGFR_genes,] ,paste0(ResFolder,set,'_EGFR_genes.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
    
    if(length(genes)>0){
      
      for(ONTO in c('BP', 'CC', 'MF')){
        # ONTO='BP' # ontology (BP, CC, MF)
        for(TD in c("over",'under')){
          # TD = "over" TD = 'under'
          print(paste0(ResFolder,set,'_',ONTO,'_',TD))
          
          CONI = FALSE #; CONI = TRUE
          params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",
                                       geneSetCollection=gsc,
                                       geneIds = genes,
                                       universeGeneIds = universe,
                                       ontology = ONTO, 
                                       pvalueCutoff = 1,
                                       conditional = CONI,
                                       testDirection = TD)
          print('Test')
          hyperGT <- hyperGTest(params)
          hyperGT.df = summary(hyperGT)
          
          hyperGT.df$BH <- p.adjust(hyperGT.df$Pvalue,method = 'BH')
          
          hyperGT.df <- hyperGT.df[,colnames(hyperGT.df)[c(1,2,8,3:7)] ]
          hyperGT.df.sig <- hyperGT.df[hyperGT.df[,'Pvalue']<=0.05,]
          
          print(as.data.table(hyperGT.df.sig))
          hyperGT.dt.sig <- as.data.table(hyperGT.df.sig);
          
          hyperGT.dt.sig=  switch (ONTO,
                             'BP' =  setkey(hyperGT.dt.sig,'GOBPID'),
                             'CC' =  setkey(hyperGT.dt.sig,'GOCCID'),
                             'MF' =  setkey(hyperGT.dt.sig,'GOMFID'),
          )
          
          ONTPO.dt=  switch (ONTO,
                         'BP' =  hyperGT.dt.sig[EGFR_BP,],
                         'CC' =  hyperGT.dt.sig[EGFR_CC,],
                         'MF' =  hyperGT.dt.sig[EGFR_MF,],
          )
          
          write.table( as.data.table(hyperGT.df)     ,paste0(ResFolder,set,'_',ONTO,'_',TD,'_all.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
          write.table( as.data.table(hyperGT.df.sig) ,paste0(ResFolder,set,'_',ONTO,'_',TD,'_sig.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
          write.table( ONTPO.dt ,paste0(ResFolder,set,'_',ONTO,'_',TD,'_EGFRterm.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
          
          htmlReport(hyperGT, file = paste0(ResFolder,set,'_',ONTO,'_',TD,'.html') )
          
          test{
            # subset the output table to get the columns of interest 
            # (GO ID, GO description, p-value)
            out <- subset(hyperGT.df.sig, select=c(1, 7, 2))
            # retrieve input genes associated with each GO identifier
            # use the org.Hs.eg data mapping to get GO terms for each ID
            goMaps <- lapply(out$GOBPID, function(x) unlist(mget(x, org.Hs.egGO2ALLEGS)))
            # subset the selected genes based on those in the mappings
            goSelected <- lapply(goMaps, function(x) selected[selected %in% x])
            # join together with a semicolon to make up the last column
            out$inGenes <- unlist(lapply(goSelected, function(x) paste(x, collapse=";")))
            # write the final data table as a tab separated file 
            write.table(out, file="go_results.tsv", sep="\t", row.names=FALSE)
            
          }
          
        }
      }     
      
      kparams <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params",
                                      geneSetCollection=KEGGgsc,
                                      geneIds = genes,
                                      universeGeneIds = universe,
                                      pvalueCutoff = 1,
                                      testDirection = "over")
      kOver <- hyperGTest(kparams)
      hyperGT.df = summary(kOver)
      
      hyperGT.df$BH <- p.adjust(hyperGT.df$Pvalue,method = 'BH')
      
      hyperGT.df <- hyperGT.df[,colnames(hyperGT.df)[c(1,2,8,3:7)] ]
      hyperGT.df.sig <- hyperGT.df[hyperGT.df[,'Pvalue']<=0.05,]
      
      hyperGT.dt.sig <- as.data.table(hyperGT.df.sig);setkey(hyperGT.dt.sig,'KEGGID')
      
      KEGGterm=hyperGT.dt.sig[EGRF_kegg,]
      print(as.data.table(hyperGT.df.sig))
      
      write.table( as.data.table(hyperGT.df)     ,paste0(ResFolder,set,'_KEGG_all.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
      write.table( as.data.table(hyperGT.df.sig) ,paste0(ResFolder,set,'_KEGG_sig.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
      write.table( KEGGterm ,paste0(ResFolder,set,'_',ONTO,'_',TD,'KEGG_EGFRterm.tsv'),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
      
      htmlReport(kOver, file = paste0(ResFolder,set,'_KEGG.html') )
      
    }else{
      print(paste0(ResFolder,set))
      print('Zero')
    }
    
  }  





}


a0<-c(1,2,3,4,5,6);   b0<-c(1,3,5);       b1<-c(2,4,6);         c0<-c(1,3,4,5);      c1<-c(2,6);         d0<-c(1,3,4,5,6);    d1<-c(2);            
addFC_log2<-function(RNs,c0,c1,D.n){
  
  if(length(c1)>1){
    FC<-2^(rowMeans(D.n$E[RNs,c1])-rowMeans(D.n$E[RNs,c0]))    
  }else{
    FC<-2^( D.n$E[RNs,c1]-rowMeans(D.n$E[RNs,c0]))          
  }
  FC<-cbind(2^(D.n$E[RNs,c0]),2^(D.n$E[RNs,c1]),FC,log2(FC))
  FC<-FC[order(abs(FC[,8]),decreasing=T),]
  colnames(FC)[8]<-"log2FC"
  return(FC)
}
FC.B = addFC_log2(RNs = rownames(D.n_L.BIC.Post_b),b0,b1,D.n)
FC.C = addFC_log2(RNs = rownames(D.n_L.BIC.Post_c),c0,c1,D.n)
FC.D = addFC_log2(RNs = rownames(D.n_L.BIC.Post_d),d0,d1,D.n)

Paper.qgenesIDs
a11=round(D.n_L.BIC.Post_c[Paper.qgenesIDs,],3);rownames(a11) = names(Paper.qgenesIDs);colnames(a11) = c('A','B','C','D')
a11


exp(log(10))

exp(D.n_L[Paper.qgenesIDs['TPR'],])
(exp(D.n_L.BIC$IC[Paper.qgenesIDs['TPR'],]))
D.n_L.BIC.Post_c.85[Paper.qgenesIDs['TPR'],]
qB<-function(B){
  return(exp(-B/2))
}
P_x_m<-sapply(D.n_L.BIC$IC[Paper.qgenesIDs['TPR'],],FUN = qB)



David.Ardell{
  ALL.MUs
  ALL.SDs
  

  RN=list(); RN.S=list()
  RN[[1]]=rownames(ALL.MUs)               ;RN.S[[1]]='ALL_Genes'
  RN[[2]]=rownames(D.n_L.BIC.Post_a)      ;RN.S[[2]]='group_a'
  RN[[3]]=rownames(D.n_L.BIC.Post_a.85)   ;RN.S[[3]]='group_a.85'
  RN[[4]]=rownames(D.n_L.BIC.Post_b)      ;RN.S[[4]]='group_b'
  RN[[5]]=rownames(D.n_L.BIC.Post_b.85)   ;RN.S[[5]]='group_b.85'
  RN[[6]]=rownames(D.n_L.BIC.Post_c)      ;RN.S[[6]]='group_c'
  RN[[7]]=rownames(D.n_L.BIC.Post_c.85)   ;RN.S[[7]]='group_c.85'
  RN[[8]]=rownames(D.n_L.BIC.Post_d)      ;RN.S[[8]]='group_d'
  RN[[9]]=rownames(D.n_L.BIC.Post_d.85)   ;RN.S[[9]]='group_d.85'
  
  for(i in c(1:9)){
    rn=RN[[i]]
    rn.S=RN.S[[i]]
    print(rn.S)
    pdf(paste0(wd,'/David_',rn.S,'.pdf'),width = 12,8)
      par(mfrow=c(1,2))
      plot(ALL.MUs[rn,1],ALL.SDs[rn,1],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu a' ,ylab='SD a',main=rn.S )
      plot(1,1)
      plot(ALL.MUs[rn,2],ALL.SDs[rn,2],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu b0',ylab='SD b',main=rn.S )
      plot(ALL.MUs[rn,3],ALL.SDs[rn,2],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu b1',ylab='SD b',main=rn.S )
      plot(ALL.MUs[rn,4],ALL.SDs[rn,3],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu c0',ylab='SD c',main=rn.S )
      plot(ALL.MUs[rn,5],ALL.SDs[rn,3],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu c1',ylab='SD c',main=rn.S )
      plot(ALL.MUs[rn,6],ALL.SDs[rn,4],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu d0',ylab='SD d',main=rn.S )
      plot(ALL.MUs[rn,7],ALL.SDs[rn,4],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(0,2.5),xlab='Mu d1',ylab='SD d',main=rn.S )
      par(mfrow=c(1,3))
      plot(ALL.MUs[rn,2],ALL.MUs[rn,3],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(4,15),xlab='Mu b0',ylab='Mu b1',main=rn.S );abline(a=0,b=1)
      plot(ALL.MUs[rn,4],ALL.MUs[rn,5],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(4,15),xlab='Mu c0',ylab='Mu c1',main=rn.S );abline(a=0,b=1)
      plot(ALL.MUs[rn,6],ALL.MUs[rn,7],col=ggplot2::alpha('black',.10),pch=20,xlim=c(4,15),ylim=c(4,15),xlab='Mu d0',ylab='Mu d1',main=rn.S );abline(a=0,b=1)
  
      # par(mfrow=c(3,1))
      plot(abs(ALL.MUs[rn,3] - ALL.MUs[rn,2]),  ALL.SDs[rn,2],col=ggplot2::alpha('black',.10),pch=20,ylim=c(0,2.5),xlab='abs( Mu b1 - Mu b0)',ylab='SD b',main=rn.S )
        abline(lm(abs(ALL.MUs[rn,3] - ALL.MUs[rn,2]) ~ ALL.SDs[rn,2]), col="red") # regression line (y~x) 
        lines(lowess(abs(ALL.MUs[rn,3] - ALL.MUs[rn,2]) ,ALL.SDs[rn,2]), col="blue") # lowess line (x,y)
      
      plot(abs(ALL.MUs[rn,5] - ALL.MUs[rn,4]),  ALL.SDs[rn,3],col=ggplot2::alpha('black',.10),pch=20,ylim=c(0,2.5),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )
        abline(lm(abs(ALL.MUs[rn,5] - ALL.MUs[rn,4]) ~ ALL.SDs[rn,3]), col="red") # regression line (y~x) 
        lines(lowess(abs(ALL.MUs[rn,5] - ALL.MUs[rn,4]) ,ALL.SDs[rn,3]), col="blue") # lowess line (x,y)
      
      plot(abs(ALL.MUs[rn,7] - ALL.MUs[rn,6]),  ALL.SDs[rn,4],col=ggplot2::alpha('black',.10),pch=20,ylim=c(0,2.5),xlab='abs( Mu d1 - Mu d0)',ylab='SD d',main=rn.S )
         abline(lm(abs(ALL.MUs[rn,7] - ALL.MUs[rn,6]) ~ ALL.SDs[rn,4]), col="red") # regression line (y~x) 
         lines(lowess(abs(ALL.MUs[rn,7] - ALL.MUs[rn,6]) ,ALL.SDs[rn,4]), col="blue") # lowess line (x,y)
           
      par(mfrow=c(1,1))
      dd=list()
      dd$'a' = density(ALL.MUs[rn,1])
      dd$'b0' = density(ALL.MUs[rn,2])
      dd$'b1' = density(ALL.MUs[rn,3])
      dd$'c0' = density(ALL.MUs[rn,4])
      dd$'c1' = density(ALL.MUs[rn,5])
      dd$'d0' = density(ALL.MUs[rn,6])
      dd$'d1' = density(ALL.MUs[rn,7])
      plot(dd$a,ylim=c(0,max(sapply(dd, function(x) max(x$y)))  ),col=  brewer.pal(12,'Paired')[12] ,main=rn.S)
      lines(dd$b0,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[1] )
      lines(dd$b1,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[2] )
      lines(dd$c0,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[3] )
      lines(dd$c1,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[4] )
      lines(dd$d0,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[5] )
      lines(dd$d1,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[6] )
      legend('topright',colnames(ALL.MUs),fill = brewer.pal(12,'Paired')[c(12,1:6)] )
      
      plot(cumsum(ALL.MUs[rn,1]),col=  brewer.pal(12,'Paired')[12],type='l' ,main=rn.S)
      lines(cumsum(ALL.MUs[rn,2]),col=  brewer.pal(12,'Paired')[1] )
      lines(cumsum(ALL.MUs[rn,3]),col=  brewer.pal(12,'Paired')[2] )
      lines(cumsum(ALL.MUs[rn,4]),col=  brewer.pal(12,'Paired')[3] )
      lines(cumsum(ALL.MUs[rn,5]),col=  brewer.pal(12,'Paired')[4] )
      lines(cumsum(ALL.MUs[rn,6]),col=  brewer.pal(12,'Paired')[5] )
      lines(cumsum(ALL.MUs[rn,7]),col=  brewer.pal(12,'Paired')[6] )
      legend('topright',colnames(ALL.MUs),fill = brewer.pal(12,'Paired')[c(12,1:6)] )
      
      par(mfrow=c(1,1))
      dd=list()
      dd$'a' = density(ALL.SDs[rn,1])
      dd$'b' = density(ALL.SDs[rn,2])
      dd$'c' = density(ALL.SDs[rn,3])
      dd$'d' = density(ALL.SDs[rn,4])
      plot(dd$a,ylim=c(0,max(sapply(dd, function(x) max(x$y)))  ),col=  brewer.pal(12,'Paired')[12] ,main=rn.S)
      lines(dd$b,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[1] )
      lines(dd$c,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[3] )
      lines(dd$d,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[5] )
      legend('topright',colnames(ALL.SDs),fill = brewer.pal(12,'Paired')[c(12,1,3,5)] )
      
      plot(cumsum(ALL.SDs[rn,1]),col=  brewer.pal(12,'Paired')[12],type='l' ,main=rn.S)
      lines(cumsum(ALL.SDs[rn,2]),col=  brewer.pal(12,'Paired')[1] )
      lines(cumsum(ALL.SDs[rn,3]),col=  brewer.pal(12,'Paired')[3] )
      lines(cumsum(ALL.SDs[rn,4]),col=  brewer.pal(12,'Paired')[5] )
      legend('topright',colnames(ALL.SDs),fill = brewer.pal(12,'Paired')[c(12,1,3,5)] )
    dev.off()

 
  }
  
  
  dd=list()
  dd$'a' = density(ALL.SDs[RN[[2]],1])
  dd$'b' = density(ALL.SDs[RN[[4]],2])
  dd$'c' = density(ALL.SDs[RN[[6]],3])
  dd$'d' = density(ALL.SDs[RN[[8]],4])
  plot(dd$a,ylim=c(0,max(sapply(dd, function(x) max(x$y)))  ),col=  brewer.pal(12,'Paired')[12] ,main=rn.S)
  lines(dd$b,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[1] )
  lines(dd$c,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[3] )
  lines(dd$d,ylim=c(0,max(sapply(dd, function(x) max(x$y)))),col=  brewer.pal(12,'Paired')[5] )
  legend('topright',colnames(ALL.SDs),fill = brewer.pal(12,'Paired')[c(12,1,3,5)] )
 
  list('a'=ALL.SDs[RN[[2]],1],'b'=ALL.SDs[RN[[4]],2])
  
  boxplot2(x=list('a'=ALL.SDs[RN[[2]],1],'b'=ALL.SDs[RN[[4]],2],'c'=ALL.SDs[RN[[6]],3] ,'d'=ALL.SDs[RN[[8]],4]), top=FALSE,ylim=c(-0.05,1.5));abline(h=0,col='gray',)

  boxplot2(x=list('a'=ALL.SDs[RN[[3]],1],'b'=ALL.SDs[RN[[5]],2],'c'=ALL.SDs[RN[[7]],3] ,'d'=ALL.SDs[RN[[9]],4]), top=FALSE,ylim=c(-0.05,1.5));abline(h=0,col='gray',)
  
  for(i in 2:9){
    boxplot(x=list('a'=ALL.SDs[RN[[i]],1],'b'=ALL.SDs[RN[[i]],2],'c'=ALL.SDs[RN[[i]],3] ,'d'=ALL.SDs[RN[[i]],4]));abline(h=0,col='gray')
    print(   table(apply(cbind('a'=ALL.SDs[RN[[i]],1],'b'=ALL.SDs[RN[[i]],2],'c'=ALL.SDs[RN[[i]],3] ,'d'=ALL.SDs[RN[[i]],4]), 1,function(x) which.min(x) )))
    boxplot(x=list(
      'b0'=ALL.MUs[RN[[i]],2]-ALL.MUs[RN[[i]],3] ,
      'c0'=ALL.MUs[RN[[i]],4]-ALL.MUs[RN[[i]],5] ,
      'd0'=ALL.MUs[RN[[i]],6]-ALL.MUs[RN[[i]],7] ),ylim=c(-5,5));abline(h=0,col='gray')
    
    print(  table(apply(cbind(
      'b0'=ALL.MUs[RN[[i]],2]-ALL.MUs[RN[[i]],3] ,
      'c0'=ALL.MUs[RN[[i]],4]-ALL.MUs[RN[[i]],5] ,
      'd0'=ALL.MUs[RN[[i]],6]-ALL.MUs[RN[[i]],7] ),1,function(x) which.max(abs(x)) )))
    }
  
  
  # tmp = melt(ALL.SDs[c(RN[[2]][1:3],RN[[4]][1:3],RN[[6]][1:3],RN[[8]][1:3]),])
  
  #In use date() -- Mon Jul 25 15:32:59 2016
  PAPERPLOTs_2016_07_25{
   ############ 
    ##BOXPLOT_all_SD
    tmp = melt(ALL.SDs[c(RN[[2]],RN[[4]],RN[[6]],RN[[8]]),])
    tmp = cbind(tmp,c(rep('A',length(RN[[2]])), rep('B',length(RN[[4]])) , rep('C',length(RN[[6]])), rep('D',length(RN[[8]]))))

    # ##BOXPLOT_all_SD_0.85
    # tmp = melt(ALL.SDs[c(RN[[3]],RN[[5]],RN[[7]],RN[[9]]),])
    # tmp = cbind(tmp,c(rep('A',length(RN[[3]])), rep('B',length(RN[[5]])) , rep('C',length(RN[[7]])), rep('D',length(RN[[9]]))))
    
    tmp = cbind(tmp,'No')
    colnames(tmp) = c('Gene','group','SD','GROUP','col')
    levels(tmp$col) = c("No","Yes")
    tmp[toupper(tmp$group) == tmp$GROUP,'col' ] = 'Yes'
    to_string <- as_labeller(c(`A` = "Genes assigned to group \n a", `B` = "Genes assigned to group \n b", `C` = 'Genes assigned to group  \n c', `D`='Genes assigned to group  \n d'))
    # to_string <- as_labeller(c(`A` = "signifcant genes of group a", `B` = "signifcant genes of group b", `C` = 'signifcant genes of group c', `D`='signifcant genes of group d'))
    
    cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6])
    g2= ggplot( aes(group, `SD`,fill=col),data = tmp)+ geom_boxplot(notch = F,varwidth = F,outlier.colour = "gray", outlier.shape = 1,coef=1.5)+facet_wrap(~ GROUP,labeller = to_string )
    g3 =g2+scale_fill_manual(values =cols,name="Assigned group") +thememap(14,0.6)+theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    tmp2= tmp[tmp$col == 'Yes',]
    g2= ggplot( aes(group, `SD`),data = tmp2)+ geom_boxplot(notch = F,varwidth = F,outlier.colour = "gray", outlier.shape = 1,coef=1.5)
    g2+scale_fill_manual(values =cols,name="Assigned group") +thememap(14,0.6)+theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))
    
    # In use date() -- Mon Jul 25 15:32:59 2016
    pdf(paste0(wd,"/List/PAPER_BOXPLOT_all_SD_2016_07_25.pdf"),height=6,width=8)
    plot(g3)
    dev.off()
    ############ 
    
    tmp = melt(ALL.MUs[c(RN[[2]],RN[[4]],RN[[6]],RN[[8]]),])
    tmp = cbind(tmp,c(rep('A',length(RN[[2]])), rep('B',length(RN[[4]])) , rep('C',length(RN[[6]])), rep('D',length(RN[[8]]))))
    tmp = cbind(tmp,'No')
    colnames(tmp) = c('Gene','group','Mean of log2 expression','GROUP','col')
    levels(tmp$col) = c("No","Yes")
    tmp[toupper(sub(0,'',tmp$group)) == tmp$GROUP,'col' ] = 'Yes'
    tmp[toupper(sub(1,'',tmp$group)) == tmp$GROUP,'col' ] = 'Yes'
    
    # to_string <- as_labeller(c(`A` = "signifcant genes of group a", `B` = "signifcant genes of group b", `C` = 'signifcant genes of group c', `D`='signifcant genes of group d'))
    # cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6])
    g2= ggplot( aes(group, `Mean of log2 expression`,fill=col),data = tmp)+ geom_boxplot(notch = F,varwidth = F,outlier.colour = "gray", outlier.shape = 1,coef=1.5)+facet_wrap(~ GROUP,labeller = to_string )
    g2+scale_fill_manual(values =cols,name="Assigned group") +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    g3 =g2+scale_fill_manual(values =cols,name="Assigned group") +thememap(14,0.6)+theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    
    # In use date() -- Mon Jul 25 15:32:59 2016
    pdf(paste0(wd,"/List/PAPER_BOXPLOT_all_MEAN_2016_07_25.pdf"),height=6,width=8)
     plot(g3)
    dev.off()
    
    summary(tmp2[tmp2$GROUP=='C','SD'])
    # tmp2= tmp[tmp$col == 'Yes',]
    col2=brewer.pal(8,'Paired')[c(7,2,6,4)] 
    g2= ggplot( aes(group, `SD`,fill=group),data = tmp2)+ geom_boxplot(notch = F,varwidth = F,outlier.colour = "gray", outlier.shape = 1,coef=1.5)
    g4=g2 +thememap(14,0.6)+theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+scale_fill_manual(values =col2 ,name="Assigned group")+theme(legend.position="none")
    
    col2=brewer.pal(8,'Paired')[c(7,1,2,5,6,3,4)] 
    tmp3= tmp[tmp$col == 'Yes',]
    g2= ggplot( aes(group, `Mean of log2 expression`,fill=group),data = tmp3)+ geom_boxplot(notch = F,varwidth = F,outlier.colour = "gray", outlier.shape = 1,coef=1.5)
    g5=g2+thememap(14,0.6)+theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+scale_fill_manual(values =col2 ,name="Assigned group")+theme(legend.position="none")
    
    library("gridExtra")
    pdf(paste0(wd,"/List/PAPER_BOXPLOT_group_SD_MEAN_2016_07_25.pdf"),height=6,width=8)
     grid.arrange(g4 ,g5 ,as.table=TRUE,nrow=2,heights=c(1,1))  
    dev.off()
    #############
      
    tmp=rbind(cbind(ALL.SDs[RN[[4]],'b'],FC.B[RN[[4]],'log2FC'],D.n_L.BIC.Post_b[RN[[4]],2]) ,
              cbind(ALL.SDs[RN[[6]],'b'],FC.C[RN[[6]],'log2FC'],D.n_L.BIC.Post_c[RN[[6]],3]) ,
              cbind(ALL.SDs[RN[[8]],'b'],FC.D[RN[[8]],'log2FC'],D.n_L.BIC.Post_d[RN[[8]],4]) )  
    tmp = as.data.frame(tmp)
    tmp = cbind(tmp,c(rep('B',length(RN[[4]])) , rep('C',length(RN[[6]])), rep('D',length(RN[[8]]))))
    colnames(tmp) = c('SD','log2FC','Posterior','GROUP')
    
    # tmp = tmp[tmp$Posterior >.7,]
    #  
    # tmp=rbind(cbind(ALL.SDs[RN[[5]],'b'],FC.B[RN[[5]],'log2FC'],D.n_L.BIC.Post_b[RN[[5]],2]) ,
    #           cbind(ALL.SDs[RN[[7]],'b'],FC.C[RN[[7]],'log2FC'],D.n_L.BIC.Post_c[RN[[7]],3]) ,
    #           cbind(ALL.SDs[RN[[9]],'b'],FC.D[RN[[9]],'log2FC'],D.n_L.BIC.Post_d[RN[[9]],4]) )  
    # tmp = as.data.frame(tmp)
    # tmp = cbind(tmp,c(rep('B',length(RN[[5]])) , rep('C',length(RN[[7]])), rep('D',length(RN[[9]]))))
    # colnames(tmp) = c('SD','log2FC','Posterior','GROUP')
    
    g2 = ggplot(tmp, aes(x=log2FC, y=SD)) + geom_point(shape=20,alpha = 1/8)+facet_wrap(~ GROUP,nrow = 3,labeller = to_string )# + geom_smooth()
    g3 = g2 + thememap(14,0.7) +theme(legend.background = element_rect(fill="grey90", size=5.5, linetype="dotted"))+theme(legend.position="bottom")
    g3
    # In use date() -- Mon Jul 25 15:32:59 2016
    pdf(paste0(wd,"/List/PAPER_ScatterPlot_all_2016_08_15.pdf"),height=6,width=8)
      plot(g3)
    dev.off()
    
        g2 = ggplot(tmp, aes(x=log2FC, y=SD)) + geom_point(shape=20,aes(colour = Posterior),alpha = 1) +xlim(c(-3,3))   + scale_colour_gradient(low = "white", high = "black",limits=c(0.25, 1)) +facet_wrap( ~GROUP,nrow = 3 ,labeller = to_string )# + geom_smooth()
        g3 = g2 + thememap(14,0.7) +theme(legend.background = element_rect(fill="grey90", size=5.5, linetype="dotted"))+theme(legend.position="bottom")
        g3
        # g3 + geom_point(data=tmp[Paper.qgenesIDs,], aes(x=log2FC, y=SD), colour="red", size=2,shape=1)
        
        g2 = ggplot(tmp, aes(x=log2FC, y= Posterior)) + geom_point(shape=20,aes(colour = SD),alpha = 1)   + scale_colour_gradient(low = "black", high = "white",limits=c(0, 1)) +facet_wrap( ~GROUP,nrow = 3 ,labeller = to_string )# + geom_smooth()
        g3 = g2 + thememap(14,1) +theme(legend.background = element_rect(fill="grey90", size=5.5, linetype="dotted"))+theme(legend.position="bottom")
        g3
        
        g2 = ggplot(tmp, aes(x=SD, y= Posterior)) + geom_point(shape=20,aes(colour = log2FC),alpha = 1)   + scale_colour_gradient2(low='black',mid='white',high='red',limits=c(-1, 1)) +facet_wrap( ~GROUP,nrow = 3 ,labeller = to_string )# + geom_smooth()
        g3 = g2 + thememap(14,0.7) +theme(legend.background = element_rect(fill="grey90", size=5.5, linetype="dotted"))+theme(legend.position="bottom")
        g3
      
    tmpC = tmp[tmp$GROUP == 'C',]
    PostCut=.75
    tmpC = tmpC[tmpC$Posterior >PostCut,] #>.85,]
    to_string <- as_labeller(c(`A` = "Genes assigned to group \n a", `B` = "Genes assigned to group \n b", `C` = 'Genes assigned to group c & \n Posterior > 0.75', `D`='Genes assigned to group  \n d'))
    g2 = ggplot(tmpC, aes(x=log2FC, y=SD))+ geom_vline(xintercept = c(-0.5,0.5),colour='gray') + geom_point(shape=20,alpha = 1) +xlim(c(-2,2)) + ylim(c(0,0.8)) +facet_wrap( ~GROUP,nrow = 3 ,labeller = to_string )# + geom_smooth()
    g3 = g2 + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom") 
    g3 = g3 + geom_point(data=tmp[Paper.qgenesIDs,], aes(x=log2FC, y=SD), colour="red", size=3,shape=21)+xlab('Log2 fold change for group c')
     
    g3
    # In use date() -- Mon Jul 25 15:32:59 2016
    pdf(paste0(wd,"/List/PAPER_ScatterPlot_groupC_2016_08_15.pdf"),height=6,width=8)
    plot(g3)
    dev.off()
    
    tmp[Paper.qgenesIDs,]
    
    pg.M=ALL.MUs[Paper.qgenesIDs,]
    pg.s=ALL.SDs[Paper.qgenesIDs,]
    
    pdf(paste0(wd,"/List/PAPER_NV_density_all_2016_07_25.pdf"),height=6,width=8)
    par(mfrow=c(2,3)  )
    for(i in 1:length(Paper.qgenesIDs)){
      # pdf(paste0(wd,"/List/PAPER_NV_density_",names(Paper.qgenesIDs)[i],"_2016_07_25.pdf"),height=6,width=8)
        cols=brewer.pal(8,'Paired')
        x <- seq(4, 12, length=10000)
        plot(0,xlim=c(4,12),ylim=c(0,4.5),pch=0,col='white',main=names(Paper.qgenesIDs)[i],ylab='density',xlab='log2 expression',cex.axis = 1.3,cex.lab=1.3)
        curve((dnorm(x,mean =  pg.M[i,'a0'], sd = pg.s[i,'a'])), add = TRUE, col = cols[7], lwd = 2)
    
        curve((dnorm(x,mean =  pg.M[i,'b0'], sd = pg.s[i,'b'])), add = TRUE, col = cols[1], lwd = 2)
        curve((dnorm(x,mean =  pg.M[i,'b1'], sd = pg.s[i,'b'])), add = TRUE, col = cols[2], lwd = 2)
        
        curve((dnorm(x,mean =  pg.M[i,'c0'], sd = pg.s[i,'c'])), add = TRUE, col = cols[5], lwd = 2)
        curve((dnorm(x,mean =  pg.M[i,'c1'], sd = pg.s[i,'c'])), add = TRUE, col = cols[6], lwd = 2)
        
        curve((dnorm(x,mean =  pg.M[i,'d0'], sd = pg.s[i,'d'])), add = TRUE, col = cols[3], lwd = 2)
        curve((dnorm(x,mean =  pg.M[i,'d1'], sd = pg.s[i,'d'])), add = TRUE, col = cols[4], lwd = 2)
        legend('topright',colnames(ALL.MUs),fill = brewer.pal(8,'Paired')[c(7,1,2,5,6,3,4)] ,cex = 1.2)
      # dev.off()
    }
    dev.off()
    # tmp =melt(pg.M)
    # gg <- ggplot(tmp , aes(x=value))
    # gg <- gg + stat_function(fun=dnorm,
    #                          color="red",
    #                          args=list(mean=0.8, 
    #                                    sd=0.5))+facet_wrap( ~X1)
    # gg
    # 
    
    
    pg.M=ALL.MUs[Paper.qgenesIDs,]
    pg.s=ALL.SDs[Paper.qgenesIDs,]
    for(i in 1:length(Paper.qgenesIDs)){
      # pdf(paste0(wd,"/List/PAPER_NV_density_",names(Paper.qgenesIDs)[i],"_2016_09_14.pdf"),height=6,width=8)
      pdf(paste0(wd,"/List/PAPER_NV_density_",names(Paper.qgenesIDs)[i],"_2016_09_22.pdf"),height=6,width=8) # optimized for TPR
      logE = D.n$E[Paper.qgenesIDs[i],]
      tmp = data.frame(logE=rep(logE,4),
                       density = rep( 0, length(rep(logE,4))),
                       group    =c(rep('a',6),rep('b',6),rep('c',6),rep('d',6)),
                       indicator= factor(c( rep(0,6) , c(0,1,0,1,0,1), c(0,1,0,0,0,1) ,c(0,1,0,0,0,0) ))
      )
      
      x <- seq(4, 12, length=1000)
      
      dM=round(max(dnorm(x,mean =  pg.M[i,'c1'], sd = pg.s[i,'c']))+0.5)
      
      col2=brewer.pal(8,'Paired')[c(2,6)] 
      col2rgb(col2)
      
      g2 = ggplot(tmp, aes(x=logE,y=density,colour=indicator))+scale_y_continuous(limits = c(0, dM),breaks = seq(0,dM,1))+scale_x_continuous(limits = c(6, 10),breaks = seq(6,10,1)) +
        geom_point(shape=20,size =0)  +facet_wrap(~ group,nrow = 4)  + 
        scale_color_manual(values =col2,name="",labels = list("g = 0","g = 1" ))+  #scale_color_manual(values =col2,name=" Indicator variable",labels = list("g = 0","g = 1" ))+
        labs(title=names(Paper.qgenesIDs)[i],y = "Density", x = "Logarithmic expression levels")
      g3 =g2+ with(tmp[tmp$group=="a",],stat_function(data=tmp[tmp$group=="a" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M[i,'a0'], sd = pg.s[i,'a']), size=1.2,colour="black") ) + #new ,colour="black"
        with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M[i,'b0'], sd = pg.s[i,'b']), size=1.2) ) +
        with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M[i,'b1'], sd = pg.s[i,'b']), size=1.2) ) +
        with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M[i,'c0'], sd = pg.s[i,'c']), size=1.2) ) +
        with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M[i,'c1'], sd = pg.s[i,'c']), size=1.2) ) +
        with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M[i,'d0'], sd = pg.s[i,'d']), size=1.2) ) +
        with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M[i,'d1'], sd = pg.s[i,'d']), size=1.2) )
      g4 = g3 + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom") 
      g4 =  
        g4 + geom_point(data=tmp, aes(x=logE,y=density),shape=20,size =3.5) +
        with( tmp[tmp$group=="a",],geom_point(data=tmp[tmp$group=="a",], aes(x=logE,y=density),shape=20,size =3.5,colour="black"))+ #new
        geom_point(data=tmp, aes(x=logE,y=density),shape=21,size =3,colour="black")
      print(g4)
      dev.off()
      
      ### do" use -> sd is not calc in the same way as in our model
      # g2 = ggplot(tmp, aes(x=logE,y=density,colour=indicator))+ylim(c(0,4))+xlim(c(5,8)) + geom_point(shape=20,size =4)  +facet_wrap(~ group,nrow = 4)  + scale_color_manual(values =col2)
      # g2+ with(tmp[tmp$group=="a",],stat_function(data=tmp[tmp$group=="a" & tmp$indicator ==0,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="a" & tmp$indicator ==0,'logE']), sd = sd(tmp[tmp$group=="a" & tmp$indicator ==0,'logE']) ))  ) +
      #     with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==0,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="b" & tmp$indicator ==0,'logE']), sd = sd(tmp[tmp$group=="b" & tmp$indicator ==0,'logE']) ))  ) +
      #     with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==1,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="b" & tmp$indicator ==1,'logE']), sd = sd(tmp[tmp$group=="b" & tmp$indicator ==1,'logE']) ))  ) +
      #     with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==0,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="c" & tmp$indicator ==0,'logE']), sd = sd(tmp[tmp$group=="c" & tmp$indicator ==0,'logE']) ))  ) +
      #     with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==1,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="c" & tmp$indicator ==1,'logE']), sd = sd(tmp[tmp$group=="c" & tmp$indicator ==1,'logE']) ))  ) +
      #     with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==0,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="d" & tmp$indicator ==0,'logE']), sd = sd(tmp[tmp$group=="d" & tmp$indicator ==0,'logE']) ))  ) +
      #     with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==1,],fun=dnorm,args=list(mean = mean(tmp[tmp$group=="d" & tmp$indicator ==1,'logE']), sd = sd(tmp[tmp$group=="d" & tmp$indicator ==1,'logE']) ))  ) 
    }
    
    Paper.qgenesIDs
    D.n_L.BIC.Post_c[Paper.qgenesIDs,]
    
    
    
    }  
  
  allgemeiner__DENSITY__PLot{
  for(i in IDs){
    pdf(paste0(wd,"/List/DENS/NV_density_",IDs,".pdf"),height=6,width=8)
    
    pg.M=ALL.MUs[IDs,]
    pg.s=ALL.SDs[IDs,]
    
    logE = D.n$E[IDs,]
    tmp = data.frame(logE=rep(logE,4),
                     density = rep( 0, length(rep(logE,4))),
                     group    =c(rep('a',6),rep('b',6),rep('c',6),rep('d',6)),
                     indicator= factor(c( rep(0,6) , c(0,1,0,1,0,1), c(0,1,0,0,0,1) ,c(0,1,0,0,0,0) ))
    )
    
    x <- seq(4, 12, length=1000)
    dM=round(max(dnorm(x,mean =  pg.M['c1'], sd = pg.s['c']))+0.5)
    
    col2=brewer.pal(8,'Paired')[c(2,6)] 
    g2 = ggplot(tmp, aes(x=logE,y=density,colour=indicator))+scale_y_continuous(limits = c(0, dM),breaks = seq(0,dM,1))+scale_x_continuous(limits = c(4.5, 10),breaks = seq(4.5,10,0.5)) + geom_point(shape=20,size =0)  +facet_wrap(~ group,nrow = 4)  + scale_color_manual(values =col2,name=" Indicator variable",labels = list("g = 0","g = 1" ))+labs(title=IDs,y = "Density", x = "Logarithmic expression levels")
    g3 =g2+ with(tmp[tmp$group=="a",],stat_function(data=tmp[tmp$group=="a" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['a0'], sd = pg.s['a']), size=1.2) ) +
      with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['b0'], sd = pg.s['b']), size=1.2) ) +
      with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M['b1'], sd = pg.s['b']), size=1.2) ) +
      with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['c0'], sd = pg.s['c']), size=1.2) ) +
      with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M['c1'], sd = pg.s['c']), size=1.2) ) +
      with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['d0'], sd = pg.s['d']), size=1.2) ) +
      with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M['d1'], sd = pg.s['d']), size=1.2) )
    g4 = g3 + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom") 
    g4 =  g4 + geom_point(data=tmp, aes(x=logE,y=density),shape=20,size =3.5) +
      geom_point(data=tmp, aes(x=logE,y=density),shape=21,size =3,colour="black")
    print(g4)
    dev.off()
    
  }  
  }
  
    
  rnC=RN[[6]]
  rnC85=RN[[7]]
  ### !!! ###
  plot((ALL.MUs[rnC,5] - ALL.MUs[rnC,4]),FC.C[rnC,'log2FC'])
  library(beanplot)
  beanplot(ALL.MUs[rnC,5],log = "")
  library(vioplot)
  vioplot(ALL.MUs[rnC,5],col="gray",drawRect = T)
  
  d <- density(ALL.MUs[rnC,5])
  plot(d, type="n", main="robbery")
  polygon(d, col="lightgray", border="gray")
  rug(ALL.MUs[rnC,5], col="red")
  
  curve(1:10,dnorm(x),lwd=3)

  x <- seq(-10, 10, length=100)
  plot(0,xlim=c(-10,10),ylim=c(0,3))
  curve((dnorm(x,mean =  6.862880 , sd = 0.01570576)), add = TRUE, col = "red", lwd = 2)
  curve((dnorm(x,mean =  6.447395 , sd = 0.01570576)), add = TRUE, col = "blue", lwd = 2)
  
  head(ALL.MUs[rnC85,])
  head(ALL.SDs[rnC85,])
  
  mean =  6.862880 ; sd = 0.01570576
  x <- seq(-10,10,length=1000)*sd + mean
  hx <- dnorm(x,mean,sd)
  
  mean =  6.447395 
  x2 <- seq(-100,100,length=1000)*sd + mean
  hx2 <- dnorm(x2,mean,sd)
  
  plot(x, hx, type="l", xlab="IQ Values", ylab="",
       main="Normal Distribution", axes=T,xlim=c(5,10))
  lines(x2, hx2,col=2)
  
  library(LSD)
  
  heatscatter( (ALL.MUs[rnC85,5] - ALL.MUs[rnC85,4]),  ALL.SDs[rnC85,3])
  
  groups = list("Green" = 1:50,"Red" = 51:150,"Blue" = 151:937)
  colors = c("darkgreen","darkred","darkblue")
  ellipsescatter((ALL.MUs[rnC85,5] - ALL.MUs[rnC85,4]),  ALL.SDs[rnC85,3],groups,colors,location = "topleft",)
  
  alpha = sapply( D.n_L.BIC.Post_c[rnC,3],function(x) round((x-.01)*100))
  heatscatter( (ALL.MUs[rnC,5] - ALL.MUs[rnC,4]),  ALL.SDs[rnC,3] ,daltonize =T,alpha=alpha)
  heatscatter( (ALL.MUs[rnC,5] - ALL.MUs[rnC,4]),  ALL.SDs[rnC,3] ,cor=FALSE,add.contour=TRUE,color.contour="red",greyscale=TRUE)

  heatscatterpoints((ALL.MUs[rnC,5] - ALL.MUs[rnC,4]),  ALL.SDs[rnC,3],greyscale=TRUE)
  
  col=sapply( D.n_L.BIC.Post_c[rnC85,3],function(x) ggplot2::alpha('black',x-0.5))
  plot((ALL.MUs[rnC85,5] - ALL.MUs[rnC85,4]),  ALL.SDs[rnC85,3],col=col,pch=20,ylim=c(0,.35),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )

  col=sapply( D.n_L.BIC.Post_c[rnC,3],function(x) ggplot2::alpha('black',x-0.3))
  plot((ALL.MUs[rnC,5] - ALL.MUs[rnC,4]),  ALL.SDs[rnC,3],col=col,pch=20,ylim=c(0,2.5),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )
  
  abline(lm((ALL.MUs[rnC,5] - ALL.MUs[rnC,4]) ~ ALL.SDs[rnC,3]), col="red") # regression line (y~x) 
  lines(lowess((ALL.MUs[rnC,5] - ALL.MUs[rnC,4]) ,ALL.SDs[rnC,3]), col="blue") # lowess line (x,y)
  
  plot((ALL.MUs[rnC,5] - ALL.MUs[rnC,4]),  D.n_L.BIC.Post_c[rnC,3],col=ggplot2::alpha('black',.10),pch=20,ylim=c(0,1),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )
  points((ALL.MUs[rnC85,5] - ALL.MUs[rnC85,4]),  D.n_L.BIC.Post_c[rnC85,3],col=ggplot2::alpha('red',.10),pch=20,ylim=c(0,1),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )
  
  
  plot(ALL.SDs[rnC,3],  D.n_L.BIC.Post_c[rnC,3],col=ggplot2::alpha('black',.10),pch=20,ylim=c(0,1),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )
  abline(h=c(0.25))
  points(ALL.SDs[rnC85,3],  D.n_L.BIC.Post_c[rnC85,3],col=ggplot2::alpha('red',.10),pch=20,ylim=c(0,1),xlab='abs( Mu c1 - Mu c0)',ylab='SD c',main=rn.S )
  
  
  
  3d{
    
        plot(ALL.SDs[rn,3],D.n_L.BIC.Post_c[rn,3],col=ggplot2::alpha('black',.10),pch=20)
        plot(ALL.MUs[rn,5],D.n_L.BIC.Post_c[rn,3],col=ggplot2::alpha('black',.10),pch=20)
        plot(D.n_L.BIC.Post_c[rn,3],ALL.SDs[rn,3],col=ggplot2::alpha('black',.10),pch=20)
        
        sca2
        xlim=c(4,15),ylim=c(0,2.5),xlab='Mu c1',ylab='SD c',main=rn.S )
    library(scatterplot3d) 
    library(rgl)
    
    s3d<-scatterplot3d(ALL.MUs[rn,4],ALL.SDs[rn,3],ALL.MUs[rn,5], pch=16, highlight.3d=T,type="p", main="3D Scatterplot")
    fit <- lm(ALL.MUs[rn,4] ~ ALL.MUs[rn,5] +ALL.SDs[rn,3]) 
    s3d$plane3d(fit)
    # plot3d(ALL.MUs[rn,4],ALL.SDs[rn,3],ALL.MUs[rn,5])
    
    s3d<-scatterplot3d(ALL.MUs[rn,4],D.n_L.BIC.Post_c[rn,3],ALL.MUs[rn,5], pch=16, highlight.3d=T,type="p", main="3D Scatterplot")
    fit <- lm(ALL.MUs[rn,4] ~ ALL.MUs[rn,5] +D.n_L.BIC.Post_c[rn,3]) 
    s3d$plane3d(fit)
    
    s3d<-scatterplot3d(ALL.MUs[rn,2],ALL.SDs[rn,2],ALL.MUs[rn,3], pch=16, highlight.3d=T,type="p", main="3D Scatterplot")
    fit <- lm(ALL.MUs[rn,2] ~ ALL.MUs[rn,3] +ALL.SDs[rn,2]) 
    s3d$plane3d(fit)
  }
  
}


xxx{
X<-read.Data()
X.n<-normalize.and.fiter.Data(X)
X.n_L<-get.ML.logLs.var.NEW(X.n$E)
X.n_L.BIC <- make.IC(X.n_L,c(2,3),6) # BIC with 6 = number of values   
X.n_L.BIC.Post<-get.posterior(X.n_L.BIC)
mX<-apply((X.n_L.BIC.Post),MARGIN=1,FUN=maxIndex)
print(table(mX))

t<-X.n_L.BIC.Post[mX==1,] ;X.n_L.BIC.Post_a<-t[order(t[,1],decreasing=T),] 
t<-X.n_L.BIC.Post[mX==2,] ;X.n_L.BIC.Post_b<-t[order(t[,2],decreasing=T),]
t<-X.n_L.BIC.Post[mX==3,] ;X.n_L.BIC.Post_c<-t[order(t[,3],decreasing=T),]
t<-X.n_L.BIC.Post[mX==4,] ;X.n_L.BIC.Post_d<-t[order(t[,4],decreasing=T),]
X.n_L.BIC.Post_a.85<-X.n_L.BIC.Post_a[X.n_L.BIC.Post_a[,1]>0.85,]
X.n_L.BIC.Post_b.85<-X.n_L.BIC.Post_b[X.n_L.BIC.Post_b[,2]>0.85,]
X.n_L.BIC.Post_c.85<-X.n_L.BIC.Post_c[X.n_L.BIC.Post_c[,3]>0.85,]
X.n_L.BIC.Post_d.85<-X.n_L.BIC.Post_d[X.n_L.BIC.Post_d[,4]>0.85,]
dim(X.n_L.BIC.Post_c.85)





intersect(EGFR_genes, IDs[IDs[,1] %in% f.input2(rownames(D.n_L.BIC.Post_c.85),rownames(X.n_L.BIC.Post_c.85))$inter,2])
f.input3(rownames(D.n_L.BIC.Post_c),rownames(X.n_L.BIC.Post_c),Paper.qgenesIDs)
f.input5(rownames(X.n_L.BIC.Post_a),rownames(X.n_L.BIC.Post_b),rownames(X.n_L.BIC.Post_c),rownames(X.n_L.BIC.Post_d),Paper.qgenesIDs)
f.input5(rownames(X.n_L.BIC.Post_a.85),rownames(X.n_L.BIC.Post_b.85),rownames(X.n_L.BIC.Post_c.85),rownames(X.n_L.BIC.Post_d.85),rownames(D.n_L.BIC.Post_c.85))


intersect(IDs[IDs[,2]=="STAT3",1],rownames(D.n_L.BIC.Post_a))
D$E[IDs[IDs[,2]=="STAT3",1],]


D[c("ILMN_1696521","ILMN_1755535","ILMN_1728858","ILMN_1798975"),]
D.n[c("ILMN_1696521","ILMN_1755535","ILMN_1728858","ILMN_1798975"),]
}

qPCR{
  ######## 
  RawQPCR <- read.csv(paste(sep="",wd,"/qPCR/MyData.csv"))
  # setwd('/Users/weinhol/GitHub/BGSC')
  # devtools::use_data(RawQPCR)
  
  
  qgenes<-c("TPR","CKAP2L","KIF5C","ROCK1","BPNT1","GALNS","GLIPR2","KEL","ALDH4A1","CDCP1","CLCA2")
  #qgenes<-c("CKAP2L", "ROCK1", "TPR","ALDH4A1", "CLCA2","GALNS")
  
  IDs<-cbind(rownames(D$E),D$genes);colnames(IDs)<-c("ILMN","SYMBOL")
  #IDs[IDs["SYMBOL"]=="TPR",]
  
  STAT3 = IDs[IDs[,2]=='STAT3',1:3]
  ERBB4 = IDs[ grepl('ERBB4',IDs[,2]),][1,]
  
  qgenesIDs<-list()
  for(i in 1:length(qgenes)){
    qgenesIDs[[i]]<-as.character(IDs[which(qgenes[i]==IDs[,"SYMBOL"]),1])
  }
  names(qgenesIDs)<-qgenes
  #intersect(c(unlist(qgenesIDs),"ILMN_1738468","ILMN_1712718","ILMN_1737949","ILMN_1737949"),EGF.predict.IDs)
  

  
  qGenes<-c()
  qGenes[[1]] <-    cbind(RawQPCR[1:6,c(3,4)]  ,RawQPCR[13:18,c(3,4)]  ,RawQPCR[25:30,c(3,4)])
  qGenes[[2]]<-     cbind(RawQPCR[1:6,c(5,6)]  ,RawQPCR[13:18,c(5,6)]  ,RawQPCR[25:30,c(5,6)])
  qGenes[[3]]<-     cbind(RawQPCR[1:6,c(7,8)]  ,RawQPCR[13:18,c(7,8)]  ,RawQPCR[25:30,c(7,8)])
  qGenes[[4]]<-     cbind(RawQPCR[1:6,c(9,10)] ,RawQPCR[13:18,c(9,10)] ,RawQPCR[25:30,c(9,10)])
  qGenes[[5]] <-    cbind(RawQPCR[1:6,c(11,12)],RawQPCR[13:18,c(11,12)],RawQPCR[25:30,c(11,12)])
  qGenes[[6]] <-    cbind(RawQPCR[1:6,c(13,14)],RawQPCR[13:18,c(13,14)],RawQPCR[25:30,c(13,14)])
  qGenes[[7]] <-    cbind(RawQPCR[1:6,c(15,16)],RawQPCR[13:18,c(15,16)],RawQPCR[25:30,c(15,16)])
  qGenes[[8]] <-    cbind(RawQPCR[1:6,c(17,18)],RawQPCR[13:18,c(17,18)],RawQPCR[25:30,c(17,18)])
  qGenes[[9]] <-    cbind(RawQPCR[1:6,c(19,20)],RawQPCR[13:18,c(19,20)],RawQPCR[25:30,c(19,20)])
  qGenes[[10]] <-   cbind(RawQPCR[1:6,c(21,22)],RawQPCR[13:18,c(21,22)],RawQPCR[25:30,c(21,22)])
  qGenes[[11]] <-   cbind(RawQPCR[1:6,c(23,24)],RawQPCR[13:18,c(23,24)],RawQPCR[25:30,c(23,24)])
  names(qGenes)<-qgenes
  
  ###
  # qGenes$'GLIPR2' wahrscheinlich mit Fehler bei GLIPR2relative.Werte 
  
  qGenes.M<-c()
  for(i in 1:length(qgenes)){
    qGenes.M[[i]]<-cbind(rowMeans((qGenes[[i]][c(1,3,5)])),rowMeans(log2(qGenes[[i]][c(2,4,6)])) )#,  log2(rowMeans(qGenes[[i]][c(2,4,6)])) )
    colnames(qGenes.M[[i]])<-c("MeanCT","log2geoMeanREL" ) #,"log2MeanREL")
  }
  names(qGenes.M)<-qgenes
  
  #qgenes<-c()
  c0=c(1,3,4,5);c1=c(2,6)
  qGenes.FC<-list()
  for(i in 1:length(qgenes)){
    qGenes.FC[[i]]<-colMeans(qGenes.M[[i]][c1,])- colMeans(qGenes.M[[i]][c0,])
  }
  names(qGenes.FC)<-qgenes
  qGenes.FC <- lapply(qGenes.FC,function(x) x* c(-1,1)) #,1))
  
  barplot(t(do.call(rbind,qGenes.FC)),beside=T,las=2)   # qGenes$'GLIPR2' wahrscheinlich mit Fehler bei GLIPR2relative.Werte 
  
  ###
  Paper.qgenes<-c("CKAP2L", "ROCK1", "TPR","ALDH4A1", "CLCA2","GALNS")
  Paper.qgenesIDs<-c(qgenesIDs[[2]],qgenesIDs[[4]][2],qgenesIDs[[1]][2],qgenesIDs[[9]][3],qgenesIDs[[11]][1],qgenesIDs[[6]][1])
  names(Paper.qgenesIDs) = Paper.qgenes
  Paper.qgenes[!Paper.qgenesIDs %in% rownames(FC.C)]
  Paper.qgenes[Paper.qgenesIDs %in% rownames(FC.C)]
  Paper.qgenes.FC<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5)])
  Paper.qgenes.FC2<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,3]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,3]))[,c(1,3,4,5)])
  
  Paper.qgenes.FC.bb<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,4,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,5)])
  Paper.qgenes.FC.cc<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5)])
  Paper.qgenes.FC.dd<-t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2)] - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5,6)])
  
  paper.1<-as.matrix(FC.C[rownames(FC.C) %in% Paper.qgenesIDs,])
  paper.1<-paper.1[Paper.qgenesIDs,]
  rownames(paper.1)<-names(Paper.qgenesIDs) #Paper.qgenes
  paper.1<-cbind(paper.1,Paper.qgenes.FC)
  paper.1<-cbind(paper.1,Paper.qgenes.FC2)

  pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new1.pdf"),height=7,width=8)
  barplot(log2(t(paper.1[,1:6])),ylim=c(0,10),beside=T,col=brewer.pal(11,"RdYlBu")[c(1:4,10,9)],legend.text=paste(sep=".","x",c(1,3,4,5,2,6)),ylab="logarithmic expression levels")#,ylab="avg microarray signal")
  dev.off()
  
  pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new2.pdf"),height=7,width=7)
  barplot(t(paper.1[,c(1,5,2,3,4,6)]),beside=T,col=brewer.pal(11,"RdYlBu")[c(1,10,2:4,9)],legend.text=paste(sep=".","x",c(1:6)),ylab="avg microarray signal")
  dev.off()
  
  barplot(t(paper.1[,c('log2FC','Paper.qgenes.FC2')]),beside = T,las=2)
  cor(paper.1[,c('log2FC','Paper.qgenes.FC2')])
  
  require("ggplot2");require("reshape");
  p2<-log2(paper.1[,c(1:6)])
  colnames(p2)<-c("x1","x3","x4","x5","x2","x6")
  df<-melt(p2)
  df<-cbind(df,c(rep("c0",24),rep("c1",12)))
  colnames(df)[4]<-"C"
  
  rns <- levels(df$X1)[c(2,5,6,1,3,4)] ;df$X1 <- factor(df$X1, levels = rns)
  rns <- levels(df$X2)[c(1,3,4,5,2,6)] ;df$X2 <- factor(df$X2, levels = rns)
  labl <- list(expression(c[0]), expression(c[1])) 
  
  thememap <- function (base_size = 12,legend_key_size=0.4, base_family = "") {
    theme_gray(base_size = base_size, base_family = base_family) %+replace% 
      theme(title = element_text(face="bold", colour=1,angle=0  ,vjust=1.0, size=base_size),
            axis.title.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
            axis.text.x  = element_text(face="bold", colour=1,angle=0, vjust=0.5, size=base_size),
            strip.text.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
            axis.title.y = element_text(face="bold", colour=1,angle=90 ,vjust=0.5,hjust=.5, size=base_size),
            axis.text.y  = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
            panel.background = element_rect(fill="white"),
            #panel.grid.minor.y = element_line(size=3),
            panel.grid.major = element_line(colour = "white"),
            legend.key.size = unit(legend_key_size, "cm"),
            legend.text = element_text(face="bold" ,colour=1,angle=0  ,vjust=0.0, size=base_size),
            legend.title = element_text(face="bold",colour=1,angle=0  ,vjust=-0.8, size=base_size)           
      )
  }

  cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6])
  cols<-c(brewer.pal(11,"RdBu")[9],brewer.pal(11,"RdBu")[2])  
  
  pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new3_limma.pdf"),height=7,width=10)
  g1<-ggplot(df,aes(x=factor(X2),y=value,fill=C)) + geom_bar(position="dodge",stat="identity",width=.5)+facet_wrap(~ X1)+  #,title = paste(sep=" ","Reads after quality control"))+
    labs(x = "", y = "Logarithmic expression levels")+scale_fill_manual(values = cols,name="Group",labels = labl) + ylim(c(0,10))+ thememap(14,0.6) +
    scale_x_discrete(labels=c(expression(x[1]),expression(x[3]),expression(x[4]),expression(x[5]),expression(x[2]),expression(x[6])))+
    theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
  
  g1+ theme(legend.justification=c(.9,-1.7), legend.position=c(1,0)) #theme(legend.position="bottom")
  dev.off()
  
  # In use date() -- Fri Jul 22 14:28:45 2016
  pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new4_limma.pdf"),height=7,width=10)
  g1<-ggplot(df,aes(x=factor(X2),y=value,fill=C)) + geom_bar(position="dodge",stat="identity",width=.5)+facet_wrap(~ X1)+  #,title = paste(sep=" ","Reads after quality control"))+
    labs(x = "", y = "Logarithmic expression levels")+scale_fill_manual(values = cols,name=" Indicator variable",labels = list("g = 0","g = 1" )) + ylim(c(0,10))+ thememap(14,0.6) +
    scale_x_discrete(labels=c(expression(x[1]),expression(x[3]),expression(x[4]),expression(x[5]),expression(x[2]),expression(x[6])))+
    theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
  
  g1 #+ theme(legend.justification=c(.9,-1.7), legend.position=c(1,0)) #theme(legend.position="bottom")
  dev.off()
  ###
  old{
  qGenes.FC.comp<-c()
  qGenes.FC.comp<-c(qGenes.FC[[1]][c(2,3)],FC.C[c(qgenesIDs[[1]])[2],"log2FC"] )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[2]][c(2,3)],FC.C[c(qgenesIDs[[2]])   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[4]][c(2,3)],FC.C[c(qgenesIDs[[4]])[2],"log2FC"] ) )
  #qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[3]][c(2,3)],log2(FC[c(qgenesIDs[[3]])[1] ])) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[5]][c(2,3)],FC.C[c(qgenesIDs[[5]])   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[6]][c(2,3)],FC.C[c(qgenesIDs[[6]])   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[7]][c(2,3)],FC.C[c(qgenesIDs[[7]])   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[8]][c(2,3)],FC.C[c(qgenesIDs[[8]])   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[9]][c(2,3)],FC.C[intersect(rownames(FC.C),qgenesIDs[[9]])[1]   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[10]][c(2,3)],FC.C[intersect(rownames(FC.C),qgenesIDs[[10]])[1]   ,"log2FC"] ) )
  qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[11]][c(2,3)],FC.C[intersect(rownames(FC.C),qgenesIDs[[11]])[1]   ,"log2FC"] ) )
  rownames(qGenes.FC.comp)<-qgenes[c(-3)]
  qGenes.FC.comp2<-qGenes.FC.comp[-6,]
  
  pdf(paste(sep="",out,"qPCR_plot.pdf"),10,6)
  barplot(t(qGenes.FC.comp[,c(1,3)]),ylim=c(-1.5,1.5),beside=T,las=2,legend.text=c("qPCR","Array"),main="log2Foldchange for group C ")
  dev.off()
  cor(qGenes.FC.comp2[,c(1,3)])
  plot(qGenes.FC.comp2[,c(1,3)],ylim=c(-2,2),xlim=c(-2,2))
  pdf(paste(sep="",out,"qPCR_plot2.pdf"),10,6)
  barplot(t(qGenes.FC.comp2[,c(1,3)]),ylim=c(-1.5,1.5),beside=T,las=2,legend.text=c("qPCR","Array"),main="log2Foldchange for group C ")
  dev.off()
  pdf(paste(sep="",out,"qPCR_plot2b.pdf"),6,6)
  plot(qGenes.FC.comp2[,c(3,1)],ylim=c(-1.5,1.5),xlim=c(-1.5,1.5),pch=19,las=1,ylab="log2 fold change RT-qPCR",xlab="log2 fold change microarray")
  #text(-1.6,2,paste("r =",round(cor(qGenes.FC.comp2[,c(1,3)]),2) ))
  dev.off()
  
  
  
  
    grC.ID<-IDs[rownames(FC.C),c(1,2)]
    #grC.ID.qPCR<-c()
    #for(g in qgenes[c(-3)]){
    #  grC.ID.qPCR<-rbind(grC.ID.qPCR,grC.ID[which(grC.ID[,2]==g),])
    #}
    #grC.ID.qPCR<-grC.ID.qPCR[-c(6,9,11),]
    grC.ID.qPCR<-grC.ID[Paper.qgenesIDs,]
    
    grC.A<-D.n$E[as.character(rownames(grC.ID.qPCR)),]
    
    grC.Q<-c()
    for(g in as.character(grC.ID.qPCR[,2])){
      grC.Q<-rbind(grC.Q,qGenes.M[[g]][,2])
    }
    rownames(grC.Q)<-rownames(grC.A)
    
    grC.Q.log2<-grC.Q #####log2(exp(grC.Q))
    grC.A.log2<-grC.A #log2(exp(grC.A))
    
    grC.Q.FC<-rowMeans(grC.Q.log2[,c1])- rowMeans(grC.Q.log2[,c0])
    grC.A.FC<-rowMeans(grC.A.log2[,c1])- rowMeans(grC.A.log2[,c0])
    
    
    paper.1
    
    FC.B3[as.character(rownames(grC.ID.qPCR)),]
    
    pdf(paste(sep="",out,"qPCR_plot2.pdf"),10,6)
    barplot(t(cbind(grC.Q.FC,grC.A.FC)),ylim=c(-1.5,1.5),beside=T,las=2,names.arg=grC.ID.qPCR[,2],legend.text=c("qPCR","Array"),main="log2Foldchange for group C ")
    #barplot(t(qGenes.FC.comp2[,c(1,3)]),ylim=c(-1.5,1.5),beside=T,las=2,legend.text=c("qPCR","Array"),main="log2Foldchange for group C ")
    dev.off()
    cbind(grC.ID.qPCR,grC.Q.FC,grC.A.FC)
    
    
    for(i in 1:length(rownames(grC.A))){
      #plot(grC.Q[i,c(1,3,4,5)],grC.A[i,c(1,3,4,5)],ylim=c(2,8),xlim=c(5,8),pch=17,ylab="Array",xlab="qPCR")
      plot(grC.A.log2[i,],grC.Q.log2[i,],pch=16,col=0,xlab="Array",ylab="qPCR",main=rownames(grC.Q.log2)[i])
      points(grC.A.log2[i,c(1,3,4,5)],grC.Q.log2[i,c(1,3,4,5)],pch=17)
      points(grC.A.log2[i,c(2,6)],grC.Q.log2[i,c(2,6)],col=2,pch=19)
      legend("topleft",c("+","-"),lty=c(0,0),col=c(2,1),pch=c(19,17))
      abline(lm(grC.Q.log2[i,]~grC.A.log2[i,]))
      print(lm(grC.Q.log2[i,]~grC.A.log2[i,])[[1]][2]) #slope #http://msenux.redwoods.edu/math/R/regression.php
      print(cor(grC.Q.log2[i,],grC.A.log2[i,]))
    }
    
    pdf(paste(sep="",out,"CorrPlot.pdf"),10,10)
    for(i in 1:length(rownames(grC.A))){
      plot(grC.A.log2[i,],grC.Q.log2[i,],pch=c(5,2,7,9,12,6),col=brewer.pal(9,"Set1")[c(2,1,3,4,5,8)],xlab="Array",ylab="qPCR",main=grC.ID.qPCR[i,2])
      legend("topleft",as.character(c(1:6)),lty=c(0,0),col=brewer.pal(9,"Set1")[c(2,1,3,4,5,8)],pch=c(5,2,7,9,12,6))
      abline(lm(grC.Q.log2[i,]~grC.A.log2[i,]))
    }
    dev.off()
  
  ### PAPER PLOT
  take<-c(1:3,5,7,9)
  grC.ID.qPCR[take,]
  grC.A.log2[take,]
  grC.Q.log2[take,]
  grC.A.FC[take]
  grC.Q.FC[take]
  
  
  pdf(paste(sep="",out,"PAPER_qPCR_Array_expression.pdf"),height=5,width=5)
  red<-brewer.pal(9,"Set1")[c(1)]
  corP<-c();corS<-c();
  for(i in take){
    plot(grC.A.log2[i,],grC.Q.log2[i,],pch=c(17,19,17,17,17,19),col=c(1,red,1,1,1,red),xlab="Microarray log2 expression",ylab="RT-qPCR log2 expression",main=grC.ID.qPCR[i,2])
    legend("topleft",as.character(c(0:1)),lty=c(0,0),col=c(1,red),pch=c(17,19))
    abline(lm(grC.Q.log2[i,]~grC.A.log2[i,]))
    print(lm(grC.Q.log2[i,]~grC.A.log2[i,])[[1]][2]) #slope #http://msenux.redwoods.edu/math/R/regression.php
    
    corP<-c(corP,cor(cbind(grC.A.log2[i,],grC.Q.log2[i,]))[2])
    corS<-c(corS,cor(cbind(grC.A.log2[i,],grC.Q.log2[i,]),method="spearman")[2])
  }
  mean(corP)
  dev.off()
  pdf(paste(sep="",out,"PAPER_qPCR_FC.pdf"),height=5,width=10)
  
  grC.FC<-cbind(grC.A.FC[take],grC.Q.FC[take])
  grC.FC<-grC.FC[order(grC.FC[,1],decreasing=T),]
  rownames(grC.FC)<-c("CKAP2L","ROCK1","TPR","ALDH4A1","CLCA2","GALNS")
  colnames(grC.FC)<-c("Microarray","RT-qPCR")
  #barplot(t(grC.FC[c(2,3,1,5,6,4),]),ylim=c(-1.5,1.5),beside=T,las=2,names.arg=grC.ID.qPCR[take,2],legend.text=c("Microarray","RT-qPCR"),ylab="log2FoldChange for group c", main=" ")
  
  barplot(t(grC.FC[c(2,3,1,5,6,4),]),ylim=c(-1.5,1.5),cex.axis=1.3,cex.names=1.3,beside=T,las=2,names.arg=c("CKAP2L","ROCK1","TPR","ALDH4A1","CLCA2","GALNS"),legend.text=c("Microarray","RT-qPCR"),ylab="log2FoldChange for group c", main=" ")
  
  #barplot(t(qGenes.FC.comp2[,c(1,3)]),ylim=c(-1.5,1.5),beside=T,las=2,legend.text=c("qPCR","Array"),main="log2Foldchange for group C ")
  dev.off()
  cor.test(grC.FC[,1],grC.FC[,2],alternative="t")
  }
  
  #QQ<-grC.Q.FC;names(QQ)<-c("CKAP2L","ROCK1","TPR","ALDH4A1","CLCA2","GALNS")
  #grC.FC<-cbind(paper.1[names(QQ),"log2FC"],QQ);colnames(grC.FC)<-c("Microarray","RT-qPCR")
  #grC.FC
  #cbind(grC.ID.qPCR,grC.Q.FC,grC.A.FC)
  #paper.1
  
  grC.FC<-paper.1[,c("log2FC","Paper.qgenes.FC")];
  rownames(grC.FC)<-c("CKAP2L","ROCK1","TPR","ALDH4A1","CLCA2","GALNS")
  colnames(grC.FC)<-c("Microarray","RT-qPCR")
  
  df2<-melt(grC.FC)
  rns <- levels(df2$X1)[c(2,5,6,1,3,4)] ;df2$X1 <- factor(df2$X1, levels = rns)
  cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6])
  #cols<-c(brewer.pal(11,"RdBu")[9],brewer.pal(11,"RdBu")[2])  
  
  # In use date() -- Fri Jul 22 14:28:45 2016
  pdf(paste0(wd,"/List/PAPER_qPCR_FC_new_Limma_2016_07_22.pdf"),height=6,width=8)
  g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
    thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
  g1  
  # g1+ theme(legend.justification=c(.9,-1.7), legend.position=c(1,0)) #theme(legend.position="bottom")
  dev.off()
  cor(grC.FC)
  plot(grC.FC[,1],grC.FC[,2],ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))
  abline(lm(grC.FC[,1]~grC.FC[,2]))
  cor(grC.FC,method="spearman")

  #### with paper PLOT!!!!
  errorbar{
    errorbar2017-01-01{
      # save.image("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2017-01-01.RData")
      load("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2017-01-01.RData")
      
      MAIL_IVO_2016-12-31{
        # computes 3 different estimators of the variance of the mean log fold change
        vars = function (K, M, N, xkm, xk, xm, x, ykn, yk, yn, y, dk, d) {
          v = rep (NA, 3);
          vk = rep (NA, K);
          vm = rep (NA, M);
          vn = rep (NA, N);
          
          
          # best
          # good for large K, good for large M and N
          v[1] = ((K*M-1) * var (as.vector (xkm)) + (K*N-1) * var (as.vector (ykn))) / (K*M-1 + K*N-1);
          # #note mit k=1                 M-werte                         N-Werte
          # v[1] = ((M-1) * var (as.vector (xkm)) + (N-1) * var (as.vector (ykn))) / (M-1 + N-1);
          
          # good
          # bad for large K, good for large M and N
          for (k in 1:K) {
            vk[k] = ((M-1) * var (xkm[k,]) + (N-1) * var (ykn[k,])) / (M-1 + N-1);
          }
          v[2] = mean (vk);
          
          # good
          # good for large K, bad for large M and N
          for (m in 1:M) {
            vm[m] = var (xkm[,m]);
          }
          for (n in 1:N) {
            vn[n] = var (ykn[,n]);
          }
          v[3] = (sum (vm) + sum (vn)) / (M + N);
          
          # # worse
          # v[4] = (var (xk) + var (yk)) / (1/M + 1/N);
          # v[5] = (M * var (xk) + N * var (yk)) / 2;
          # v[6] = ((M-1) * var (xm) + (N-1) * var (yn)) / (M-1 + N-1) * K;
          # 
          # ############### v[8] ist wie v[6] nur ohne das K --- ohne den korrigierten Vorfactor 
          # v[8] = ((M-1) * var (xm) + (N-1) * var (yn)) / (M-1 + N-1 )
          # #### NEU
          # ############
          # 
          # ###############
          # # worst
          # v[7] = var (dk) / (1/M + 1/N);
          
          return (v * (1/M + 1/N) / K);
        }
        
        main_1 = function () {
          var_min = 0.1;
          var_max = 5.0;
          var_delta = 0.1;
          sig_vec = sqrt (seq (var_min, var_max, var_delta));
          L = length (sig_vec);
          A = 3; # A different estimators of sigma^2 in function vars
          
          dmean = rep (NA, L);
          dmean_theo = rep (NA, L);
          dvar = rep (NA, L);
          dvar_theo = rep (NA, L);
          
          dvar_mean = matrix (NA, A, L);
          dvar_var = matrix (NA, A, L);
          
          Es = 1e4; # ensemble size
          K = 3; # number of replicates
          M = 4; # number of x values
          N = 2; # number of y values
          mux = 3.0; # \mu_X
          muy = 5.0; # \mu_Y
          
          for (l in 1:L) {
            sig = sig_vec[l];
            d = rep (NA, Es); # vector of Es log fold changes
            dvaraes = matrix (NA, A, Es);
            for (es in 1:Es) {
              
              # generate data
              xkm = matrix (rnorm (K * M, mux, sig), K, M);
              ykn = matrix (rnorm (K * N, muy, sig), K, N);
              
              # compute means
              xk = rowMeans (xkm);
              xm = colMeans (xkm)
              x = mean (xk);
              yk = rowMeans (ykn);
              yn = colMeans (ykn);
              y = mean (yk);
              
              # compute K log fold changes
              dk = xk - yk;
              
              # compute mean log fold change
              d[es] = mean (dk);
              
              # sanity check
              d1 = x - y;
              if (abs (d[es] - d1) > 1e-14) {
                print (d[es] - d1);
              }
              
              # compute A variance estimators for d by function vars
              dvaraes[,es] = vars (K, M, N, xkm, xk, xm, x, ykn, yk, yn, y, dk, d);
            }
            
            # compute empirical and theoretical means and variances 
            dmean[l] = mean (d);
            dvar[l] = var (d);
            dmean_theo[l] = mux - muy;
            dvar_theo[l] = (1 / M + 1 / N) / K * sig^2;
            
            # compute means and variances of A variance estimators for d
            for (a in 1:A) {
              dvar_mean[a,l] = mean (dvaraes[a,]);
              dvar_var[a,l] = var (dvaraes[a,]);
            }		 
          }
          
          # plot empirical and theoretical means and variances
          plot (dmean_theo, dmean, xlab = "Theoretical mean of log fold changes", ylab = "Mean of log fold changes", type = "b");
          lines (dmean_theo, dmean_theo, col = 2);
          # quartz ();
          plot (dvar_theo, dvar, xlab = "Theoretical variance of log fold changes", ylab = "Variance of log fold changes", type = "b");
          lines (dvar_theo, dvar_theo, col = 2);
          
          # plot means of A variance estimators for d
          # quartz ();
          plot (dvar_theo, dvar_theo, xlab = "Theoretical variance of log fold changes", ylab = "Mean variance of log fold changes", type = "l", col = 2);
          for (a in 1:A) {
            lines (dvar_theo, dvar_mean[a,], col = a+2);
          }
          
          # plot variances of A variance estimators for d
          # quartz ();
          plot (dvar_theo, dvar_var[1,], xlab = "Theoretical variance of log fold changes", ylab = "Variance of variances of log fold changes", type = "l", col = 3);
          for (a in 2:A) {
            lines (dvar_theo, dvar_var[a,], col = a+2);
          }
        }
        
        main_1 ();
        
        
        MAIL{
          # 31.12.2016
          # Hi Claus,
          # 
          # ich habe 7 verschiedene Varianzschätzer untersucht.  Alle haben die beweisbare - und numerisch überprüfte - Eigenschaft, dass ihr Erwartungswert gleich der richtigen Varianz des mittleren Log Fold Changes ist.  Einer dieser Schätzer ist die empirische Varianz der K Replikate der Log Fold Changes dividiert durch K.  Ein weiterer Schätzer ist der, den Du vorschlugst, mit korrigiertem Vorfaktor.  Die 5 weiteren Schätzer kannst Du dem - quick and dirty - R-Skript entnehmen.  Oder ich beschreibe Dir sie alle gern beim nächsten Gespräch.
          # 
          # Der aus statistischer Sicht wichtige Unterschied dieser 7 Schätzer ist - bei gleichem und korrektem Erwartungswert - deren Varianz.  Hier bevorzugen wir selbstverständlich den oder die Schätzer mit der geringsten Varianz.  Der in dieser Hinsicht schlechteste der sieben Schätzer - mit der grössten Varianz - ist die empirische Varianz der K Replikate der Log Fold Changes dividiert durch K.  D. h. dieser Schätzer eignet sich als sanity check oder gold standard, aber nicht in der Praxis.  Der beste der 7 Schätzer ist v[1] (s. R-Skript).
          # 
          # v[1] = ((K*M-1) * var (as.vector (xkm)) + (K*N-1) * var (as.vector (ykn))) / (K*M-1 + K*N-1);
          # 
          # D. h. wir berechnen den gemeinsamen x-Mittelwert der K=3 * M=4 x-Werte und die Varianz der K=3 * M=4 x-Werte um diesen gemeinsamen x-Mittelwert.  Wir berechnen den gemeinsamen y-Mittelwert der K=3 * M=2 y-Werte und die Varianz der K=3 * M=4 y-Werte um diesen gemeinsamen y-Mittelwert.  Wir berechnen die gewichtete Summe der beiden Varianzen, wobei die Gewichte gleich K=3 * M=4 - 1 sowie K=3 * M=2 - 1 sind. Und wir dividieren diese gewichtete Summe durch K=3 * M=4 - 1 + K=3 * M=2 - 1.  Dies ist ein erwartungstreuer Schätzer für \sigma^2.
          # 
          # Diesen Schätzwert multiplizieren wir abschliessend mit (1/M + 1/N) / K.  So erhalten wir einen erwartungstreuen Schätzer für die Varianz des mittleren log Fold Changes.  D. h. der Erwartungswert dieses Schätzers ist gleich der richtigen Varianz des mittleren Log Fold Changes.  Von den 7 untersuchten erwartungstreuen Varianzschätzern hat dieser Schätzer die kleinste Varianz.  D. h. die sich daraus ergebenden Fehlerbalken sind beweisbar im Mittel korrekt und ihre quadratischen Abweichungen von den wahren Werten sind - von den 7 untersuchten Varianten - am geringsten.
          # 
          # Diese Variante hat ebenfalls genau wie Deine Variante die Eigenschaft, dass sie auch für K=1 funktioniert, was z. B. für die empirische Varianz der K Replikate der Log Fold Changes dividiert durch K nicht gilt.  D. h. diesen Schätzer können wir gleichermassen für die Schätzung der Fehlerbalken der M+N Array-Werte sowie der K*(M+N) qPCR-Werte verwenden.  Für die Array-Werte - mit K=1 - kommt übrigens genau unser Schätzer von neulich - mit korrigiertem Vorfaktor - raus.  Es wäre super, wenn Du dieses Ergebnis noch einmal überprüfen könntest, damit wir dann diesen Schätzer für die realen Daten nutzen können.
        }
      }
      
      PAPER_PLOT{  
        ### sdt.err of pooled varianz 2016-12-31" 
        #### Array 
        K=1
        xkm = D.n$E[Paper.qgenesIDs,c0] ; M=length(c0)
        ykn = D.n$E[Paper.qgenesIDs,c1] ; N=length(c1)
       
        v1=rep(NA,length(rownames(xkm))); names(v1) <- rownames(xkm)
        for(g in rownames(xkm)){
          v1[g] = ((K*M-1) * var (as.vector (xkm[g,])) + (K*N-1) * var (as.vector (ykn[g,]))) / (K*M-1 + K*N-1);
        }
        
        VarFC_v1 <- (v1 * (1/M + 1/N) / K)
        
        VarFC_v1 <- VarFC_v1[Paper.qgenesIDs]
        names(VarFC_v1) <- names(Paper.qgenesIDs)
        #### qPCR
        ####### relVal 
          qGenes.relVal   <- lapply(qGenes,function(x) x[c(2,4,6)])
          qGenes.relVal.1 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,1]) ); 
          qGenes.relVal.2 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,2]) )
          qGenes.relVal.3 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,3]) )
          colnames(D.n$E) -> colnames(qGenes.relVal.1) -> colnames(qGenes.relVal.2) -> colnames(qGenes.relVal.3) 
          
          qGenes.rel=list()
          for(g in names(Paper.qgenesIDs)){
            qGenes.rel[[g]] <- rbind(qGenes.relVal.1[g,],
                                    qGenes.relVal.2[g,],
                                    qGenes.relVal.3[g,])
          }
          
          qGenes.rel.FC = rep(NA,length(names(qGenes.rel))); names(qGenes.rel.FC) <- names(qGenes.rel)
          for(g in names(qGenes.ct)){
            qGenes.rel.FC[g] <- mean(colMeans(log2(qGenes.rel[[g]][,c1] ))) -  mean(colMeans(log2(qGenes.rel[[g]][,c0] )))
          }
        
          K = nrow(qGenes.rel[[1]])
          M=length(c0)
          N=length(c1)
          v1q.rel=rep(NA,length(names(qGenes.rel))); names(v1q.rel) <- names(qGenes.rel)
          for(g in names(qGenes.rel)){
            
            xkm <- log2(qGenes.rel[[g]][,c0])
            ykn <- log2(qGenes.rel[[g]][,c1])
            
            v1q.rel[g] = ((K*M-1) * var (as.vector (xkm)) + (K*N-1) * var (as.vector (ykn))) / (K*M-1 + K*N-1);
          }
          VarFC_v1q.rel <- (v1q.rel * (1/M + 1/N) / K)
          
        # ###### ctGen.ctRef
        #   qGenes.ctGen.ctRef   <- lapply(qGenes,function(x) x[c(1,3,5)])
        #   qGenes.ctGen.ctRef.1 <- do.call(rbind,lapply(qGenes.ctGen.ctRef,function(x) x[,1]) ); 
        #   qGenes.ctGen.ctRef.2 <- do.call(rbind,lapply(qGenes.ctGen.ctRef,function(x) x[,2]) )
        #   qGenes.ctGen.ctRef.3 <- do.call(rbind,lapply(qGenes.ctGen.ctRef,function(x) x[,3]) )
        #   colnames(D.n$E) -> colnames(qGenes.ctGen.ctRef.1) -> colnames(qGenes.ctGen.ctRef.2) -> colnames(qGenes.ctGen.ctRef.3) 
        #   
        #   qGenes.ct=list()
        #   for(g in names(Paper.qgenesIDs)){
        #     qGenes.ct[[g]] <- rbind(qGenes.ctGen.ctRef.1[g,],
        #                           qGenes.ctGen.ctRef.2[g,],
        #                           qGenes.ctGen.ctRef.3[g,])
        #   }
        #   
        #   qGenes.ct.FC = rep(NA,length(names(qGenes.ct))); names(qGenes.ct.FC) <- names(qGenes.ct)
        #   for(g in names(qGenes.ct)){
        #     qGenes.ct.FC[g] <- mean(colMeans(qGenes.ct[[g]][,c0] )) -  mean(colMeans(qGenes.ct[[g]][,c1] ))
        #   }
        #   cbind(grC.FC,qGenes.ct.FC)
        #   
        #   K = nrow(qGenes.ct[[1]])
        #   M=length(c0)
        #   N=length(c1)
        #   v1q.ct=rep(NA,length(names(qGenes.ct))); names(v1q.ct) <- names(qGenes.ct)
        #   for(g in names(qGenes.ct)){
        #     
        #     xkm <- qGenes.ct[[g]][,c0]
        #     ykn <- qGenes.ct[[g]][,c1]
        #     
        #     v1q.ct[g] = ((K*M-1) * var (as.vector (xkm)) + (K*N-1) * var (as.vector (ykn))) / (K*M-1 + K*N-1);
        #   }
        #   VarFC_v1q.ct <- (v1q.ct * (1/M + 1/N) / K)
        #   VarFC_v1q.ct- VarFC_v1q.rel
        # 
          
        df2$VarFC <- c(VarFC_v1,VarFC_v1q.rel)
        limits <- aes(ymax = df2$value + df2$VarFC, ymin=df2$value - df2$VarFC)
        pdf(paste0(wd,"/List/PAPER_qPCR_FC_new_Limma_2016-12-31__VarFC_IVO.pdf"),height=6,width=8)
          g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
            thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
          g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) 
        dev.off() 
      
        
       
      }
      
      
      
    }
  ##### 
  errorbar2016-12-20{ 
  SF_std.error.xy=c()
  for(i in Paper.qgenes){
    x=log2(paper.1[i,c('SF 767 1','SF 767 3','SF 767 4','SF 767 5')])
    y=log2(paper.1[i,c('SF 767 2','SF 767 6')])
    
    require(plotrix)
    require(psych)
    # describe(x)$se
    # 
    # std.error(x)
    # std.error(y)
    # sd(x, na.rm=TRUE) /  sqrt(length(x)) 
    # sd(y, na.rm=TRUE) /  sqrt(length(y)) 
    std.error.xy = sqrt( std.error(x)^2 +  std.error(y)^2)
    SF_std.error.xy =c(SF_std.error.xy,std.error.xy)
    
    # https://www.researchgate.net/post/Can_anyone_help_with_calculating_error_in_RT-qPCRs_fold-change_data
  }
  
  qP_std.error.xy=c()
  for(i in Paper.qgenes){
    x <- qGenes.M[[i]][,2][c(1,3,4,5)]
    y <- qGenes.M[[i]][,2][c(2,6)]
    
    require(plotrix)
    require(psych)
    std.error.xy = sqrt( std.error(x)^2 +  std.error(y)^2)
    qP_std.error.xy =c(qP_std.error.xy,std.error.xy)
  }
  
  df2$sr=c(SF_std.error.xy,qP_std.error.xy)
  limits <- aes(ymax = df2$value + df2$sr, ymin=df2$value - df2$sr)
  dodge <- position_dodge(width=0.9)
  
  pdf(paste0(wd,"/List/PAPER_qPCR_FC_new_Limma_2016-12-20__errorbar.pdf"),height=6,width=8)
    g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
      thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25)
  dev.off()
  ##### errorbar 2016-12-20" 
  }
    
  OLD_PAPER_PLOT{  
  ### sdt.err of pooled varianz 2016-12-22" 
    #### Array 
      xx = D.n$E[Paper.qgenesIDs,c1]
      yy = D.n$E[Paper.qgenesIDs,c0]
      xxM = rowMeans(xx)
      yyM = rowMeans(yy)
      poolVARa = (           rowSums((xx-xxM)^2) + rowSums((yy-yyM)^2))             / ( length(c1)-1 +length(c0)-1 )
      ### nach IVO 2016-12-22"
                                                                                         # !!!! 4 nicht 5
      poolVAR_ivo = ((length(c1)-1)*apply(xx,1,var) + (length(c0)-1)*apply(yy,1,var) ) / ( length(c1)-1+length(c0)-1 )
      # poolVAR_ivo == poolVARa
      pooledERR_ivo = sqrt(poolVAR_ivo) #/ sqrt( length(c0)+length(c1) )
    
    ### qPCR
      Qxx <- do.call(rbind,lapply(qGenes.M,function(x) x[,2][c1]))
      Qyy <- do.call(rbind,lapply(qGenes.M,function(x) x[,2][c0]))
      
      QpoolVAR_ivo = ((length(c1)-1)*apply(Qxx,1,var) + (length(c0)-1)*apply(Qyy,1,var) ) / ( length(c1)-1+length(c0)-1 )
      QpoolVAR_ivo <- QpoolVAR_ivo[names(Paper.qgenesIDs)]
      QpooledERR_ivo = sqrt(QpoolVAR_ivo) #/ sqrt( length(c0)+length(c1) )
      
      df2$POOLsrIVO =c(pooledERR_ivo,QpooledERR_ivo)
    
      
    #### PAPER PLOT   
      limits <- aes(ymax = df2$value + df2$POOLsrIVO, ymin=df2$value - df2$POOLsrIVO)
        pdf(paste0(wd,"/List/PAPER_qPCR_FC_new_Limma_2016-12-22__POOLEDerrorbar_IVO.pdf"),height=6,width=8)
        g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
          thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
        g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) 
      dev.off()
    
    # pdf(paste0(wd,"/List/PAPER_qPCR_FC_new_Limma_2016-12-22__errorbar_Compare.pdf"),height=6,width=8)
    #   limits <- aes(ymax = df2$value + df2$sr, ymin=df2$value - df2$sr)
    #   g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
    #     thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    #   g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) + labs( title='sqrt( std.error(x)^2 +  std.error(y)^2)')
    #   
    #   limits <- aes(ymax = df2$value + df2$POOLsrIVO, ymin=df2$value - df2$POOLsrIVO)
    #   g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
    #     thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    #   g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) + labs( title='( len(x)-1 *var(x) + len(y)-1*var(y) ) / ( len(x)-1 + len(y)-1) ')
    #   
    #   limits <- aes(ymax = df2$value + df2$POOLsr, ymin=df2$value - df2$POOLsr)
    #   g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
    #     thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    #   g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) + labs( title='( len(x)-1 *var(x) + len(y)-1*var(y) ) / ( [len(x)+len(y)] -1) ')
    # dev.off()
  ### sdt.err of pooled varianz 2016-12-22" 
  }  
  
  OLD_simulation_POOLDEVAR{
    
    ######### simulation ########
    # init 
    mu0 = 12;mu1 = 8 ; pSD = 4
    
    # K = # replicates ; M = # genes of x ; N = # genes of y
    K = 3 ; M=4 ;N=2
    
    c() -> fcs -> sterr -> sterrREP
    for(i in 1:1001){
      c() -> Tmp -> Tx -> Ty
      for(j in 1:K){
        x = rnorm(M,mu0,pSD)
        y = rnorm(N,mu1,pSD)
        
        Tmp[j] = mean(x) - mean(y)       
        Tx = rbind(Tx,x)
        Ty = rbind(Ty,y)
        
      }
      
      ### Mean of replicates 
      xM <- colMeans(Tx) ### K=1, dann ist xM=colMeans(Tx) das selbe wie x 
      yM <- colMeans(Ty)
      
      #### NEU 
        ##### Bin mir nicht sicher ob du mir: '' \sigma^2 die Varianz der x- und y-Werte ist.'' die gepoolte Var of x und y meinst.
        poolVAR = ((M-1)*var(xM) + (N-1)*var(yM) ) / ( M-1+N-1 )
        VarFC = ( 1/K * (1/M + 1/N) ) * poolVAR                     # Varianz des log2 fold changes !NEU!
      #### NEU 
        
      
      sterrREP <- c(sterrREP, sd(Tmp) / sqrt(2) )
      fcs <- c(fcs , mean(xM) - mean(yM) )  

      #### NEU 
        sterr <- c(sterr,  sqrt(VarFC)    ) ## Ist das richtig ???
      #### NEU 
      
    }

    median(sterr)
    sd(fcs)
    median(sterrREP)

    hist(fcs,20)
    hist(sterr,20)
    hist(sterrREP,20)
    
    ######### t-test ########
    #### data 
      grC.FCs <-cbind(c(1.2114816 ,0.7842151 ,1.2203032 ,-1.2364964 ,-1.2509517 ,-0.8778586),c(1.0941976  ,0.8206860  ,1.0850383 ,-0.8431039 ,-0.9681039 ,-0.4347706))
      dimnames(grC.FCs) = list(c('CKAP2L', 'ROCK1', 'TPR', 'ALDH4A1', 'CLCA2', 'GALNS'),c('Microarray' , 'RT-qPCR'))
      grC.FCs <- cbind(grC.FCs,
      c(1.3269167, 1.0385890, 1.4098864, -0.7729494, -0.8479494, -0.2979494), 
      c(1.1427910, 0.9497030, 0.9350965, -0.7234841, -0.8984841, -0.3234841), 
      c(0.8128850, 0.4737661, 0.9101322, -1.0328784, -1.1578784, -0.6828784) )
      colnames(grC.FCs)[3:5] <- c( 'RT-qPCR_R1','RT-qPCR_R2','RT-qPCR_R3')
      
    plot(Microarray ~ RT.qPCR, data = data.frame(grC.FCs))
    abline(lm(Microarray ~ RT.qPCR, data = data.frame(grC.FCs)),col="red")
    # summary <- summary(lm(Microarray ~ RT.qPCR, data = data.frame(grC.FCs)))
    # adjRsq <- summary$adj.r.squared
    # fStat <- summary$fstatistic
    # pValue <- pf(fStat[1], fStat[2], fStat[3], lower.tail = F)
    p <- cor.test(grC.FCs[,1],grC.FC[,2])$p.value
    
    
    ## hattest du das im Kopf beim T-test ?
    t.test(grC.FCs['ROCK1',1 ], grC.FCs['ROCK1',c(3:5) ], var.equal = T)
    
    tt <- sapply(rownames(grC.FCs), function(y) {
      t.test(grC.FCs[y,1 ], grC.FCs[y,c(3:5) ], var.equal = T)$p.val
    })
    tt
    
    t.test(grC.FCs['ROCK1',1 ], grC.FCs['ROCK1',2 ], paired = T,var.equal = T)
    
    
    OLD{
    fcs <- c()
    sterr <- c()
    sterrREP <- c()
    for(i in 1:1001){
      c() -> Tmp -> Tx -> Ty
      for(j in 1:3){
        x = rnorm(4,mu0,pSD)
        y = rnorm(2,mu1,pSD)
        
        Tmp[j] = mean(x) - mean(y)       
        Tx = rbind(Tx,x)
        Ty = rbind(Ty,y)
        
      }
      xM <- colMeans(Tx)
      yM <- colMeans(Ty)
      
      # TMPfc <- mean(xM) - mean(yM)
      # TMPfc2 <- mean(Tmp)
      # TMPfc - TMPfc2
      
      sterrREP <- c(sterrREP, sd(Tmp) / sqrt(2) )
      
      fcs <- c(fcs , mean(xM) - mean(yM) )  
      sterr <- c(sterr, sqrt( ( (length(xM)-1)*var(xM) + (length(yM)-1)*var(yM) ) / ( length(xM)-1 + length(yM)-1) )   )
      
    }
    median(sterr)
    # mean(sterr)
    sd(fcs)
    
    median(sterrREP)
    # mean(sterrREP)
    
    hist(fcs)
    hist(sterr,20)
    hist(sterrREP,20)
    }
    
  }
  
  ########## test mit replicates einzeln  
  tt{  
  c1 
  qGenes.relVal <- lapply(qGenes,function(x) x[c(2,4,6)])
  qGenes.relVal.1 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,1]) ); 
  qGenes.relVal.2 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,2]) )
  qGenes.relVal.3 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,3]) )
  colnames(D.n$E) -> colnames(qGenes.relVal.1) -> colnames(qGenes.relVal.2) -> colnames(qGenes.relVal.3) 
  
  log2FC.qGenes.relVal.1 <- rowMeans(log2(qGenes.relVal.1[,c1])) - rowMeans(log2(qGenes.relVal.1[,c0]))
  log2FC.qGenes.relVal.2 <- rowMeans(log2(qGenes.relVal.2[,c1])) - rowMeans(log2(qGenes.relVal.2[,c0]))
  log2FC.qGenes.relVal.3 <- rowMeans(log2(qGenes.relVal.3[,c1])) - rowMeans(log2(qGenes.relVal.3[,c0]))
  
  log2FC.qGenes.relVal <- rbind(log2FC.qGenes.relVal.1,log2FC.qGenes.relVal.2,log2FC.qGenes.relVal.3)
  
  cbind(grC.FC, 'meanQPCR' = colMeans(log2FC.qGenes.relVal)[names(Paper.qgenesIDs)])
  
  df2<-melt(grC.FC)
  rns <- levels(df2$X1)[c(2,5,6,1,3,4)] ;df2$X1 <- factor(df2$X1, levels = rns)
  cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6],'green')

  df3 <-  df2 
  df3 <- rbind(df3,df3[7:12,])
  df3$POOLsrIVO[13:18] <- ( (apply(log2FC.qGenes.relVal,2,sd)[names(Paper.qgenesIDs)]) / sqrt(3) )
  levels(df3$X2) <- c(levels( df3$X2),'NEW')
  df3$X2[13:18] <- factor('NEW')
  
  limits <- aes(ymax = df3$value + df3$POOLsrIVO, ymin=df3$value - df3$POOLsrIVO)
  pdf(paste0(wd,"/List/IVO.pdf"),height=6,width=8)
  g1<-ggplot(df3,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-2,2))+ scale_fill_manual(values = cols,name="") +
    thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
  g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) 
  dev.off()
  }
  
  }
  
  2017-01-05__Pooled_standard_error{
    ### IDEA !!! http://wolfweb.unr.edu/~ldyer/classes/396/PSE.pdf --- use  Satterthwaite Approximation
    
    
    ## save.image("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2017-01-05.RData")
    load("~/Kappler/Kappler_Wichmann_Medizin/Auswertung/WS2017-01-05.RData")
    
    #### Array 
      K=1
      xkm = D.n$E[Paper.qgenesIDs,c0] ; M=length(c0)
      ykn = D.n$E[Paper.qgenesIDs,c1] ; N=length(c1)
      
      v1=rep(NA,length(rownames(xkm))); names(v1) <- rownames(xkm)
      for(g in rownames(xkm)){
        v1[g] = ((K*M-1) * var (as.vector (xkm[g,])) + (K*N-1) * var (as.vector (ykn[g,]))) / (K*M-1 + K*N-1);
      }
    
    # Pooled standard error PSR  ## https://onlinecourses.science.psu.edu/stat200/node/60
      ## https://classroom.udacity.com/courses/ud201/lessons/1330208559/concepts/1548634750923#
      
      PSRv1 = sqrt(v1) * sqrt( (1/(M*K)) + (1/(N*K)) )  #### sqrt( (v1 / (M*K) ) + (v1 / (N*K)) )
      df=(M*K)+(N*K)-2
      
      xm = rowMeans(xkm)
      yn = rowMeans(ykn)
      FC <- yn -xm
      t= ( FC ) / PSRv1
    
    # Confidence Intervals
      alpha=0.05
      E = qt( (1-alpha/2), df=df)*PSRv1
      ConfInt <- matrix( (yn -xm) + c(-E, E) ,length(rownames(xkm)), dimnames =list(rownames(xkm) ,c('ConfIntLow','ConfIntUp')) )
    # t.test(ykn[g,],xkm[g,],var.equal = T)
    
    gg <- data.frame(
      "Gene" = names(Paper.qgenesIDs),
      "FC" = yn -xm,
      "E" = E,
      "PSR" = PSRv1,
      "Set" ="Microarray"
    )

    # limits <- aes(ymax = gg$FC + gg$E, ymin=gg$FC - gg$E)
    # dodge <- position_dodge(width=0.9)
    # g1<-ggplot(gg,aes(x=factor(Gene),y=FC,fill=Gene)) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change") + ylim(c(-3,3))
    # # + scale_fill_manual(values = cols,name="") + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    # g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) 


    #### qPCR
    ####### relVal 
      qGenes.relVal   <- lapply(qGenes,function(x) x[c(2,4,6)])
      qGenes.relVal.1 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,1]) ); 
      qGenes.relVal.2 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,2]) )
      qGenes.relVal.3 <- do.call(rbind,lapply(qGenes.relVal,function(x) x[,3]) )
      colnames(D.n$E) -> colnames(qGenes.relVal.1) -> colnames(qGenes.relVal.2) -> colnames(qGenes.relVal.3) 
      
      qGenes.rel=list()
      for(g in names(Paper.qgenesIDs)){
        qGenes.rel[[g]] <- rbind(qGenes.relVal.1[g,],
                                 qGenes.relVal.2[g,],
                                 qGenes.relVal.3[g,])
      }
      
      qGenes.rel.FC = rep(NA,length(names(qGenes.rel))); names(qGenes.rel.FC) <- names(qGenes.rel)
      for(g in names(qGenes.rel)){
        qGenes.rel.FC[g] <- mean(colMeans(log2(qGenes.rel[[g]][,c1] ))) -  mean(colMeans(log2(qGenes.rel[[g]][,c0] )))
      }
      
      K = nrow(qGenes.rel[[1]])
      M=length(c0)
      N=length(c1)
      v1q.rel=rep(NA,length(names(qGenes.rel))); names(v1q.rel) <- names(qGenes.rel)
      for(g in names(qGenes.rel)){
        
        xkm <- log2(qGenes.rel[[g]][,c0])
        ykn <- log2(qGenes.rel[[g]][,c1])
        
        v1q.rel[g] = ((K*M-1) * var (as.vector (xkm)) + (K*N-1) * var (as.vector (ykn))) / (K*M-1 + K*N-1);
      }

    # Pooled standard error PSR  ## https://onlinecourses.science.psu.edu/stat200/node/60
      PSRv1q.rel = sqrt(v1q.rel) * sqrt( (1/(M*K)) + (1/(N*K)) )
      dfv1q.rel=(M*K)+(N*K)-2
      
      xm <- do.call(c,lapply(qGenes.rel , function(x) mean(as.vector(log2(x[,c0])))))
      yn <- do.call(c,lapply(qGenes.rel , function(x) mean(as.vector(log2(x[,c1])))))
      FC <- yn -xm
      t <- FC / PSRv1q.rel
      
    # 95% Confidence Interval
      alpha = 0.05
      E = qt( (1-alpha/2), df=dfv1q.rel)*PSRv1q.rel
      ConfInt <- matrix( (yn -xm) + c(-E, E) ,length(names(Paper.qgenesIDs)), dimnames =list(names(Paper.qgenesIDs) ,c('ConfIntLow','ConfIntUp')) )
      
      t[g]
      ConfInt[g,]
      t.test(as.vector(log2(qGenes.rel[[g]][,c1])),as.vector(log2(qGenes.rel[[g]][,c0])),var.equal = T)
      
      ggQ <- data.frame(
        "Gene" = names(Paper.qgenesIDs),
        "FC" = FC,
        "E" = E,
        "PSR" = PSRv1q.rel,
        "Set" ="qPCR"
      )
    
      
    ######### PLOT   
    TMP <- rbind(gg,ggQ)
    rns <- levels(TMP$Gene)[c(2,5,6,1,3,4)] ;TMP$Gene <- factor(TMP$Gene, levels = rns)
    cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6])
    
    thememap <- function (base_size = 12,legend_key_size=0.4, base_family = "") {
      theme_gray(base_size = base_size, base_family = base_family) %+replace% 
        theme(title = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
              axis.title.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
              axis.text.x  = element_text(face="bold", colour=1,angle=0, vjust=0.5, size=base_size),
              strip.text.x = element_text(face="bold", colour=1,angle=0  ,vjust=0.5, size=base_size),
              axis.title.y = element_text(face="bold", colour=1,angle=90 ,vjust=0.5,hjust=.5, size=base_size),
              axis.text.y  = element_text(face="bold", colour=1,angle=0  ,vjust=0.0, size=base_size),
              panel.background = element_rect(fill="white"),
              #panel.grid.minor.y = element_line(size=3),
              panel.grid.major = element_line(colour = "white"),
              legend.key.size = unit(legend_key_size, "cm"),
              legend.text = element_text(face="bold" ,colour=1,angle=0  ,vjust=0.0, size=base_size),
              legend.title = element_text(face="bold",colour=1,angle=0  ,vjust=-0.8, size=base_size)           
        )
    }
    
    
    pdf(paste0(wd,"/List/2017-01-05_PSE_ConfInter95.pdf"),height=6,width=8)
    
      limits <- aes(ymax = TMP$FC + TMP$E , ymin=TMP$FC - TMP$E)
      dodge <- position_dodge(width=0.9)
      g1<-ggplot(TMP,aes(x=factor(Gene),y=FC,fill=factor(Set))) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change") + ylim(c(-2.5 , 2.5)) +
         scale_fill_manual(values = cols,name="") + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
      g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25)+ ggtitle("95% Confidence Interval")
      
      limits <- aes(ymax = TMP$FC + TMP$PSR , ymin=TMP$FC - TMP$PSR)
      dodge <- position_dodge(width=0.9)
      g1<-ggplot(TMP,aes(x=factor(Gene),y=FC,fill=factor(Set))) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change") + ylim(c(-2.5 , 2.5)) +
        scale_fill_manual(values = cols,name="") + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
      g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25)+ ggtitle("Pooled standard error")
   
    dev.off()  
    
    ##### use in PAPER 2017-01-27
    pdf(paste0(wd,"/List/2017-01-27_PSE.pdf"),height=6,width=8)
    require(ggplot2)
    
    limits <- aes(ymax = TMP$FC + TMP$PSR , ymin=TMP$FC - TMP$PSR)
    dodge <- position_dodge(width=0.9)
    
    g1<-ggplot(TMP,aes(x=factor(Gene),y=FC,fill=factor(Set))) + geom_bar(position="dodge",stat="identity",width=.8)+labs(x = "", y = "Log2-fold change") + ylim(c(-2.5 , 2.5)) +
      scale_fill_manual(values = cols,name="") + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
    g1 + geom_bar(position="dodge", stat="identity")+ geom_errorbar(limits, position=dodge, width=0.25) #+ ggtitle("Pooled standard error")
    
    dev.off()  
    save.image(paste0(wd,"/List/2017-01-27_PSE.Rdata"))
    ##### use in PAPER 2017-01-27
    
    
  }
  
  
  
  Posterior_for_qPCR{ #2016_09_22
  
    get.logLs.qPCR<-function(qq){
        # qq= t(qGenes$TPR[c(1,3,5)])
      N<-c()
      a0<-c(1,2,3,4,5,6);   N[1]  <- length(a0);
      b0<-c(1,3,5);         N[2] <- length(b0);
      b1<-c(2,4,6);         N[3] <- length(b1);
      c0<-c(1,3,5,4);       N[4] <- length(c0);
      c1<-c(2,6);           N[5] <- length(c1);
      d0<-c(1,3,5,4,6);     N[6] <- length(d0);
      d1<-c(2);             N[7] <- length(d1);
      sub<-list(a0=a0,b0=b0,b1=b1,c0=c0,c1=c1,d0=d0,d1=d1)
      
      get.var.qq<-function(d,sub,MUs,N,i0,i1){
        # N=N3;i0=4;i1=5
        return(  ( sum((qq[,sub[[i0]]]-MUs[i0])^2) + sum((qq[,sub[[i1]]]-MUs[i1])^2)) / ( N[i0]+N[i1]-1 ) ) 
      }
      
      
      calc.MUs.qq<-function(qq,sub){  
        MUs<-c()
        for(i in 1:length(sub)){
          MUs[i]<-mean(qq[,sub[[i]]])
                  # mean(colMeans(qq[,sub[[i]]]))
        }
        return(MUs)
      }
     
      MUs <- calc.MUs.qq(qq,sub)
      
      VARs<-c()
      N3=N*3
      VARs[1]<- 1/(N3[1]-1)*sum( (qq[,a0]-MUs[1])^2)
      VARs[2]<-get.var.qq(qq,sub,MUs,N3,2,3)
      VARs[3]<-get.var.qq(qq,sub,MUs,N3,4,5)
      VARs[4]<-get.var.qq(qq,sub,MUs,N3,6,7)
      
      # d=colMeans(qq)
      # VARs2<-c()
      # VARs2[1]<- 1/(N[1]-1)*sum( (d[a0]-MUs[1])^2)
      # VARs2[2]<-get.var(d,sub,MUs,N,2,3)
      # VARs2[3]<-get.var(d,sub,MUs,N,4,5)
      # VARs2[4]<-get.var(d,sub,MUs,N,6,7)
      
      
      normprob = function (x,m,v) { (1/sqrt(2*pi*v))*exp(-((x-m)^2)/(2*v))  }   
      
      logL2<-c()  
      logL2[1] <- log(prod(c(normprob(qq[,sub[[1]]][,1],MUs[1],VARs[1]),
                             normprob(qq[,sub[[1]]][,2],MUs[1],VARs[1]),
                             normprob(qq[,sub[[1]]][,3],MUs[1],VARs[1]),
                             normprob(qq[,sub[[1]]][,4],MUs[1],VARs[1]),
                             normprob(qq[,sub[[1]]][,5],MUs[1],VARs[1]),
                             normprob(qq[,sub[[1]]][,6],MUs[1],VARs[1]))))
      
      logL2[2] <- log(prod(c(normprob(qq[,sub[[2]]][,1],MUs[2],VARs[2]),
                             normprob(qq[,sub[[2]]][,2],MUs[2],VARs[2]),
                             normprob(qq[,sub[[2]]][,3],MUs[2],VARs[2]),
                             normprob(qq[,sub[[3]]][,1],MUs[3],VARs[2]),
                             normprob(qq[,sub[[3]]][,2],MUs[3],VARs[2]),
                             normprob(qq[,sub[[3]]][,3],MUs[3],VARs[2]))))
      
      logL2[3] <- log(prod(c(normprob(qq[,sub[[4]]][,1],MUs[4],VARs[3]),
                             normprob(qq[,sub[[4]]][,2],MUs[4],VARs[3]),
                             normprob(qq[,sub[[4]]][,3],MUs[4],VARs[3]),
                             normprob(qq[,sub[[4]]][,4],MUs[4],VARs[3]),
                             normprob(qq[,sub[[5]]][,1],MUs[5],VARs[3]),
                             normprob(qq[,sub[[5]]][,2],MUs[5],VARs[3]))))
      
      logL2[4] <- log(prod(c(normprob(qq[,sub[[6]]][,1],MUs[6],VARs[4]),
                             normprob(qq[,sub[[6]]][,2],MUs[6],VARs[4]),
                             normprob(qq[,sub[[6]]][,3],MUs[6],VARs[4]),
                             normprob(qq[,sub[[6]]][,4],MUs[6],VARs[4]),
                             normprob(qq[,sub[[6]]][,5],MUs[6],VARs[4]),
                             normprob(qq[,sub[[7]]][1:3],MUs[7],VARs[4]))))
      
      
      # barplot(exp(logL2))
      # barplot(exp(logL))
      # ### Test for 
      # x=D.n$E
      # g="ILMN_1730999" #TPR
      # d=x[g,]
      # logL2
      # [1] -5.534718 -3.177679  1.601310 -4.269458
      # > logL
      # [1] -5.534718 -3.177679  1.601310 -4.269458
      #
        return(logL2)
      }
    
    logL2.qq=c()
    for(i in names(qGenes)){
      # qq=t(qGenes[[i]][c(1,3,5)])
      qq=t(qGenes[[i]][c(2,4,6)])
      logL2.qq = rbind(logL2.qq,get.logLs.qPCR(qq))
    };rownames(logL2.qq) = names(qGenes);colnames(logL2.qq) <- c('a','b','c','d')

    logL2.qq[names(qGenes) [c(1,2,4,6,9,11)],]
    qq.BIC <- make.IC(logL2.qq,c(2,3),18) # BIC with 6 = number of values   
    qq.BIC.Post<-get.posterior(qq.BIC,Pis = c(0.7,0.1,0.1,0.1))
    m<-apply((qq.BIC.Post),MARGIN=1,FUN=maxIndex)

    apply(logL2.qq,1,which.max)
    
    dq=c()
    for(i in names(qGenes)){
      # qq=colMeans(t(qGenes[[i]][c(1,3,5)]))
      #qq=colMeans(t(qGenes[[i]][c(2,4,6)]))
      qq=t(rowMeans(log2(qGenes[[i]][c(2,4,6)])))
      dq=rbind(dq,qq)
    }  ;rownames(dq) = names(qGenes)
    D.n_L.dq<-get.ML.logLs.var.NEW(dq)
    D.n_L.dq.BIC <- make.IC(D.n_L.dq,c(2,3),6) # BIC with 6 = number of values   
    D.n_L.dq.BIC$IC
    D.n_L.dq.BIC.Post<-get.posterior(D.n_L.dq.BIC,Pis = c(0.7,0.1,0.1,0.1))
    m<-apply((D.n_L.dq.BIC.Post),MARGIN=1,FUN=maxIndex)
    D.n_L.dq.BIC.Post[names(qGenes) [c(1,2,4,6,9,11)],]    
  
    qq.BIC.Post[names(qGenes) [c(1,2,4,6,9,11)],]    
    
    D.n_L.BIC.Post[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],]
    D.n_L.BIC$IC[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],]
    
    dq 
    D$E[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],]
    
    pdf(paste0(wd,"/List/Posterior_for_qPCR_2016_10_19.pdf"))
    cols<-c(brewer.pal(9,"RdGy")[8],brewer.pal(9,"Reds")[6])
    barplot(rbind( D.n_L.dq.BIC.Post[names(qGenes) [c(1,2,4,6,9,11)],1],
                   D.n_L.BIC.Post[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],1]  ),main='posterior a',beside = T,col=cols,ylim=c(0,1),legend.text = c("qPCR","microArray")) 
   
    barplot(rbind( D.n_L.dq.BIC.Post[names(qGenes) [c(1,2,4,6,9,11)],2],
                   D.n_L.BIC.Post[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],2]  ),main='posterior b',beside = T,col=cols,ylim=c(0,1),legend.text = c("qPCR","microArray")) 
  
    barplot(rbind( D.n_L.dq.BIC.Post[names(qGenes) [c(1,2,4,6,9,11)],3],
                   D.n_L.BIC.Post[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],3]  ),main='posterior c',beside = T,col=cols,ylim=c(0,1),legend.text = c("qPCR","microArray")) 
    
    barplot(rbind( D.n_L.dq.BIC.Post[names(qGenes) [c(1,2,4,6,9,11)],4],
                   D.n_L.BIC.Post[Paper.qgenesIDs[names(qGenes) [c(1,2,4,6,9,11)]],4]  ),main='posterior d',beside = T,col=cols,ylim=c(0,1),legend.text = c("qPCR","microArray")) 
    
   dev.off()
    
    D.n_L.dq.BIC.Post
    
  # qGenes.M<-c()
  # for(i in 1:length(qgenes)){
  #   qGenes.M[[i]]<-cbind(rowMeans(log2(qGenes[[i]][c(1,3,5)])),rowMeans(log2(qGenes[[i]][c(2,4,6)])),  log2(rowMeans(qGenes[[i]][c(2,4,6)])) )
  #   colnames(qGenes.M[[i]])<-c("log2geoMeanCT","log2geoMeanREL","log2MeanREL")
  # }
  # names(qGenes.M)<-qgenes
  # Paper.qgenes.FC<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5)])
  # 
  }
}

TEEMU{
#### TEEMU
D.n_L.BIC.Post.TEEMU<-get.posterior(D.n_L.BIC,Pis = c(.25,.25,.25,.25))
m.TEEMU<-apply((D.n_L.BIC.Post.TEEMU),MARGIN=1,FUN=maxIndex)
print(table(m.TEEMU))
print(table(m))


t<-D.n_L.BIC.Post.TEEMU[m.TEEMU==1,] ;D.n_L.BIC.Post_a.TEEMU<-t[order(t[,1],decreasing=T),] 
t<-D.n_L.BIC.Post.TEEMU[m.TEEMU==2,] ;D.n_L.BIC.Post_b.TEEMU<-t[order(t[,2],decreasing=T),]
t<-D.n_L.BIC.Post.TEEMU[m.TEEMU==3,] ;D.n_L.BIC.Post_c.TEEMU<-t[order(t[,3],decreasing=T),]
t<-D.n_L.BIC.Post.TEEMU[m.TEEMU==4,] ;D.n_L.BIC.Post_d.TEEMU<-t[order(t[,4],decreasing=T),]
D.n_L.BIC.Post_a.85.TEEMU<-D.n_L.BIC.Post_a.TEEMU[D.n_L.BIC.Post_a.TEEMU[,1]>0.85,]
D.n_L.BIC.Post_b.85.TEEMU<-D.n_L.BIC.Post_b.TEEMU[D.n_L.BIC.Post_b.TEEMU[,2]>0.85,]
D.n_L.BIC.Post_c.85.TEEMU<-D.n_L.BIC.Post_c.TEEMU[D.n_L.BIC.Post_c.TEEMU[,3]>0.85,]
D.n_L.BIC.Post_d.85.TEEMU<-D.n_L.BIC.Post_d.TEEMU[D.n_L.BIC.Post_d.TEEMU[,4]>0.85,]


sapply(f.input2(rownames(D.n_L.BIC.Post_a.TEEMU),rownames(D.n_L.BIC.Post_a)),length)/nrow(D.n_L.BIC.Post_a)
sapply(f.input2(rownames(D.n_L.BIC.Post_b.TEEMU),rownames(D.n_L.BIC.Post_b)),length)/nrow(D.n_L.BIC.Post_b)
sapply(f.input2(rownames(D.n_L.BIC.Post_c.TEEMU),rownames(D.n_L.BIC.Post_c)),length)/nrow(D.n_L.BIC.Post_a)
sapply(f.input2(rownames(D.n_L.BIC.Post_d.TEEMU),rownames(D.n_L.BIC.Post_d)),length)/nrow(D.n_L.BIC.Post_d)

cbind(m.TEEMU[unlist(qgenesIDs)],m[unlist(qgenesIDs)],names(unlist(qgenesIDs)))
cbind(m.TEEMU[Paper.qgenesIDs],m[Paper.qgenesIDs],Paper.qgenesIDs)



#### 
# In usedate() Tue Aug  2 16:31:42 2016" 2016_08_02
pdf(paste0(wd,"/List/PAPER_compare_Pis_2016_08_02.pdf"),height=12,width=8)

par(mfrow= c(2,1))
  PIS = rbind( c(.7,.1,.1,.1),
               c(.25,.25,.25,.25),
              c(.9, .01,.08,.01),
              c(.8, .05,.1,.05),
              c(.7, .1,.2,.1) )
  for(k in 1:nrow(PIS)){
  pis=PIS[k,]
  piTMP<-get.posterior(D.n_L.BIC,Pis = pis)
  m.piTMP<-apply((piTMP),MARGIN=1,FUN=maxIndex)
  a = cbind( unlist(table(m.piTMP)) , dim(D.n_L.BIC$IC)[1] * pis) 
  # m.piTMP[Paper.qgenesIDs]
  # m.piTMP[RNs]
  # piTMP[RNs,]
  a
  a /dim(D.n_L.BIC$IC)[1]
  
  qq1 = density(piTMP[m.piTMP==1,1])
  qq2 = density(piTMP[m.piTMP==2,2])
  qq3 = density(piTMP[m.piTMP==3,3])
  qq4 = density(piTMP[m.piTMP==4,4])
  plot(qq1,xlim=c(0.25,1),ylim=c(0,max( max(qq1$y),max(qq2$y),max(qq3$y),max(qq4$y))) , main = paste(c('a','b','c','d') , pis,round(a[,1] /dim(D.n_L.BIC$IC)[1],2),sep=" => "))
  lines(qq2 ,col=2,lwd=2)
  lines(qq3 ,col=3,lwd=2)
  lines(qq4 ,col=6,lwd=2)
  legend('top',fill =c(1,2,3,6),c('a','b','c','d'))
  
  mat = matrix(0, 4,4 )
  colnames(mat)= c("1","2","3","4")
  for( i in 1:4){
  tmpTab=table(apply(cbind('a'=ALL.SDs[rownames(piTMP[m.piTMP==i,]),1],
                        'b'=ALL.SDs[rownames(piTMP[m.piTMP==i,]),2],
                        'c'=ALL.SDs[rownames(piTMP[m.piTMP==i,]),3] ,
                        'd'=ALL.SDs[rownames(piTMP[m.piTMP==i,]),4]), 1,function(x) which.min(x) ))
  mat[i,names(tmpTab)]= tmpTab
  }
  barplot(t(mat),beside =T , names.arg = c('a','b','c','d') , legend.text = c('a','b','c','d'), main = "number of genes with smallest SD")
  }
  dev.off()
  
allgemeiner__DENSITY__PLot{
  RNs = rownames(piTMP[m.piTMP==4,])
  RNs = rownames(piTMP[m.piTMP==1,])
  RNs = rownames(piTMP[m.piTMP==2,])
  RNs = rownames(piTMP[m.piTMP==3,])
  
  RNs = as.character(STAT3[1:3,1])
  RNs = as.character(ERBB4[1,1])
  for(IDs in RNs){
    # IDs='ILMN_1782788'
    print(IDs)
    pdf(paste0(wd,"/List/DENS/NV_density_ERBB4",IDs,".pdf"),height=6,width=8)
    
    pg.M=ALL.MUs[IDs,]
    pg.s=ALL.SDs[IDs,]
    
    logE = D.n$E[IDs,]
    tmp = data.frame(logE=rep(logE,4),
                     density = rep( 0, length(rep(logE,4))),
                     group    =c(rep('a',6),rep('b',6),rep('c',6),rep('d',6)),
                     indicator= factor(c( rep(0,6) , c(0,1,0,1,0,1), c(0,1,0,0,0,1) ,c(0,1,0,0,0,0) ))
    )
    
    x <- seq(4, 20, length=1000)
    dM=round(max(
      max(dnorm(x,mean =  pg.M['a0'], sd = pg.s['a'])),
      max(dnorm(x,mean =  pg.M['b1'], sd = pg.s['b'])),
      max(dnorm(x,mean =  pg.M['c1'], sd = pg.s['c'])),
      max(dnorm(x,mean =  pg.M['d1'], sd = pg.s['d'])) )+0.5) 
      
    col2=brewer.pal(8,'Paired')[c(2,6)] 
    g2 = ggplot(tmp, aes(x=logE,y=density,colour=indicator))+scale_y_continuous(limits = c(0, dM),breaks = seq(0,dM,1))+scale_x_continuous(limits = c(4.5, 15),breaks = seq(4.5,15,1)) + geom_point(shape=20,size =0)  +facet_wrap(~ group,nrow = 4)  + scale_color_manual(values =col2,name=" Indicator variable",labels = list("g = 0","g = 1" ))+labs(title=IDs,y = "Density", x = "Logarithmic expression levels")
    g3 =g2+ with(tmp[tmp$group=="a",],stat_function(data=tmp[tmp$group=="a" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['a0'], sd = pg.s['a']), size=1.2) ) +
      with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['b0'], sd = pg.s['b']), size=1.2) ) +
      with(tmp[tmp$group=="b",],stat_function(data=tmp[tmp$group=="b" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M['b1'], sd = pg.s['b']), size=1.2) ) +
      with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['c0'], sd = pg.s['c']), size=1.2) ) +
      with(tmp[tmp$group=="c",],stat_function(data=tmp[tmp$group=="c" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M['c1'], sd = pg.s['c']), size=1.2) ) +
      with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==0,],fun=dnorm,args=list(mean =  pg.M['d0'], sd = pg.s['d']), size=1.2) ) +
      with(tmp[tmp$group=="d",],stat_function(data=tmp[tmp$group=="d" & tmp$indicator ==1,],fun=dnorm,args=list(mean =  pg.M['d1'], sd = pg.s['d']), size=1.2) )
    g4 = g3 + thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom") 
    g4 =  g4 + geom_point(data=tmp, aes(x=logE,y=density),shape=20,size =3.5) +
      geom_point(data=tmp, aes(x=logE,y=density),shape=21,size =3,colour="black")
    print(g4)
    dev.off()
    
  }  
  
  STAT3
  
  pg.s[RNs,]
  RNs
  
}



#### TEEMU 

qu.bg <- read.delim( paste0(wd,"/../KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_quantile_bkgd_SF767-1-799-6.txt"),sep="\t")
FILE<-qu.bg
FILE.N<-"qu_bkgd"
path="~/Promotion/Kappler_Wichmann_Medizin/Auswertung/NEW_BG/"
IDs<-FILE[,c(1:2,54,55)]
rownames(IDs)<-as.character(FILE[,1])
head(IDs)
EGFR_genes<-c("ABL1","ABL2","ACTA1","ADRBK1","ADRM1","AGTR1","AHNAK","AHSA1","AIP","AKAP12","ALCAM","ALDOA","AMH","ANKS1A","ANKS1B","ANXA1","APBA3","APBB1","APBB2","APBB3","APPL1","AR","AREG","ARF4","ARF6","ATIC","ATP1A1","ATP1B1","ATP5C1","BLK","BMX","BTC","CALM1","CALM2","CALM3","CAMK2A","CAMK2G","CAMLG","CASP1","CAV1","CAV2","CAV3","CBL","CBLB","CBLC","CBR1","CD44","CD59","CD82","CDC25A","CDH1","CEACAM1","CEBPB","CEND1","CFL1","CISH","CLNK","CLTA","CLTCL1","CMTM8","CNTN2","COL9A3","COX2","CRK","CSRP1","CTNNB1","CTNND1","DCN","DCTN2","DEGS1","DNAJC4","DOK2","DOK4","DOK5","DOK6","DUSP3","DYNC1H1","EEF1G","EGF","ELF3","EPB41","EPHA2","EPPK1","EPS15","EPS8","ERBB2","ERBB3","ERBB4","EREG","ERRFI1","ESR1","EZR","FAH","FAS","FBXO25","FER","FES","FKBP4","FKBP8","FN1","GAB1","GABARAPL2","GAPDH","GNAI2","GOT1","GPM6B","GRB10","GRB14","GRB2","GRB7","HBEGF","HEXIM1","HGS","HIST3H3","HOXC10","HSP90B1","HSPA1A","HSPA4","HSPA8","HTT","ICAM1","INPPL1","IRS1","IRS4","ITGA5","JAK2","JUP","KCTD9P2","KRT17","KRT18","KRT7","KRT8","LYN","MAP2K1","MAP3K12","MAP3K14","MAP4K1","MAPK1","MAPK14","MAPK8IP1","MAPK8IP2","MAST1","MATR3","MDH1","MEGF6","MET","MGARP","MOB4","MUC1","NCAM1","NCK1","NCK2","NDN","NEDD4","NUMB","NUMBL","OAZ1","OLFM1","PAK1","PCNA","PDCD6IP","PDGFRB","PDK3","PIK3C2A","PIK3C2B","PIK3R1","PIK3R2","PIK3R3","PIN4","PITPNA","PKIA","PLCG1","PLCG2","PLD1","PLD2","PLEC","PLSCR1","PPP5C","PRCC","PRDX1","PRKACA","PRKAR1A","PRKAR1B","PRKCA","PRKCB","PRKCD","PRKD1","PSMA7","PSMD4","PTGDS","PTK2","PTK2B","PTK6","PTPN1","PTPN11","PTPN12","PTPN2","PTPN6","PTPRJ","PTPRS","RAB3A","RAP1GDS1","RASA1","RGS16","RIN2","RIPK1","RNF115","RNF126","RQCD1","RUSC2","S100A7","S100A9","SCAMP1","SCAMP3","SEC13","SEPP1","SFN","SGSM2","SH2B1","SH2B3","SH2D1A","SH2D2A","SH2D3A","SH3BGRL3","SHC1","SHC2","SHC3","SLA","SLC3A2","SLC9A3R1","SMURF2","SNRPD2","SNX1","SNX2","SNX4","SNX6","SOCS1","SOCS3","SOCS4","SOCS5","SOCS7","SORBS2","SOS1","SOS2","SPARCL1","SPCS2","SRC","STAM2","STAT1","STAT3","STAT5A","STAT5B","STIP1","STUB1","SYK","TGFA","TJP1","TLN1","TMCO3","TNC","TNK2","TNS4","TPI1","TPM1","TRPV1","TUBA1A","TUBA4A","UBB","UBE2V2","UCHL1","UROD","VAPA","VAV1","VAV2","VAV3","XRCC6","YWHAZ","ZAP70","ZNF259")
# from http://biograph.be/concept/show/C1414313
Paper_EGFR_gene <- read.delim(paste0(wd,"/Paper_EGFR_gene.txt"))
EGFR_genes2<-as.character(Paper_EGFR_gene[,1])
#http://en.wikipedia.org/wiki/Epidermal_growth_factor_receptor
wiki<-c("AR","ARF4","CAV1","CAV3","CBL","CBLB","CBLC","CDC25A","CRK","CTNNB1","DCN","EGF","GRB14","Grb2","JAK2","MUC1","NCK1","NCK2","PKCA","PLCG1","PLSCR1","PTPN1","PTPN11","PTPN6","PTPRK","SH2D3A","SH3KBP1","SHC1","SOS1","Src","STAT1","STAT3","STAT5A","UBC","WAS")
#http://hintdb.hgc.jp/htp/proteins/P00533.html
HitPredict<-c("1433Z","ERBB2","CBL","EGF","P3C2B","TGFA","STAM2","RIN1","EGFR","CRK","DUS3","MK14","ARF4","M4K1","PLS1","PTN1","CASP1","PTPRS","P85B","PAK1","ACTS","MET","CEACAM1","TNR6","MP2K1","MK01","ELF3","PLD1","KPCD1","ANXA1","EZRI","DCN","HBEGF","CAV2","CLH2","HXC10","LYN","S10A7","S10A9","SEC13","AT1A1","ICAM1","CD166","CMTM8","SOCS5","P3C2A","SH3BGRL3","ADRBK1","PDC6I","NHRF1","ESR1","ITA5","GNAI2","KU70","MUC1","P85A","KPCD","TNC","COX2","PIPNA","EPS8","BTC","AMH","CAV3","SNX1","GRB7","GRB14","EREG","KCC2G","RGS16","ZPR1","PKIA","SCAM1","SCAMP3","KCC2A","PTN2","ERRFI1","DEGS1","AGTR1","CYH2","PTN12","HDAC6","GELS","HSP71","ARPC5","SH3L1","UB2V2","SPCS2","DYHC1","M3K12","PRDX1","AKAP12","FKBP4","CNTN2","GBRL2","TPIS","OAZ1","MATR3","MDHC","HSP74","STIP1","PUR9","EF1G","COF1","RAB3A","AATC","FAAA","ENPL","UROD","ALDOA","AHSA1","SRBS2","HEXI1","Q59EJ3","PPP5","ADRM1","PDK3","FKBP8","CO9A3","DCTN2","GPM6B","AHNAK","ARF6","PSMD4","SEPP1","PTGDS","PRKAR1B","CBR1","CSRP1","UCHL1","PSA7","MEGF6","AIP","MOBL3","PIN4","PLPL2","DNJC4","NECD","CD049","RUSC2","KCTD9","HS90A","Q59G22","Q53GZ6","SNP25","Q4QQI8","PLCH1","MAST1","STUB1","NOE1","AP2M1","PRCC","Q8N4S1","CEND","SH321","SYN1","Q59GR8","TMCO3","LST2","ATPG","SMD2","GRB2","YES","SHC1","STAT3","CDC2","FAK1","EPHA2","FAK2","BCAR1","BCAR3","ACK1","ERBB3","ERBB4","CD45","PTPRJ","STA5A","PTN11","CSK","KPCA","SHIP2","VAV","VAV2","STAT1","SRC","JAK2","PTK6","Q6PID4","NCK1","IGF1R","KAPCA","NCK2","PTN6","RASA1","GRB10","PLCG1","SH3K1","SOCS1","STA5B","RIPK1","M3K14","PGFRB","SHC3","FER","PLEC1","IQGA1","UBS3A","PTPRB","VAV3","SH23A","SOCS3","H31","CD59","S10AB","1433S","FHL2","1433T","CTND1","1433E","HSPB1","CTNA1","AT2A2","CBLB","EFNB1","IMB1","K2C1","HNRPF","H2B1C","PAXI","THIO","K2C6A","HNRPK","DNJA3","HSP7C","HNRH1","H14","K2C5","GRP78","K1C17","GRP75","PKP2","HS90B","EF1A1","DNJA1","CH60","CAV1","K1C14","EPS15","GAB2","CALM","G3P","4F2","CTNB1","PLAK","ANDR","LRSM1","CAMLG","ZO1","CBLC","CADH1","HD","DOK2","SOS1","SOS2","MPIP1","KAP0","PLD2","q02297-7","CD44","CD82","SNX6","SNX2","SNX4","CEBPB","GAB1","UBIQ","HGS","H31T","K2C7","K2C8","K1C18","EPIPL","UBC","NEDD8","SH3G2","LRIG1","PEX19","EI2BG","TBRG4","AT1B1","CLCA","SC5A1","AREG","AP2A1","VAPA","FINC","LEG3","TCTP","CCD50","UBP8","PCNA","EPN1","SPY2","NDFIP1","NDFIP2","LEG1","NEDD4","ACTA","SMUF2","OTU7B","TBA1A","ADT2","SGSM2","RABX5","TBA4A","WASP")
EGFR_genes<-sort(unique(c(wiki,EGFR_genes2)))
EGFR_genes<-sort(unique(c(wiki)))
EGFR_genes<-sort(unique(c(wiki,EGFR_genes2,HitPredict)))
length(EGFR_genes)
# EGFR_genes=unique(EGFR_genes2)

f.input2(EGFR_genes,IDs[,2])
EGFR_genes.IDs<-IDs[IDs[,2]==EGFR_genes[1],]
for(i in 2:length(EGFR_genes)) EGFR_genes.IDs<-rbind(EGFR_genes.IDs,IDs[IDs[,2]==EGFR_genes[i],])
inter=f.input2(rownames(EGFR_genes.IDs),rownames(D.n$E))$inter

plot(D.n$E[inter,][1,],ylim=c(0,20),type="l")
for(i in 2:length(inter)){
  lines(D.n$E[inter,][i,],ylim=c(0,20),type="l")
}
matplot(t(D.n$E[inter,c(1,3,5,4,2,6)]),type="l",lty=1)
re
pheatmap(D.n$E[inter,c(1,3,5,4,2,6)],kmeans_k = 5)

a0<-c(1,2,3,4,5,6);   b0<-c(1,3,5);       b1<-c(2,4,6);         c0<-c(1,3,4,5);      c1<-c(2,6);         d0<-c(1,3,4,5,6);    d1<-c(2);            
FC=(rowMeans(D.n$E[inter,c0])-rowMeans(D.n$E[inter,c1]))
inter2=names(FC)[abs(FC)>0.5]

print(table(m.TEEMU[inter]))
print(table(m[inter]))

print(table(m.TEEMU[inter2]))
print(table(m[inter2]))


print(table(m.TEEMU[Paper.qgenesIDs]))
print(table(m[Paper.qgenesIDs]))

f.input5(rownames(D.n_L.BIC.Post_a.85.TEEMU),
        rownames(D.n_L.BIC.Post_b.85.TEEMU),
        rownames(D.n_L.BIC.Post_c.85.TEEMU),
        rownames(D.n_L.BIC.Post_d.85.TEEMU),
inter)

f.input5(rownames(D.n_L.BIC.Post_a.85),
         rownames(D.n_L.BIC.Post_b.85),
         rownames(D.n_L.BIC.Post_c.85),
         rownames(D.n_L.BIC.Post_d.85),
         inter)

f.input5(rownames(D.n_L.BIC.Post_a.85.TEEMU),
         rownames(D.n_L.BIC.Post_b.85.TEEMU),
         rownames(D.n_L.BIC.Post_c.85.TEEMU),
         rownames(D.n_L.BIC.Post_d.85.TEEMU),
         inter2)

f.input5(rownames(D.n_L.BIC.Post_a.85),
         rownames(D.n_L.BIC.Post_b.85),
         rownames(D.n_L.BIC.Post_c.85),
         rownames(D.n_L.BIC.Post_d.85),
         inter2)
}