f.row<-function(a){
  rownames(a)<-c(as.character(FILE[,1]))
  a
}

T120 <- function(x) {1.2*x}     
T80 <- function(x) {0.8*x}
T125 <- function(x) {1.25*x} #== {(1/0.8)*x}     
T83 <- function(x) {(1/1.2)*x}
logF_B.A <- function(x) {log(x[2]/x[1],base=2)}
FC_B.A <- function(x) {(x[2]/x[1])}
logF_A.B <- function(x) {log(x[1]/x[2],base=2)}
f.dist_a.b<-function(X,a,b){
  S<-X[,c(a,b)]
  
  qq<-c()                                           #rows (MARGIN=1), columns (MARGIN=2), or both (MARGIN=c(1,2));
  qq<-apply(S,MARGIN=1,FUN=dist)  
  qq<-cbind(qq,apply(S,MARGIN=1,FUN=mean)  )
  qq<-cbind(qq,apply(t(S[,1]),MARGIN=1,FUN=T80))
  qq<-cbind(qq,apply(t(S[,1]),MARGIN=1,FUN=T120))
  qq<-cbind(qq,apply(t(S[,2]),MARGIN=1,FUN=T80))
  qq<-cbind(qq,apply(t(S[,2]),MARGIN=1,FUN=T120))
  qq<-cbind(qq,apply(S,MARGIN=1,FUN=logF_B.A))
  qq<-cbind(qq,apply(S,MARGIN=1,FUN=logF_A.B))
  colnames(qq)<-c("eucl","mean","V1*80","V1*120","V2*80","V2*120","log(+EGF/-EGF)","log(-EGF/+EGF)")
  qq<-f.row(qq)
  qq
}

f.vgl_a.b<-function(X,a,b){
  S<-X[,c(a,b)]
  
  qq<-c()                                           #rows (MARGIN=1), columns (MARGIN=2), or both (MARGIN=c(1,2));
  qq<-apply(t(S[,1]),MARGIN=1,FUN=T80)
  qq<-cbind(qq,apply(t(S[,1]),MARGIN=1,FUN=T120))
  qq<-cbind(qq,apply(t(S[,1]),MARGIN=1,FUN=T83))
  qq<-cbind(qq,apply(t(S[,1]),MARGIN=1,FUN=T125))
  qq<-cbind(qq,apply(S,MARGIN=1,FUN=logF_A.B))
  colnames(qq)<-c("-EGF*80","-EGF*120","-EGF*83","-EGF*125","log(-EGF/+EGF)")
  #qq<-f.row(qq)
  qq
}

f.FOLD_Henri<-function(x) {x[2]/x[1]}
#f.FOLD_Henri<-function(x) {x[1]/x[2]}
f.vgl_henri<-function(X,a,b,Change,LOW,UP){
 #Probe+EGF/ Probe-EGF = Faktor:            Keine Änderung wenn Faktor zw. 0,8-1,2 liegt
  S<-X[,c(a,b)]
  qq<-apply(S,MARGIN=1,FUN=f.FOLD_Henri)
  DA<-cbind(qq,LOW,UP)

  if(Change==0){
    print("between == 1 ")
    w<-apply(DA,MARGIN=1,FUN=is.between)
  }else{
    print("no.between == 1 ")
    w<-apply(DA,MARGIN=1,FUN=is.no.between)
  }
  
  w
}


f.noChange<-function(x){
  if (x[1]<=1){ return(1)  ## keine Änderung zwischen a und b == 1
                #if (log(x[1])<=0){ return(1) 
  }else{return(0)}             ## Änderung zwischen a und b == 0
}
### Comment für 1 gegen 2 if (x[1]<=2){ return(1) ca 10000 drin

f.Change<-function(x){
  if (x[1]>1){ return(1)  ## Änderung zwischen a und b == 1
  }else{return(0)}             ## keine Änderung zwischen a und b == 0
}

f.3<-function(x){
  if (sum(x)==3){ return(1)  ## Änderung zwischen a und b == 1
  } else{return(NA)} 
}

f.3_names <- function(x){
  x1<-apply(x,MARGIN=1,FUN=f.3)
  x2 <- x1[!is.na(x1)]
  print(length(names(x2)))
  return(names(x2))
}

f.1<-function(x){
  if (x==1){ return(1)  ## Änderung zwischen a und b == 1
  } else{return(NA)} 
}

f.1_names <- function(x){
  x1<-  apply(as.matrix(x),MARGIN=1,FUN=f.1)
  x2 <- x1[!is.na(x1)]
  print(length(names(x2)))
  return(names(x2))
}

f.input4 = function (p,q,r,s,name){
  input  <-list(A=p,B=q,C=r,D=s)
  names(input)<-name
  #print(input)
  input
}

f.input3 = function (p,q,r,name){
  input  <-list(A=p,B=q,C=r)
  names(input)<-name
  #print(input)
  input
}

f.input2 = function (p,q,name){
  input  <-list(A=p,B=q)
  names(input)<-name
  #print(input)
  input
}

#is.between <- function(x, a, b) {
#  (x - a)  *  (b - x) > 0
#}

is.between <- function(X) {
  #X[2] <= X[1] <= X[3]
  
  if( ((X[1] - X[2])  *  (X[3] - X[1])) >= 0){
    #print(c("IN",X[2],X[1],X[3]))
    return (1)
  }else {    
    #print(c("OUT",X[2],X[1],X[3]));
    return (0)}
}

is.no.between <- function(X) {
  #X[2] > X[1] || X[1] > X[3]
  if( ((X[1] - X[2])  *  (X[3] - X[1])) <= 0){
    return (1)
  }else {return (0)}
}

f.intervall<-function(A,D,a,b,COL,Change){ 
  
  #A == DATA
  #D == from f.dist_a.b
  #a == a.th column of A  
  #b == b.th column of A
  #COL == 1 first column of A between D[,5:6]; 
  #COL == 2 2ed column of A between D[,3:4];
  #Change == 0 no.Change
  #Change == 1 Change
  
  print(COL)
  
  S<-A[,c(a,b)]
  if(COL==1){
    DA<-cbind(S[,1],D[,5:6])
  }else if (COL==2){
    DA<-cbind(S[,2],D[,3:4])
  }else{exit }
  
  if(Change==0){
    print("between == 1 ")
    w<-apply(DA,MARGIN=1,FUN=is.between)
  }else{
    print("no.between == 1 ")
    w<-apply(DA,MARGIN=1,FUN=is.no.between)
  }
  return (w)
}

f.between<-function(A,D,Change){ 
  
  #A == DATA
  #D == from f.vgl_a.b
  #Change == "no.C" ->  no.Change
  #Change == C -> Change
  
  DA<-cbind(A,D)

  if(Change=="C"){
    #print("between == 1 ")
    w<-apply(DA,MARGIN=1,FUN=is.between)
  }else{
    #print("no.between == 1 ")
    w<-apply(DA,MARGIN=1,FUN=is.no.between)
  }
  return (w)
}

f.logF.Change<-function(x){
  if (abs(x)>0.2){ return(1)  
  }else{return(0)}             
}

f.logF.no.Change<-function(x){
  if (abs(x)<0.2){ return(1)  
  }else{return(0)}             
}

f.logF.Change_a<-function(x,PARA_FOL){
  #print("|logF|>T -> Genes with change ")
  if (abs(x)>PARA_FOL){ return(1)  
  }else{return(0)}             
}
#CALL
#V14a<-apply(as.matrix(D[,7]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=0.1)
# parameter muss selben namen haben wie in Funktion
f.logF.no.Change_a<-function(x,PARA_FOL){
  #print("|logF|<T -> Genes with no change ")
  if (abs(x)<=PARA_FOL){ return(1)  
  }else{return(0)}             
}

f.ven <- function(G,g.Label,MAIN){
  colnames(G)<-g.Label
  c<-vennCounts(G)
  vennDiagram(c,main=MAIN)
  g.e.nam<-f.3_names(G)
  print(g.e.nam[1])
  return (g.e.nam)
}

f.qq_NV_plo<-function(LogFoldChange,STR){
  print(STR)
  QQSample<-sort(LogFoldChange)
  M<-length(QQSample)
  m1<-mean(QQSample)
  s1<-sd(QQSample)
  plot(QQSample,1:M/M, xlab="log2(fold-change)",ylab="cumulative distribution function", xlim=c(-5,5), main=paste(sep="",STR), type="l", pch=20)
  lines(QQSample, pnorm(QQSample, mean=m1, sd=s1), lwd=1, col="red")
  legend("bottomright", c("log2(fold-change)", "Gaussian distribution"), col=c("black", "red"), pch=20)
  
  plot(QQSample,1:M/M, xlab="log2(fold-change)",ylab="cumulative distribution function", xlim=c(-5,5), main=paste(sep="",STR), type="p", pch=20)
  lines(QQSample, pnorm(QQSample, mean=m1, sd=s1), lwd=3, col="red")
  legend("bottomright", c("log2(fold-change)", "Gaussian distribution"), col=c("black", "red"), pch=20)
  
  qqnorm(QQSample, main=paste(sep="","Normal Q-Q Plot \n",STR), ylim=c(-5,5), xlim=c(-5,5))
  qqline(QQSample,col="red")
}

########### remove rows 
#=< --le   || >=	--ge
f.compare_ge<-function(A,T){
  if(A>=T){return(1)}
  else{return(0)}
}
f.compare_le<-function(A,T){
  if(A<=T){return(1)}
  else{return(0)}
}
f.compare_eq<-function(A,T){
  if(A == T){return(1)}
  else{return(0)}
}

f.remove_pval_row<-function(A,T){
  
  B<-(as.matrix(A))
  q<-apply(B,MARGIN=1,FUN=f.compare_ge,T=T)
  
  if(sum(q)==0){return(TRUE)}
  else{FALSE}
}

f.get_remove_pval_row_matrix<-function(A,Names,T){
  
  #A is NxM a matrix of pval N=row=GenID M=col=sample/tread   
  #Names are the subset rownames ( default Names=rownames(A)) 
  # T = threshold
  TRUE_FALSE_VAL<-apply(A[Names,],MARGIN=1,FUN=f.remove_pval_row,T=T)
  m<-A[Names,]
  
  return(m[TRUE_FALSE_VAL,])  
}
########### remove rows 

########### return 2 set names 
f.ex<-function(b,A){
  q<-apply(as.matrix(A),MARGIN=1,FUN=f.compare_eq,T=b)
  if(sum(q)==0){return(TRUE)}
  else{FALSE}
}

f.two_name_list_set <-function(a,b,set,show.venn=FALSE){
  # input muss be a vector with names
  #eg
  #a<-rownames(over.names.1_Remov)
  #b<-he
  if (show.venn){
    input=f.input2( a, b,c("1. argument","2. argument"))
    venn(input,simplify=F)
  }
  
  if(set==1){
    print("Only in 1. argument ")
    exa<-apply(as.matrix(a),MARGIN=1,FUN=f.ex,A=b)
    only.a<-a[exa]
    OUT<-only.a
  }else if(set==2){
    print("Only in 2. argument ")
    exb<-apply(as.matrix(b),MARGIN=1,FUN=f.ex,A=a)
    only.b<-b[exb]
    OUT<-only.b
  }else if(set==3){
    print("in both -  intersect ")
    both<-intersect(a,b)
    OUT<-both
  }else if(set==4){
    print("in all -  intersect ")
    all<-union(a,b) 
    OUT<-all
  }
  return(OUT)
}
########### return 2 set names


if(0){ # old version 
    ###################### HEATMAP with correlation ##############################
    #### Markus Boenn
    CorrPWsignif <- function(data){
      M <- matrix(NA, nrow=nrow(data), ncol=nrow(data), dimnames=list(rownames(data), rownames(data)))
      N <- nrow(data)
      for(n1 in 1:N){
        for(n2 in 1:N){
          M[n1, n2] <- cor.test(data[n1, ], data[n2, ])$p.value
        }
      }
      return(round(M,4))
    }
    getCorHeatColors <- function(data){
      COL <- heat.colors(length(seq(-1,1,0.01)))
      CC <- round(cor(t(data)), 2); alle <- as.vector((CC+1.01)*100); COL0 <- COL[ min(alle):max(alle) ]
    }
    CorrHeat <- function(data, main, cluster=FALSE, getp=FALSE,round=TRUE){
      COL0 <- getCorHeatColors(data)
      if(round==TRUE){CC <- round(cor(t(data)), 2)};
      if(round==FALSE){CC <- cor(t(data))}
      M <- matrix('', nrow=nrow(data), ncol=nrow(data))
      if(getp==TRUE){ M <- CorrPWsignif(data); print(M) }
      if(cluster==FALSE){heatmap.2(CC, Rowv=NA, Colv=NA, dendrogram=NULL, trace='none', main=main, col=COL0, margins=c(16,16), cellnote=M)}
      if(cluster==TRUE){heatmap.2(CC, trace='none', main=main, col=COL0, margins=c(16,16))}
    }
    ####################################################
}
###################### HEATMAP with correlation ##############################
#### Markus Boenn
CorrPWsignif <- function(data,method){
  M <- matrix(NA, nrow=nrow(data), ncol=nrow(data), dimnames=list(rownames(data), rownames(data)))
  N <- nrow(data)
  for(n1 in 1:N){
    for(n2 in 1:N){
      M[n1, n2] <- cor.test(data[n1, ], data[n2, ],method=method)$p.value
    }
  }
  return(round(M,4))
}
getCorHeatColors <- function(data,method){
  COL <- heat.colors(length(seq(-1,1,0.01)))
  CC <- round(cor(t(data),method=method), 2); alle <- as.vector((CC+1.01)*100); COL0 <- COL[ min(alle):max(alle) ]
}
CorrHeat <- function(data, main, cluster=FALSE, getp=FALSE,round=TRUE,method="pearson"){
  COL0 <- getCorHeatColors(data,method)
  if(round==TRUE){CC <- round(cor(t(data),method=method), 2)};
  if(round==FALSE){CC <- cor(t(data),method=method)}
  M <- matrix('', nrow=nrow(data), ncol=nrow(data))
  if(getp==TRUE){ M <- CorrPWsignif(data); print(M) }
  if(cluster==FALSE){heatmap.2(CC, Rowv=NA, Colv=NA, dendrogram=NULL, trace='none', main=main, col=COL0, margins=c(16,16), cellnote=M)}
  if(cluster==TRUE){heatmap.2(CC, trace='none', main=main, col=COL0, margins=c(16,16))}
}
####################################################
