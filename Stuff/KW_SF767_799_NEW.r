library("limma")
library("VennDiagram")
library("venneuler")
library(hexbin)
library(ggplot2)
library(gplots)
#install.packages("ggdendro")
library(ggdendro)
#install.packages("rpart")
library(rpart)
library(hexbin)
library(calibrate)
library(MASS)

FOLDER="~/Promotion"
#FOLDER="~/SSHFS/SOL_home/weinhol/Promotion"

#no.bg <- read.delim("~/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_bkgd_SF767-1-799-6.txt",sep="\t")
#no.nobg <- read.delim("~/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_nobkgd_SF767-1-799-6.txt",sep="\t") #,encoding="latin1")
qu.bg <- read.delim( paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_quantile_bkgd_SF767-1-799-6.txt"),sep="\t")
#qu.nobg <- read.delim("~/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_quantile_nobkgd_SF767-1-799-6.txt",sep="\t",encoding="latin1")

#qu.nobg2 <- read.table("~/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_nobkgd_SF767-1-799-6_adj.txt",sep="\t",header=TRUE,fill=TRUE)
#qu.nobg2 <-read.csv("~/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_nobkgd_SF767-1-799-6.csv",sep="\t")
#qu.nobg2[1512,]
#qu.nobg["PROBE_ID"=="ILMN_1685079",]
#?read.delim
#SF2["ILMN_2101810",]
#FILE[1511,]

#FILE<-no.nobg
#FILE.N<-"nonorm_nobkgd"
#FILE<-qu.nobg
#FILE.N<-"qu_nobkgd"

FILE<-qu.bg
#FILE.N<-"qu_bkgd_corr10"
FILE.N<-"qu_bkgd"

henri <- read.delim("~/Promotion/Kappler_Wichmann_Medizin/Auswertung/Genes_Henri.txt",sep="\t")
path="~/Promotion/Kappler_Wichmann_Medizin/Auswertung/NEW_BG/"
#path="~/Promotion/Kappler_Wichmann_Medizin/Auswertung/"
source("~/Promotion/Kappler_Wichmann_Medizin/Auswertung/KW_functions.r")

SF1<-FILE[,3:6]        # SF.767.1  KO
SF2<-FILE[,7:10]       # SF.767.2  KO+EGF
SF3<-FILE[,11:14]      # SF.767.3  ALL
SF4<-FILE[,15:18]      # SF.767.4  ALL+EGF
SF5<-FILE[,19:22]      # SF.767.5  14
SF6<-FILE[,23:26]      # SF.767.6  14+EGF
rownames(FILE)<-as.character(FILE[,1])


na.1=c("SF767.1","SF767.2","SF767.3","SF767.4","SF767.5","SF767.6")

X1<-FILE[,27:30]      #X799.1  KO
X2<-FILE[,31:34]      #X799.2  KO+EGF
X3<-FILE[,35:38]      #X799.3  ALL
X4<-FILE[,39:42]      #X799.4  ALL+EGF
X5<-FILE[,43:46]      #X799.5  14
X6<-FILE[,47:50]      #X799.6  14+EGF

na.2=c("799.1","799.2","799.3","799.4","799.5","799.6")
na.3=c("SF_767","799")


SF1<-f.row(SF1);SF2<-f.row(SF2);SF3<-f.row(SF3);SF4<-f.row(SF4);SF5<-f.row(SF5);SF6<-f.row(SF6);
X1<-f.row(X1)  ;X2<-f.row(X2)  ;X3<-f.row(X3)  ;X4<-f.row(X4)  ;X5<-f.row(X5)  ;X6<-f.row(X6);
#qu.nobg<-f.row(SF1);

X<-cbind(X1[,1],X2[,1],X3[,1],X4[,1],X5[,1],X6[,1])
colnames(X)<-c(na.2)
X<-f.row(X)
X.pval<-cbind(X1[,2],X2[,2],X3[,2],X4[,2],X5[,2],X6[,2])
colnames(X.pval)<-c(na.2)
X.pval<-f.row(X.pval)

SF<-cbind(SF1[,1],SF2[,1],SF3[,1],SF4[,1],SF5[,1],SF6[,1])
colnames(SF)<-c(na.1)
SF<-f.row(SF)
SF.pval<-cbind(SF1[,2],SF2[,2],SF3[,2],SF4[,2],SF5[,2],SF6[,2])
colnames(SF.pval)<-c(na.1)
SF.pval<-f.row(SF.pval)

#new for bg Dataset
  corr=0 #10
  SF<-SF-corr
  pVal<-apply(SF,1,remove.neg.val)
  print(sum(pVal))
  SF<-SF[pVal,]
  SF.pval<-SF.pval[pVal,]

  pValX<-apply(X,1,remove.neg.val)
  print(sum(pValX))
  X<-X[pValX,]
  X.pval<-X.pval[pValX,]

##################################
############ Pval estimation   ### 
##################################
SF.pval.max<-apply(SF.pval,MARGIN=1,FUN=max )
f.anz_le_pval<-function(pval,LIST){
  length(SF.pval[LIST<pval,1])
}
f.anz_ge_pval<-function(pval,LIST){
  length(SF.pval[LIST>pval,1])
}
#pvals<-c(1/1000000,1/100000,1/10000,1/1000,1/800,1/600,1/400,1/200,1/100,1/10)
#pvals<-c(1/1000000,1/100000,1/10000,1/1000,1/500,1/100,1/10)

pvals<-seq(from=1/1000000,to=.10,by=.00005)
pdf(paste(sep="",path,FILE.N,"_Pval_",na.3[1],".pdf"))
a<-apply(as.matrix(pvals),MARGIN=1,FUN=f.anz_le_pval,SF.pval.max)
plot(pvals,a,log="xy",type="b",ylim=c(11000,18000),ylab="number of genes", main="pval < threshold")
plot(pvals,a,type="b",ylim=c(11000,18000),ylab="number of genes", main="pval < threshold")

b<-apply(as.matrix(pvals),MARGIN=1,FUN=f.anz_ge_pval,SF.pval.max)
plot(pvals,b,log="xy",type="b",ylim=range(b),ylab="number of genes", main="pval > threshold")
plot(pvals,b,type="b",ylim=range(b),ylab="number of genes", main="pval > threshold")
dev.off()
##################################
############ Pval estimation   ### 
##################################



##################################
############ Heatmap      ######## 
##################################
#plot(X[,1],X[,2], log="",xlim=range(X[,1]),ylim=range(X[,2]),main=paste(sep="","a"))
#lines(lowess(X[,1],X[,2]), col="blue",lty=2) 
#abline(lm(X[,1]~X[,2]))
#abline(0,1)
#abline(h=0,v=0)

meth="euclidean";
#meth="maximum" 
#meth="manhattan"
pdf(paste(sep="",path,FILE.N,"_Heatmap_799_SF.767.pdf"))
  dists <- dist( t( X ) , method = meth)
  heatmap( as.matrix( dists ),symm=TRUE, scale="none", margins=c(6,6),main=(paste(sep="",na.3[2] )));    
  
  dists <- dist( t( SF ) , method = meth)
  heatmap( as.matrix( dists ),symm=TRUE, scale="none", margins=c(6,6),main=(paste(sep="",na.3[1])));  	
  
  a<-cbind(X,SF)
  dists <- dist( t( a ) , method = meth)
  heatmap( as.matrix( dists ),symm=TRUE, scale="none", margins=c(6,6),main=(paste(sep="",na.3[2]," and ",na.3[1] )));		
dev.off()

pdf(paste(sep="",path,FILE.N,"_Heatmap.pavl_corr_799_SF.767.pdf"))
  tmp<-f.get_remove_pval_row_matrix(SF.pval,rownames(SF.pval),0.001)
  SF_pcor<-SF[rownames(tmp),]
  dists <- dist( t( SF_pcor ) , method = meth)
  heatmap( as.matrix( dists ),symm=TRUE, scale="none", margins=c(6,6),main=(paste(sep="",na.3[1])));    
  
  tmp<-f.get_remove_pval_row_matrix(X.pval,rownames(X.pval),0.001)
  X_pcor<-X[rownames(tmp),]
  dists <- dist( t( X_pcor ) , method = meth)
  heatmap( as.matrix( dists ),symm=TRUE, scale="none", margins=c(6,6),main=(paste(sep="",na.2[1])));
  
  XSF<-cbind(X,SF)
  XSF.pval<-cbind(X.pval,SF.pval)
  tmp<-f.get_remove_pval_row_matrix(XSF.pval,rownames(XSF.pval),0.001)
  XSF_pcor<-XSF[rownames(tmp),]
  dists <- dist( t( XSF_pcor ) , method = meth)
  heatmap( as.matrix( dists ),symm=TRUE, scale="none", margins=c(6,6),main=(paste(sep="",na.3[2]," and ",na.3[1] )));   
  
  heatmap.2(cor(SF_pcor),trace="none")
  heatmap.2(cor(X_pcor),trace="none")
  heatmap.2(cor(XSF_pcor),trace="none")
  
  CorrHeat(t(SF_pcor),na.3[1],cluster=T, getp=F,round=F)
  CorrHeat(t(X_pcor),na.2[1],cluster=T, getp=F,round=F)
  CorrHeat(t(XSF_pcor),"both",cluster=T, getp=F,round=F)

dev.off()
##################################
############ Heatmap      ######## 
##################################

##################################
############ Vergleich     ####### 
##################################

##################################
############  VGL        ######## 
##################################
main_VGL_FOLDC <-function (PARA_FOL,A,a1,a2,b1,b2,c1,c2){
  
  V.KO<-f.vgl_a.b(A,a1,a2)
  V.ALL<-f.vgl_a.b(A,b1,b2)
  V.14<-f.vgl_a.b(A,c1,c2)
  
  print(paste(sep="",a1," ",a2," --- |logF|>T -> Genes with change "))
  FC.KO<-apply(as.matrix(V.KO[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL)
  print(paste(sep="",b1," ",b2," --- |logF|<T -> Genes with no change "))
  FC.ALL<-apply(as.matrix(V.ALL[,5]),MARGIN=1,FUN=f.logF.no.Change_a,PARA_FOL=PARA_FOL)
  print(paste(sep="",c1," ",c2," --- |logF|>T -> Genes with change "))
  FC.14<-apply(as.matrix(V.14[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL)
  
  FC.g<-cbind(FC.KO,FC.ALL,FC.14)
  g.nam<-f.3_names(FC.g)
  c<-vennCounts(FC.g)
  vennDiagram(c)
  return (g.nam)
}

main_VGL_FOLDC_2 <-function (PARA_FOL.C,PARA_FOL.noC,A,a1,a2,b1,b2,c1,c2){
  
  V.KO<-f.vgl_a.b(A,a1,a2)
  V.ALL<-f.vgl_a.b(A,b1,b2)
  V.14<-f.vgl_a.b(A,c1,c2)
  
  print(paste(sep="",a1," ",a2," --- |logF|>T -> Genes with change "))
  FC.KO<-apply(as.matrix(V.KO[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL.C)
  print(paste(sep="",b1," ",b2," --- |logF|<T -> Genes with no change "))
  FC.ALL<-apply(as.matrix(V.ALL[,5]),MARGIN=1,FUN=f.logF.no.Change_a,PARA_FOL=PARA_FOL.noC)
  print(paste(sep="",c1," ",c2," --- |logF|>T -> Genes with change "))
  FC.14<-apply(as.matrix(V.14[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL.C)
  
  FC.g<-cbind(FC.KO,FC.ALL,FC.14)
  g.nam<-f.3_names(FC.g)
  c<-vennCounts(FC.g)
  vennDiagram(c)
  return (g.nam)
}

main_VGL_FOLDC_3 <-function (PARA_FOL.a,PARA_FOL.b,PARA_FOL.c,A,a1,a2,b1,b2,c1,c2){
  
  V.KO<-f.vgl_a.b(A,a1,a2)
  V.ALL<-f.vgl_a.b(A,b1,b2)
  V.14<-f.vgl_a.b(A,c1,c2)
  
  print(paste(sep="",a1," ",a2," --- |logF|>T -> Genes with change "))
  FC.KO<-apply(as.matrix(V.KO[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL.a)
  print(paste(sep="",b1," ",b2," --- |logF|<T -> Genes with no change "))
  FC.ALL<-apply(as.matrix(V.ALL[,5]),MARGIN=1,FUN=f.logF.no.Change_a,PARA_FOL=PARA_FOL.b)
  print(paste(sep="",c1," ",c2," --- |logF|>T -> Genes with change "))
  FC.14<-apply(as.matrix(V.14[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL.c)
  
  FC.g<-cbind(FC.KO,FC.ALL,FC.14)
  g.nam<-f.3_names(FC.g)
  c<-vennCounts(FC.g)
  vennDiagram(c)
  return (g.nam)
}

hist(log(SF_pcor),freq=F)
lines(density(log(SF_pcor)))
m=mean(log(SF_pcor))
s=sd(log(SF_pcor))
x=seq(0,10,length=20000)
y=rnorm(x,mean=m,sd=s[1])
hist(y,freq=F)
lines(density(y))
rug(log(SF_pcor))


##                          Change and no.change
MEGF.1<-main_VGL_FOLDC_2(log2(1.5),log2(1.5),SF_pcor,1,2,3,4,5,6)
MEGF.2<-main_VGL_FOLDC_2(log2(1.5),log2(1.2),SF_pcor,1,2,3,4,5,6)
MEGF.12<-f.two_name_list_set(MEGF.1, MEGF.2,3,TRUE)
##                        KO         , All       , 14 
MEGF.1<-main_VGL_FOLDC_3(log2(1.25),log2(1.25),log2(1.25),SF_pcor,1,2,3,4,5,6)
MEGF.1<-main_VGL_FOLDC_3(log2(1.5),log2(1.5),log2(1.5),SF_pcor,1,2,3,4,5,6)
MEGF.2<-main_VGL_FOLDC_3(log2(1.5),log2(1.2),log2(1.5),SF_pcor,1,2,3,4,5,6)
MEGF.12<-f.two_name_list_set(MEGF.1, MEGF.2,3,TRUE)



f.ven_over4_b<-function(e1,b,FOLDC){
  input=f.input4( b[[e1[1]]], b[[e1[2]]],b[[e1[3]]],b[[e1[4]]],c(FOLDC[e1[1]],FOLDC[e1[2]],FOLDC[e1[3]],FOLDC[e1[4]]))
  h1<-venn(input,simplify=T)
}
f.get_log2F_B.A<-function(SF){
  
  qq<-c()                                           #rows (MARGIN=1), columns (MARGIN=2), or both (MARGIN=c(1,2));
  qq<-apply(SF[,c(1,2)],MARGIN=1,FUN=logF_B.A)
  qq<-cbind(qq,apply(SF[,c(3,4)],MARGIN=1,FUN=logF_B.A))
  qq<-cbind(qq,apply(SF[,c(5,6)],MARGIN=1,FUN=logF_B.A))
  
  colnames(qq)<-c("KO_+/-","ALL_+/-","14_+/-")
  return(qq)
}
f.get_FC_B.A<-function(SF){
  
  qq<-c()                                           #rows (MARGIN=1), columns (MARGIN=2), or both (MARGIN=c(1,2));
  qq<-apply(SF[,c(1,2)],MARGIN=1,FUN=FC_B.A)
  qq<-cbind(qq,apply(SF[,c(3,4)],MARGIN=1,FUN=FC_B.A))
  qq<-cbind(qq,apply(SF[,c(5,6)],MARGIN=1,FUN=FC_B.A))
  
  colnames(qq)<-c("KO_+/-","ALL_+/-","14_+/-")
  return(qq)
}
f.list.to.txt<-function(SF,genN_1.25,log2ed=TRUE,path="",FILE.N="",pval=""){

genes_1.25<-SF[genN_1.25,]
genes_1.25.LFC<-f.get_log2F_B.A(genes_1.25)
if(log2ed==TRUE){genes_1.25<-f.get_log2F_B.A(genes_1.25)
                 Str="log2.FoldChange"}
if(log2ed==FALSE){genes_1.25<-f.get_FC_B.A(genes_1.25)
                  Str="FoldChange"}
genes_1.25<-cbind(as.character(FILE[genN_1.25,"PROBE_ID"]),genes_1.25)
genes_1.25<-cbind(genes_1.25,as.character(FILE[genN_1.25,"SYMBOL"]))
genes_1.25<-cbind(genes_1.25,as.character(FILE[genN_1.25,"SEARCH_KEY"]))
genes_1.25<-cbind(genes_1.25,as.character(FILE[genN_1.25,"ILMN_GENE"]))
genes_1.25<-cbind(genes_1.25,as.character(FILE[genN_1.25,"CHROMOSOME"]))
genes_1.25<-cbind(genes_1.25,as.character(FILE[genN_1.25,"DEFINITION"]))
genes_1.25<-cbind(genes_1.25,as.character(FILE[genN_1.25,"SYNONYMS"]))

genes_1.25<-genes_1.25[order(abs(as.double(genes_1.25.LFC[,"14_+/-"])),decreasing=T),]
genes_1.25[,2]<-signif(as.double(genes_1.25[,2]),3)
genes_1.25[,3]<-signif(as.double(genes_1.25[,3]),3)
genes_1.25[,4]<-signif(as.double(genes_1.25[,4]),3)

if(log2ed==TRUE){colnames(genes_1.25)<-c("PROBE_ID","KO_log2(EGF+/EGF-)","ALL_log2(EGF+/EGF-)","14_log2(EGF+/EGF-)","SYMBOL","SEARCH_KEY","ILMN_GENE","CHROMOSOME","DEFINITION","SYNONYMS")}
if(log2ed==FALSE){colnames(genes_1.25)<-c("PROBE_ID","KO_(EGF+/EGF-)","ALL_(EGF+/EGF-)","14_(EGF+/EGF-)","SYMBOL","SEARCH_KEY","ILMN_GENE","CHROMOSOME","DEFINITION","SYNONYMS")}

#2^log2(1.25) == 1.25

write.table(genes_1.25,file=paste( path,FILE.N,'.',Str,'.T=',pval,'.txt', sep='' ),sep="\t",append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE);
S="#open with UFT-8 and language = English"
write.table(S,file=paste( path,FILE.N,'.',Str,'.T=',pval,'.csv', sep='' ),sep="\t",append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE);
write.table(genes_1.25,file=paste( path,FILE.N,'.',Str,'.T=',pval,'.csv', sep='' ),sep="\t",append=TRUE,quote=FALSE,col.names=TRUE,row.names=FALSE);

return(genes_1.25)
}

pdf(paste(sep="",path,FILE.N,"_FoldChange_[CG]_SF.767.pdf"),width=11,height=11)
#FOLDC<-seq(from=1.1,to=2.6,by=.1)
FOLDC<-seq(from=1,to=2.25,by=.05)

b<-apply(as.matrix(log2(FOLDC)),MARGIN=1,FUN=main_VGL_FOLDC,SF_pcor,1,2,3,4,5,6)
b2<-sapply(b,length)
b3<-t(as.matrix(b2))
colnames(b3)<-FOLDC
t(b3)
plot(FOLDC,b2,type="b",ylim=range(b2),xlim=c(1,2.5),xlab="foldChange",ylab="number regulated genes [CG]", main="log2(foldChange)")
plot(log2(FOLDC),b2,type="b",ylim=range(b2),xlab="log2(foldChange)",ylab="number regulated genes [CG]", main="log2(foldChange)")
barplot(b2,FOLDC,names.arg=(FOLDC))
barplot(b3,xlab="threshold",ylab=" number of Genes")

f.ven_over4_b(c(1,2,3,4),b,FOLDC)
f.ven_over4_b(c(3,4,5,6),b,FOLDC)
f.ven_over4_b(c(5,6,7,8),b,FOLDC)
f.ven_over4_b(c(10,12,14,16),b,FOLDC)
f.ven_over4_b(c(16,18,22,24),b,FOLDC)

e1<-c(3,6)
he.only.36<-f.two_name_list_set(b[[e1[1]]], b[[e1[2]]],3,TRUE)
e1<-c(6,10)
he.only.610<-f.two_name_list_set(b[[e1[1]]], b[[e1[2]]],3,TRUE)

e1<-c(6,11)
he.only.611<-f.two_name_list_set(b[[e1[1]]], b[[e1[2]]],3,TRUE)

dev.off()

l1<-f.list.to.txt(SF,b[[6]],log2ed=TRUE,path,FILE.N,FOLDC[6])
l2<-f.list.to.txt(SF,b[[6]],log2ed=FALSE,path,FILE.N,FOLDC[6])

l3<-f.list.to.txt(SF,b[[10]],log2ed=TRUE,path,FILE.N,FOLDC[10])
l4<-f.list.to.txt(SF,b[[10]],log2ed=FALSE,path,FILE.N,FOLDC[10])

l5<-f.list.to.txt(SF,b[[11]],log2ed=TRUE,path,FILE.N,FOLDC[11])
l6<-f.list.to.txt(SF,b[[11]],log2ed=FALSE,path,FILE.N,FOLDC[11])

l7<-f.list.to.txt(SF,he.only.36,log2ed=TRUE,path,FILE.N,paste(seq="",FOLDC[3],"&",FOLDC[6]))
l8<-f.list.to.txt(SF,he.only.36,log2ed=FALSE,path,FILE.N,paste(seq="",FOLDC[3],"&",FOLDC[6]))

l9<-f.list.to.txt(SF,he.only.610,log2ed=TRUE,path,FILE.N,paste(seq="",FOLDC[6],"&",FOLDC[10]))
l10<-f.list.to.txt(SF,he.only.610,log2ed=FALSE,path,FILE.N,paste(seq="",FOLDC[6],"&",FOLDC[10]))

l9<-f.list.to.txt(SF,he.only.611,log2ed=TRUE,path,FILE.N,paste(seq="",FOLDC[6],"&",FOLDC[11]))
l10<-f.list.to.txt(SF,he.only.611,log2ed=FALSE,path,FILE.N,paste(seq="",FOLDC[6],"&",FOLDC[11]))

MEGF.1<-main_VGL_FOLDC_2(log2(FOLDC[6]),log2(FOLDC[6]),SF_pcor,1,2,3,4,5,6)
MEGF.2<-main_VGL_FOLDC_2(log2(FOLDC[11]),log2(FOLDC[11]),SF_pcor,1,2,3,4,5,6)
pdf(paste(sep="",path,FILE.N,"_SF.767_FoldChange_[CG]_Henri.pdf"),width=11,height=11)
barplot(b3,xlab="threshold",ylab=" number of Genes")
input<-f.input2(MEGF.1,MEGF.2,c(FOLDC[6],FOLDC[11]))
venn(input,simplify=F)
dev.off()


main_VGL_E <-function (PARA_FOL,A,a1,a2,b1,b2,c1,c2){

MEGF.1<-main_VGL_E(log2(1.25),SF_pcor,1,2,3,4,5,6)

pdf(paste(sep="",path,FILE.N,"_FoldChange_[E]_SF.767.pdf"),width=9,height=9)
FOLDC<-seq(from=1.1,to=2.6,by=.1)
b<-apply(as.matrix(log2(FOLDC)),MARGIN=1,FUN=main_VGL_E,SF_pcor,1,2,3,4,5,6)
plot(FOLDC,b,type="b",ylim=range(b),xlim=c(1.1,2.5),xlab="foldChange",ylab="number regulated genes [CG]", main="log2(foldChange)")
plot(log2(FOLDC),b,type="b",ylim=range(b),xlab="log2(foldChange)",ylab="number regulated genes [CG]", main="log2(foldChange)")
dev.off()


main_VGL <-function (A,a1,a2,b1,b2,c1,c2,PARA_FOL){
  
  V.KO<-f.vgl_a.b(A,a1,a2)
  V.ALL<-f.vgl_a.b(A,b1,b2)
  V.14<-f.vgl_a.b(A,c1,c2)
  
  STR=paste(as.character(a1),as.character(a2),sep="")
  f.qq_NV_plo(V.KO[,5],STR)
  
  STR=paste(as.character(b1),as.character(b2),sep="")
  f.qq_NV_plo(V.ALL[,5],STR)
  
  STR=paste(as.character(c1),as.character(c2),sep="")
  f.qq_NV_plo(V.14[,5],STR)
  
  print(paste(sep="",a1," ",a2," --- |logF|>T -> Genes with change "))
  FC.KO<-apply(as.matrix(V.KO[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL)
  print(paste(sep="",b1," ",b2," --- |logF|<T -> Genes with no change "))
  FC.ALL<-apply(as.matrix(V.ALL[,5]),MARGIN=1,FUN=f.logF.no.Change_a,PARA_FOL=PARA_FOL)
  print(paste(sep="",c1," ",c2," --- |logF|>T -> Genes with change "))
  FC.14<-apply(as.matrix(V.14[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL)
  
  FC.KO.n<-f.1_names(FC.KO)
  FC.ALL.n<-f.1_names(FC.ALL)
  FC.14.n<-f.1_names(FC.14)
  
  FC.g<-cbind(FC.KO,FC.ALL,FC.14)
  g.nam<-f.3_names(FC.g)
  c<-vennCounts(FC.g)
  vennDiagram(c)
  return (cbind(FC.g))
  
  #FC.ALL.C<-apply(as.matrix(V.ALL[,5]),MARGIN=1,FUN=f.logF.Change_a,PARA_FOL=PARA_FOL)
  #FC.g.C<-cbind(FC.ALL.C,FC.14)
  #c<-vennCounts(FC.g.C)
  #vennDiagram(c)
  
  #FC.g<-cbind(FC.KO,FC.ALL.C,FC.14)
  #g.nam<-f.3_names(FC.g)
  #c<-vennCounts(FC.g)
  #vennDiagram(c)
}

f.plot_venn<-function(SF,What,FILE.N,path,henri){ 
  ##### VENN PLOT
  PARA_FOL=log2(1.25)
  PARA_FOL2=log2(1.2)
  pdf(paste(sep="",path,FILE.N,"QQ_vgl_",What,".767_v1.2.pdf"))
  KO.1<-main_VGL(SF,1,2,1,4,1,6,PARA_FOL)
  KO.2<-main_VGL(SF,1,2,1,4,1,6,PARA_FOL2)
  
  MEGF.1<-main_VGL(SF,1,2,3,4,5,6,PARA_FOL)
  MEGF.2<-main_VGL(SF,1,2,3,4,5,6,PARA_FOL2)
  dev.off()
   
  pdf(paste(sep="",path,FILE.N,"_Venn_vgl_",What,".767_v1.2.pdf"))
  
  g.Label<-c("12 up","14 eq","16 up")
  over.names.1<-f.ven(KO.1,g.Label,"log2(-EFG/+EGF) with threshold log2(1.25) -> KO <1>")
  over.names.2<-f.ven(KO.2,g.Label,"log2(-EFG/+EGF) with threshold log2(1.2) -> KO <1>")
  
  g.Label<-c("12 up","34 eq","56 up")
  over.names.3<-f.ven(MEGF.1,g.Label,"log2(-EFG/+EGF) with threshold log2(1.25) -> -EFG <1,3,5>")
  over.names.4<-f.ven(MEGF.2,g.Label,"log2(-EFG/+EGF) with threshold log2(1.2) -> -EFG <1,3,5>")
  
  input=f.input2( over.names.1, over.names.2,c("KO <1> T=log2(1.25)","log fold change\nKO <1> T=log2(1.2)"))
  venn(input,simplify=F)
  input=f.input2( over.names.3, over.names.4,c("-EFG <1,3,5> T=log2(1.25)","log fold change\n-EFG <1,3,5> T=log2(1.2)"))
  venn(input,simplify=F)
  
  input=f.input4( over.names.1, over.names.2,over.names.3, over.names.4,c("KO <1> \nT=log2(1.25)","KO <1> T=log2(1.2)","-EFG <1,3,5>\n T=log2(1.25)","log fold change\n-EFG <1,3,5> T=log2(1.2)"))
  #venn(input,simplify=F)
  venn(input,simplify=T)
  
  he<-c(as.character(as.matrix(henri)))
  input=f.input2( he, over.names.1,c("H","log fold change\nKO <1> T=log2(1.25)"))
  h1<-venn(input,simplify=F)
  input=f.input2( he, over.names.2,c("Henri","log fold change\nKO <1> T=log2(1.2)"))
  h2<-venn(input,simplify=F)
  input=f.input2( he, over.names.3,c("H","log fold change\n-EFG <1,3,5> T=log2(1.25)"))
  h3<-venn(input,simplify=F)
  input=f.input2( he, over.names.4,c("H","log fold change\n-EFG <1,3,5> T=log2(1.2)"))
  h4<-venn(input,simplify=F)
  
  input=f.input3( he,  over.names.1, over.names.2,c("Henri","KO <1> T=log2(1.25)","overlap with Henri \nKO <1> T=log2(1.2)"))
  h3<-venn(input,simplify=F)
  
  over.names.2_Remov<-f.get_remove_pval_row_matrix(SF.pval,over.names.2,0.001)
  over.names.1_Remov<-f.get_remove_pval_row_matrix(SF.pval,over.names.1,0.001)
  he_Remov<-f.get_remove_pval_row_matrix(SF.pval,he,0.001)
  print(c("Henri",length(he),length(rownames(he_Remov))))
  print(c("T=log2(1.25)",length(over.names.1),length(rownames(over.names.1_Remov))))
  print(c("T=log2(1.2)",length(over.names.2),length(rownames(over.names.2_Remov))))
  input=f.input2( rownames(he_Remov), rownames(over.names.2_Remov),c("Henri","log fold change - pval correct 0.001\nKO <1> T=log2(1.2)"))
  venn(input,simplify=F)
  input=f.input3( rownames(he_Remov), rownames(over.names.1_Remov),rownames(over.names.2_Remov),c("Henri","KO <1> T=log2(1.25)","overlap with Henri - pval correct 0.001 \nKO <1> T=log2(1.2)"))
  h3<-venn(input,simplify=F)
  
  dev.off()
  ##### VENN PLOT
}

f.plot_venn(SF,"SF",FILE.N,path,henri)
f.plot_venn(X,"X",FILE.N,path,henri)


##################################
############  VGL         ########
##################################

main_Henri <-function (A,a1,a2,b1,b2,c1,c2){
  
  V.KO<-f.vgl_henri(A,a1,a2,1,0.84,1.16)
  V.ALL<-f.vgl_henri(A,b1,b2,0,0.80,1.2)
  V.14<-f.vgl_henri(A,c1,c2,1,0.84,1.16)
  
  FC.KO.n<-f.1_names(V.KO)
  FC.ALL.n<-f.1_names(V.ALL)
  FC.14.n<-f.1_names(V.14)
  
  FC.g<-cbind(V.KO,V.ALL,V.14)
  g.nam<-f.3_names(FC.g)
  c<-vennCounts(FC.g)
  vennDiagram(c)
  
  return (cbind(FC.g))
}
KO.H<-main_Henri(SF,1,2,1,4,1,6)
over.names.1<-f.ven(KO.H,g.Label,"SF- 0.8 < +EFG/-EGF < 1.2")
input=f.input2( he, over.names.1,c("Henri","0.8 < +EFG/-EGF < 1.2\nKO <1> "))
h1<-venn(input,simplify=F)

he<-c(as.character(as.matrix(henri)))
What="My"
pdf(paste(sep="",path,FILE.N,"_Venn_compare_with_HENRI_",What,".767_v1.2.pdf"))

KO.H<-main_Henri(SF,1,2,1,4,1,6)
g.Label<-c("12 up","14 eq","16 up")
over.names.1<-f.ven(KO.H,g.Label,"SF- 0.8 < +EFG/-EGF < 1.2")
input=f.input2( he, over.names.1,c("Henri","0.8 < +EFG/-EGF < 1.2\nKO <1> "))
h1<-venn(input,simplify=F)

over.names.1_Remov<-f.get_remove_pval_row_matrix(SF.pval,over.names.1,0.1)
print(c("MY",length(over.names.1),length(rownames(over.names.1_Remov))))
tmp_KO<-rownames(over.names.1_Remov)
input=f.input2( he, rownames(over.names.1_Remov),c("Henri","SF 0.8 <= +EFG/-EGF <= 1.2 - pval correct 0.1\nKO <1>"))
h1<-venn(input,simplify=F)

over.names.1_Remov<-f.get_remove_pval_row_matrix(SF.pval,over.names.1,0.001)
print(c("MY",length(over.names.1),length(rownames(over.names.1_Remov))))
he_Remov<-f.get_remove_pval_row_matrix(SF.pval,he,0.001)
print(c("Henri",length(he),length(rownames(he_Remov))))

input=f.input2( rownames(he_Remov), rownames(over.names.1_Remov),c("Henri","SF 0.8 <= +EFG/-EGF <= 1.2 - pval correct 0.001\nKO <1>"))
venn(input,simplify=F)

##### return names lists !!!!! don't remove!!!
#he.only<-f.two_name_list_set(rownames(over.names.1_Remov),he,2)
#SF[he.only,]
#he.only.2<-f.two_name_list_set(rownames(over.names.1_Remov), rownames(he_Remov),2)
#SF[he.only.2,]
##### return names lists !!!!! don't remove!!!


KO.H<-main_Henri(SF,1,2,3,4,5,6)
g.Label<-c("12 up","34 eq","56 up")
over.names.1<-f.ven(KO.H,g.Label,"SF - 0.8 < +EFG/-EGF < 1.2")
input=f.input2( he, over.names.1,c("Henri","0.8 < +EFG/-EGF < 1.2\n-EFG <1,3,5> "))
h1<-venn(input,simplify=F)

over.names.1_Remov<-f.get_remove_pval_row_matrix(SF.pval,over.names.1,0.1)
print(c("MY",length(over.names.1),length(rownames(over.names.1_Remov))))
tmp_135<-rownames(over.names.1_Remov)
input=f.input2( he, rownames(over.names.1_Remov),c("Henri","SF 0.8 <= +EFG/-EGF <= 1.2 - pval correct 0.1\n -EFG <1,3,5>"))
h1<-venn(input,simplify=F)

over.names.1_Remov<-f.get_remove_pval_row_matrix(SF.pval,over.names.1,0.001)
print(c("MY",length(over.names.1),length(rownames(over.names.1_Remov))))
he_Remov<-f.get_remove_pval_row_matrix(SF.pval,he,0.001)
print(c("Henri",length(he),length(rownames(he_Remov))))

input=f.input2( rownames(he_Remov), rownames(over.names.1_Remov),c("Henri","SF 0.8 <= +EFG/-EGF <= 1.2 - pval correct 0.001\n-EFG <1,3,5>"))
venn(input,simplify=F)

input=f.input2( tmp_KO, tmp_135,c("KO <1>","SF 0.8 <= +EFG/-EGF <= 1.2 - pval correct 0.1\n-EFG <1,3,5>"))
venn(input,simplify=F)

input=f.input3(he,tmp_KO,tmp_135,c("Henri","KO <1>","overlap with Henri - pval correct 0.1 \n-EFG <1,3,5>"))
h3<-venn(input,simplify=F)


dev.off()

#he.only.2<-f.two_name_list_set(rownames(over.names.1_Remov), tmp_KO,2,TRUE)

