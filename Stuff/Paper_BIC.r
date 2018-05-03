on='/Volumes/sol.informatik.uni-halle.de'

FOLDER="~/Promotion"
#FOLDER="~/SSHFS/SOL_home/weinhol/Promotion"
FOLDER=paste0(on,"/Promotion")

source(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/KW_functions.r"))
source(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/MY_DE_Bayes.functions.r"))




qu.bg <- read.delim( paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_quantile_bkgd_SF767-1-799-6.txt"),sep="\t")
FILE<-qu.bg
FILE.N<-"qu_bkgd"
path="~/Promotion/Kappler_Wichmann_Medizin/Auswertung/NEW_BG/"
IDs<-FILE[,c(1:2,54,55)]
rownames(IDs)<-as.character(FILE[,1])
head(IDs)
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

## pVal Filter
S<-0.01
goodPV<-(rowSums(cbind(SF.pval[,1]<S,SF.pval[,2]<S,SF.pval[,3]<S,SF.pval[,4]<S,SF.pval[,5]<S,SF.pval[,6]<S))==6)
RNsPV<-rownames(SF.pval)
length(RNsPV[goodPV])
RNs.sig<-intersect(rownames(SF),RNsPV[goodPV])
## pVal Filter

########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### 
############# use Limma to normalize !!!!
########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### 
on<-"~";#on<-"~/SSHFS/SOL_home/weinhol"

no.nobg <- read.ilmn(files=paste(sep="",on,"/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_nonorm_nobkgd_SF767-1-799-6.txt")  ,sep="\t",
                     ctrlfiles=paste(sep="",on,"/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/ControlProbeProfile_SF767-1-799-6.txt"),
                     other.columns=c("Detection","BEAD_STDERR","Avg_NBEADS")
                    ) #,encoding="latin1")
x<-no.nobg
##### 2. ### reduce to our 6 samples !!!! 
x2<-x
targets<-c()
targets$CellType=c("KO","KO_EGF","ALL","ALL_EGF","14","14_EGF")
x2$E<-x2$E[,c(1:6)]
x2$other$Detection<-x2$other$Detection[,c(1:6)]
x2$other$BEAD_STDERR<-x2$other$BEAD_STDERR[,c(1:6)]
x2$other$Avg_NBEADS<-x2$other$Avg_NBEADS[,c(1:6)]

pe2 <- propexpr(x2);dim(pe2) <- c(2,3);dimnames(pe2) <- list(CellType=c("Control","EGF"),Donor=c(1,2,3));pe2 ;mean(pe2[1:2,])
#The proportion of probes that are expressed varies from 49-54%. The average is 51,9 %.

##### 3.1 ### background correction followed by quantile normalization !!!! 
#neqc performs background correction followed by quantile normalization, using negative control probes for background correction and both negative and positive controls for normalization. 
#nec is similar but performs background correction only.
#After normalization, the intensities are log2 transformed and the control probes are removed.
table(x2$genes$Status)
y2 <- neqc(x2,negctrl="NEGATIVE") ; NORM<-"Limma_BG_QN"     #,robust=TRUE)
dim(y2)
expressed <- rowSums(y2$other$Detection < 0.05)  ==6 #>= 3 #
y2 <- y2[expressed,]
dim(y2)


########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### 
############# use Limma to normalize !!!!
########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### ########### 

#EGFR_genes<-c("ABL1","ABL2","ACTA1","ADRBK1","ADRM1","AGTR1","AHNAK","AHSA1","AIP","AKAP12","ALCAM","ALDOA","AMH","ANKS1A","ANKS1B","ANXA1","APBA3","APBB1","APBB2","APBB3","APPL1","AR","AREG","ARF4","ARF6","ATIC","ATP1A1","ATP1B1","ATP5C1","BLK","BMX","BTC","CALM1","CALM2","CALM3","CAMK2A","CAMK2G","CAMLG","CASP1","CAV1","CAV2","CAV3","CBL","CBLB","CBLC","CBR1","CD44","CD59","CD82","CDC25A","CDH1","CEACAM1","CEBPB","CEND1","CFL1","CISH","CLNK","CLTA","CLTCL1","CMTM8","CNTN2","COL9A3","COX2","CRK","CSRP1","CTNNB1","CTNND1","DCN","DCTN2","DEGS1","DNAJC4","DOK2","DOK4","DOK5","DOK6","DUSP3","DYNC1H1","EEF1G","EGF","ELF3","EPB41","EPHA2","EPPK1","EPS15","EPS8","ERBB2","ERBB3","ERBB4","EREG","ERRFI1","ESR1","EZR","FAH","FAS","FBXO25","FER","FES","FKBP4","FKBP8","FN1","GAB1","GABARAPL2","GAPDH","GNAI2","GOT1","GPM6B","GRB10","GRB14","GRB2","GRB7","HBEGF","HEXIM1","HGS","HIST3H3","HOXC10","HSP90B1","HSPA1A","HSPA4","HSPA8","HTT","ICAM1","INPPL1","IRS1","IRS4","ITGA5","JAK2","JUP","KCTD9P2","KRT17","KRT18","KRT7","KRT8","LYN","MAP2K1","MAP3K12","MAP3K14","MAP4K1","MAPK1","MAPK14","MAPK8IP1","MAPK8IP2","MAST1","MATR3","MDH1","MEGF6","MET","MGARP","MOB4","MUC1","NCAM1","NCK1","NCK2","NDN","NEDD4","NUMB","NUMBL","OAZ1","OLFM1","PAK1","PCNA","PDCD6IP","PDGFRB","PDK3","PIK3C2A","PIK3C2B","PIK3R1","PIK3R2","PIK3R3","PIN4","PITPNA","PKIA","PLCG1","PLCG2","PLD1","PLD2","PLEC","PLSCR1","PPP5C","PRCC","PRDX1","PRKACA","PRKAR1A","PRKAR1B","PRKCA","PRKCB","PRKCD","PRKD1","PSMA7","PSMD4","PTGDS","PTK2","PTK2B","PTK6","PTPN1","PTPN11","PTPN12","PTPN2","PTPN6","PTPRJ","PTPRS","RAB3A","RAP1GDS1","RASA1","RGS16","RIN2","RIPK1","RNF115","RNF126","RQCD1","RUSC2","S100A7","S100A9","SCAMP1","SCAMP3","SEC13","SEPP1","SFN","SGSM2","SH2B1","SH2B3","SH2D1A","SH2D2A","SH2D3A","SH3BGRL3","SHC1","SHC2","SHC3","SLA","SLC3A2","SLC9A3R1","SMURF2","SNRPD2","SNX1","SNX2","SNX4","SNX6","SOCS1","SOCS3","SOCS4","SOCS5","SOCS7","SORBS2","SOS1","SOS2","SPARCL1","SPCS2","SRC","STAM2","STAT1","STAT3","STAT5A","STAT5B","STIP1","STUB1","SYK","TGFA","TJP1","TLN1","TMCO3","TNC","TNK2","TNS4","TPI1","TPM1","TRPV1","TUBA1A","TUBA4A","UBB","UBE2V2","UCHL1","UROD","VAPA","VAV1","VAV2","VAV3","XRCC6","YWHAZ","ZAP70","ZNF259")
# from http://biograph.be/concept/show/C1414313
Paper_EGFR_gene <- read.delim(paste(FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/Paper_EGFR_gene.txt",sep=""))
EGFR_genes2<-as.character(Paper_EGFR_gene[,1])
#http://en.wikipedia.org/wiki/Epidermal_growth_factor_receptor
wiki<-c("AR","ARF4","CAV1","CAV3","CBL","CBLB","CBLC","CDC25A","CRK","CTNNB1","DCN","EGF","GRB14","Grb2","JAK2","MUC1","NCK1","NCK2","PKCA","PLCG1","PLSCR1","PTPN1","PTPN11","PTPN6","PTPRK","SH2D3A","SH3KBP1","SHC1","SOS1","Src","STAT1","STAT3","STAT5A","UBC","WAS")
#http://hintdb.hgc.jp/htp/proteins/P00533.html
HitPredict<-c("1433Z","ERBB2","CBL","EGF","P3C2B","TGFA","STAM2","RIN1","EGFR","CRK","DUS3","MK14","ARF4","M4K1","PLS1","PTN1","CASP1","PTPRS","P85B","PAK1","ACTS","MET","CEACAM1","TNR6","MP2K1","MK01","ELF3","PLD1","KPCD1","ANXA1","EZRI","DCN","HBEGF","CAV2","CLH2","HXC10","LYN","S10A7","S10A9","SEC13","AT1A1","ICAM1","CD166","CMTM8","SOCS5","P3C2A","SH3BGRL3","ADRBK1","PDC6I","NHRF1","ESR1","ITA5","GNAI2","KU70","MUC1","P85A","KPCD","TNC","COX2","PIPNA","EPS8","BTC","AMH","CAV3","SNX1","GRB7","GRB14","EREG","KCC2G","RGS16","ZPR1","PKIA","SCAM1","SCAMP3","KCC2A","PTN2","ERRFI1","DEGS1","AGTR1","CYH2","PTN12","HDAC6","GELS","HSP71","ARPC5","SH3L1","UB2V2","SPCS2","DYHC1","M3K12","PRDX1","AKAP12","FKBP4","CNTN2","GBRL2","TPIS","OAZ1","MATR3","MDHC","HSP74","STIP1","PUR9","EF1G","COF1","RAB3A","AATC","FAAA","ENPL","UROD","ALDOA","AHSA1","SRBS2","HEXI1","Q59EJ3","PPP5","ADRM1","PDK3","FKBP8","CO9A3","DCTN2","GPM6B","AHNAK","ARF6","PSMD4","SEPP1","PTGDS","PRKAR1B","CBR1","CSRP1","UCHL1","PSA7","MEGF6","AIP","MOBL3","PIN4","PLPL2","DNJC4","NECD","CD049","RUSC2","KCTD9","HS90A","Q59G22","Q53GZ6","SNP25","Q4QQI8","PLCH1","MAST1","STUB1","NOE1","AP2M1","PRCC","Q8N4S1","CEND","SH321","SYN1","Q59GR8","TMCO3","LST2","ATPG","SMD2","GRB2","YES","SHC1","STAT3","CDC2","FAK1","EPHA2","FAK2","BCAR1","BCAR3","ACK1","ERBB3","ERBB4","CD45","PTPRJ","STA5A","PTN11","CSK","KPCA","SHIP2","VAV","VAV2","STAT1","SRC","JAK2","PTK6","Q6PID4","NCK1","IGF1R","KAPCA","NCK2","PTN6","RASA1","GRB10","PLCG1","SH3K1","SOCS1","STA5B","RIPK1","M3K14","PGFRB","SHC3","FER","PLEC1","IQGA1","UBS3A","PTPRB","VAV3","SH23A","SOCS3","H31","CD59","S10AB","1433S","FHL2","1433T","CTND1","1433E","HSPB1","CTNA1","AT2A2","CBLB","EFNB1","IMB1","K2C1","HNRPF","H2B1C","PAXI","THIO","K2C6A","HNRPK","DNJA3","HSP7C","HNRH1","H14","K2C5","GRP78","K1C17","GRP75","PKP2","HS90B","EF1A1","DNJA1","CH60","CAV1","K1C14","EPS15","GAB2","CALM","G3P","4F2","CTNB1","PLAK","ANDR","LRSM1","CAMLG","ZO1","CBLC","CADH1","HD","DOK2","SOS1","SOS2","MPIP1","KAP0","PLD2","q02297-7","CD44","CD82","SNX6","SNX2","SNX4","CEBPB","GAB1","UBIQ","HGS","H31T","K2C7","K2C8","K1C18","EPIPL","UBC","NEDD8","SH3G2","LRIG1","PEX19","EI2BG","TBRG4","AT1B1","CLCA","SC5A1","AREG","AP2A1","VAPA","FINC","LEG3","TCTP","CCD50","UBP8","PCNA","EPN1","SPY2","NDFIP1","NDFIP2","LEG1","NEDD4","ACTA","SMUF2","OTU7B","TBA1A","ADT2","SGSM2","RABX5","TBA4A","WASP")
EGFR_genes<-sort(unique(c(wiki,EGFR_genes2)))
EGFR_genes<-sort(unique(c(wiki)))
EGFR_genes<-sort(unique(c(wiki,EGFR_genes2,HitPredict)))
length(EGFR_genes)

EGFR_genes.IDs<-IDs[IDs[,2]==EGFR_genes[1],]
for(i in 2:length(EGFR_genes)) EGFR_genes.IDs<-rbind(EGFR_genes.IDs,IDs[IDs[,2]==EGFR_genes[i],])


out=paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/List/")
#cor2=0
#X.TEST<-get.ln.data(SF[goodPV,]-cor2)     ### log 
#X.TEST<-get.log2.data(SF[goodPV,]-cor2)

X.TEST<-y2$E;colnames(X.TEST)<-targets$CellType  ####log2 !!!!
head(X.TEST)

get.ML.logLs.var.NEW<-function(x){  # with var-estimator
  
  NV.old<-function(d,mu,tau){
    #tau is precision
    return(dnorm(d,mu,sqrt(1/tau))) # f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))  
    #                                     return((tau/(2*pi))^(1/2)*exp(-0.5*tau*(x-mu)^2)) 
  }
  
  NV<-function(d,mu,tau){
    #tau is precision
    # f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))  
    return(-0.5*tau*(d-mu)^2)
  }
  
  get.var<-function(d,sub,MUs,N,i0,i1){
   # return(  ( sum((d[sub[[i0]]]-MUs[i0])^2) + sum((d[sub[[i1]]]-MUs[i1])^2)) / ( N[i0]+N[i1]-2 ) ) 
    return(  ( sum((d[sub[[i0]]]-MUs[i0])^2) + sum((d[sub[[i1]]]-MUs[i1])^2)) / ( N[i0]+N[i1]-1 ) ) 
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
    
    VARs<-c()
    VARs[1]<- 1/(N[1]-1)*sum( (d[a0]-MUs[1])^2)
    VARs[2]<-get.var(d,sub,MUs,N,2,3)
    VARs[3]<-get.var(d,sub,MUs,N,4,5)
    VARs[4]<-get.var(d,sub,MUs,N,6,7)
    SDs<-sqrt(VARs)
    TAUs<-1/VARs
    
    logL.old<-c()   
    logL.old[1]<-sum( log(NV.old(d[sub[[1]]],MUs[1],TAUs[1])) )
    logL.old[2]<-sum( log(NV.old(d[sub[[2]]],MUs[2],TAUs[2])) , log(NV.old(d[sub[[3]]],MUs[3],TAUs[2])) ) 
    logL.old[3]<-sum( log(NV.old(d[sub[[4]]],MUs[4],TAUs[3])) , log(NV.old(d[sub[[5]]],MUs[5],TAUs[3])) )  
    logL.old[4]<-sum( log(NV.old(d[sub[[6]]],MUs[6],TAUs[4])) , log(NV.old(d[sub[[7]]],MUs[7],TAUs[4])) )  
    
    #nach : http://de.wikipedia.org/wiki/Maximum-Likelihood-Methode#Allgemeine_Tests
    #logL<-c()   
    #logL[1]<-log( (TAUs[1]/(2*pi))^(N[1]/2) * exp(sum( (NV(d[sub[[1]]],MUs[1],TAUs[1]))  ) ) )    
    #logL[2]<-log( (TAUs[2]/(2*pi))^(N[1]/2) * exp(sum( (NV(d[sub[[2]]],MUs[2],TAUs[2])) , (NV(d[sub[[3]]],MUs[3],TAUs[2])) ) ) )    
    #logL[3]<-log( (TAUs[3]/(2*pi))^(N[1]/2) * exp(sum( (NV(d[sub[[4]]],MUs[4],TAUs[3])) , (NV(d[sub[[5]]],MUs[5],TAUs[3])) ) ) )    
    #logL[4]<-log( (TAUs[4]/(2*pi))^(N[1]/2) * exp(sum( (NV(d[sub[[6]]],MUs[6],TAUs[4])) , (NV(d[sub[[7]]],MUs[7],TAUs[4])) ) ) )

    logL<-c()   
    logL[1]<-( -(N[1]/2) * log(VARs[1]*2*pi) +(sum( (NV(d[sub[[1]]],MUs[1],TAUs[1]))  ) ) )    
    logL[2]<-( -(N[1]/2) * log(VARs[2]*2*pi) +(sum( (NV(d[sub[[2]]],MUs[2],TAUs[2])) , (NV(d[sub[[3]]],MUs[3],TAUs[2])) ) ) )    
    logL[3]<-( -(N[1]/2) * log(VARs[3]*2*pi) +(sum( (NV(d[sub[[4]]],MUs[4],TAUs[3])) , (NV(d[sub[[5]]],MUs[5],TAUs[3])) ) ) )    
    logL[4]<-( -(N[1]/2) * log(VARs[4]*2*pi) +(sum( (NV(d[sub[[6]]],MUs[6],TAUs[4])) , (NV(d[sub[[7]]],MUs[7],TAUs[4])) ) ) )
    
    #sum( 0.5*1/VARs[3]*(d[sub[[4]]]-MUs[4])^2 , 0.5*1/VARs*(d[sub[[5]]]-MUs[5])^2)
    
    #sum(TAUs[1]*(d[sub[[1]]]-MUs[1])^2)*log(6)*2
    
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
  
  logL<-apply(x,1,logLs,sub,N)
  logL<-restructur(logL)
  #print(head(logL))
  
  return(logL)
}

X.TEST_L<-get.ML.logLs.var.NEW(X.TEST)
#A <- make.IC(X.TEST_L,c(2,3),2) # AIC
#A.o<-A$IC[A$minV==3,]
#or<-order(A.o[,3],decreasing=F)
#A.o2<-A.o[or,]
#A.N<-as.character(rownames(A.o2))
#length(A.N)
range(X.TEST_L)

X.TEST["ILMN_1730999",]
X.TEST_L["ILMN_1730999",]
B$IC["ILMN_1730999",]

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

B <- make.IC(X.TEST_L,c(2,3),6) # BIC with 6 = number of values   
singleB<-function(){
    B3<-B$IC[B$minV==3,]
    or2<-order(B3[,3],decreasing=F)
    B3<-B3[or2,]
    B3.N<-as.character(rownames(B3))
    length(B3.N)    
    B4<-B$IC[B$minV==4,]
    or2<-order(B4[,4],decreasing=F)
    B4<-B4[or2,]
    B4.N<-as.character(rownames(B4))
    length(B4.N)     
    B1<-B$IC[B$minV==1,]
    or2<-order(B1[,1],decreasing=F)
    B1<-B1[or2,]
    B1.N<-as.character(rownames(B1))
    length(B1.N)    
    B2<-B$IC[B$minV==2,]
    or2<-order(B2[,2],decreasing=F)
    B2<-B2[or2,]
    B2.N<-as.character(rownames(B2))
    length(B2.N)     
}
singleB()
sort(IDs[B3.N,2])
sort(IDs[B4.N,2])
sort(IDs[B2.N,2])
B$IC[c(as.character(IDs[IDs[,2]=="TPR",1])),]
B$IC["ILMN_1730999",]
EGFR_genes.IDs[intersect(B3.N,EGFR_genes.IDs[,1]),2]
EGFR_genes.IDs[intersect(B4.N,EGFR_genes.IDs[,1]),2]


#########################
Prob<-function(){
  # Posterior P(m|x ) aus BIC berechnen 
  #-> BIC=-2 ln(P(x|m)) =>  P(x|m)= exp(-BIC/2)       => function q
  # P(m|x) = (P(x|m) * P(m) [bzw pi]) / sum_m Zähler  => function w

  qB<-function(B){
    return(exp(-B/2))
  }

  P_x_m<-apply(B$IC,MARGIN=1,FUN=qB)
  P_x_m<-t(log(P_x_m))

  #
  #old #pis<-c(0.90,0.008,0.07,0.022)
  #new#
  pis<-c(0.90,0.01,0.08,0.01)
  #pis<-c(0.90,0.02,0.06,0.02) ### can't be used with >0.85 => lose of  ILMN_1803236 CLCA2 with 0.8399344  !!!!!!!
  
  #pis<-c(0.25,0.25,0.25,0.25)

  w<-function(B,pis){
    t<-c()
    for(i in 1:4){
      t[i]<-B[i]*pis[i]
    }
    return(t/sum(t))
  }

  P_m_x<-apply(exp(P_x_m),MARGIN=1,FUN=w,pis   )
  P_m_x<-t(P_m_x)
  
  m1<-apply((P_x_m),MARGIN=1,FUN=maxIndex)
  m2<-apply((P_m_x),MARGIN=1,FUN=maxIndex)
  #m3<-apply(X.TEST_L,MARGIN=1,FUN=maxIndex)
  #m4<-apply(B$IC,MARGIN=1,FUN=minIndex)


  #table(apply((X.TEST_L),MARGIN=1,FUN=maxIndex))
  print(table(m1))
  print(table(m2))
  #print(table(m3))
  #print(table(m4))
  #input<-f.input3(ALT_125N,FC_125_BG.N,B.c.SIG$N,c("FC_1.25","FC_1.25_BG","BIC.var"))  ;venn(input) ;
  #input<-f.input2(rownames(P_m_x_c),B.N,c("FC_1.25_BG","BIC.var"))  ;venn(input) ;
  
  P_m_x_t<-P_m_x[m2==1,]
  P_m_x_a<-P_m_x_t[order(P_m_x_t[,1],decreasing=T),]
  P_m_x_t<-P_m_x[m2==2,]
  P_m_x_b<-P_m_x_t[order(P_m_x_t[,2],decreasing=T),]
  P_m_x_t<-P_m_x[m2==3,]
  P_m_x_c<-P_m_x_t[order(P_m_x_t[,3],decreasing=T),]
  P_m_x_t<-P_m_x[m2==4,]
  P_m_x_d<-P_m_x_t[order(P_m_x_t[,4],decreasing=T),]

  sum(P_m_x_c[,3]>0.85)
}

dim(P_m_x_a)
dim(P_m_x_b)
dim(P_m_x_c)
dim(P_m_x_d)

P_m_x_a.85<-P_m_x_a[P_m_x_a[,1]>0.85,]
P_m_x_b.85<-P_m_x_b[P_m_x_b[,2]>0.85,]
P_m_x_c.85<-P_m_x_c[P_m_x_c[,3]>0.85,]
P_m_x_d.85<-P_m_x_d[P_m_x_d[,4]>0.85,]
dim(P_m_x_c.85)

hist(P_m_x[,3])
####################################
a0<-c(1,2,3,4,5,6);   b0<-c(1,3,5);       b1<-c(2,4,6);         c0<-c(1,3,4,5);      c1<-c(2,6);         d0<-c(1,3,4,5,6);    d1<-c(2);            
ADD_FC<-function(){
  a0<-c(1,2,3,4,5,6);   b0<-c(1,3,5);       b1<-c(2,4,6);         c0<-c(1,3,4,5);      c1<-c(2,6);         d0<-c(1,3,4,5,6);    d1<-c(2);            
    
  addFC<-function(RNs,c0,c1){
    if(length(c1)>1){
      FC<-exp(rowMeans(X.TEST[RNs,c1])-rowMeans(X.TEST[RNs,c0]))    
    }else{
      FC<-exp( X.TEST[RNs,c1]-rowMeans(X.TEST[RNs,c0]))          
    }
    FC<-cbind(exp(X.TEST[RNs,c0]),exp(X.TEST[RNs,c1]),FC,log2(FC))
    FC<-FC[order(abs(FC[,8]),decreasing=T),]
    colnames(FC)[8]<-"log2FC"
    return(FC)
  }
  addFC_log2<-function(RNs,c0,c1){
    if(length(c1)>1){
      FC<-2^(rowMeans(X.TEST[RNs,c1])-rowMeans(X.TEST[RNs,c0]))    
    }else{
      FC<-2^( X.TEST[RNs,c1]-rowMeans(X.TEST[RNs,c0]))          
    }
    FC<-cbind(2^(X.TEST[RNs,c0]),2^(X.TEST[RNs,c1]),FC,log2(FC))
    FC<-FC[order(abs(FC[,8]),decreasing=T),]
    colnames(FC)[8]<-"log2FC"
    return(FC)
  }
  #FC.B2<-addFC(rownames(P_m_x_b),b0,b1)
  #FC.B3<-addFC(rownames(P_m_x_c),c0,c1)
  #FC.B4<-addFC(rownames(P_m_x_d),d0,d1)
  FC.B2<-addFC_log2(rownames(P_m_x_b),b0,b1)
  FC.B3<-addFC_log2(rownames(P_m_x_c),c0,c1)
  FC.B4<-addFC_log2(rownames(P_m_x_d),d0,d1)
  
  FC.B2<-cbind(FC.B2,P_m_x_b[rownames(FC.B2),2])
  colnames(FC.B2)[9]<-"Prob_for_group"
  FC.B3<-cbind(FC.B3,P_m_x_c[rownames(FC.B3),3])
  colnames(FC.B3)[9]<-"Prob_for_group"
  FC.B4<-cbind(FC.B4,P_m_x_d[rownames(FC.B4),4])
  colnames(FC.B4)[9]<-"Prob_for_group"
  B2.N<-rownames(FC.B2)
  B3.N<-rownames(FC.B3)
  B4.N<-rownames(FC.B4)
  
  YYY<-cbind(FC.B3[order(FC.B3[,8]),],P_m_x_c[rownames(FC.B3[order(FC.B3[,8]),]),])
  yyy<-rownames(mfc.c.gene3.c)
  YYY<-YYY[yyy,]
  YYY[order(YYY[,8]),]
  
  head(FC.B3[order(FC.B3[,"Prob_for_group"],decreasing=T), ],700)
  
    #### out 
  write_out<-function(){
    
    FC.B2.85<-FC.B2[FC.B2[,9]>0.85,]
    FC.B3.85<-FC.B3[FC.B3[,9]>0.85,]; FC.B3.85<-FC.B3.85[order(FC.B3.85[,9],decreasing=T),]
    FC.B4.85<-FC.B4[FC.B4[,9]>0.85,]
    length(FC.B3[,1])
    
    FC.B2.ID<-as.matrix(cbind(rownames(FC.B2), as.character(IDs[rownames(FC.B2),2]), round(FC.B2,4)))
    colnames(FC.B2.ID)[c(1,2)]<-c("ID","Sym")
    FC.B3.ID<-as.matrix(cbind(rownames(FC.B3), as.character(IDs[rownames(FC.B3),2]), round(FC.B3,4)))
    colnames(FC.B3.ID)[c(1,2)]<-c("ID","Sym")
    FC.B4.ID<-as.matrix(cbind(rownames(FC.B4), as.character(IDs[rownames(FC.B4),2]), round(FC.B4,4)))
    colnames(FC.B4.ID)[c(1,2)]<-c("ID","Sym")
    
    out=paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/List/")
    write.table(FC.B2.ID,paste(sep="",out,"PAPER_N1_newBIC_LIMMA_qu_bkgd_BIC_grB_AllInfo",paste(sep="_",pis[1],pis[2],pis[3],pis[4]),".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    write.table(FC.B3.ID,paste(sep="",out,"PAPER_N1_newBIC_LIMMA_qu_bkgd_BIC_grC_AllInfo",paste(sep="_",pis[1],pis[2],pis[3],pis[4]),".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    write.table(FC.B4.ID,paste(sep="",out,"PAPER_N1_newBIC_LIMMA_qu_bkgd_BIC_grD_AllInfo",paste(sep="_",pis[1],pis[2],pis[3],pis[4]),".txt"),row.names=F,col.names=T,quote=F,sep="\t")
  
    write.table(FC.B2.ID,paste(sep="",out,"PAPER_N1_newBIC_LIMMA_qu_bkgd_BIC_grB_AllInfo_prob85",paste(sep="_",pis[1],pis[2],pis[3],pis[4]),".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    write.table(FC.B3.ID,paste(sep="",out,"PAPER_N1_newBIC_LIMMA_qu_bkgd_BIC_grC_AllInfo_prob85",paste(sep="_",pis[1],pis[2],pis[3],pis[4]),".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    write.table(FC.B4.ID,paste(sep="",out,"PAPER_N1_newBIC_LIMMA_qu_bkgd_BIC_grD_AllInfo_prob85",paste(sep="_",pis[1],pis[2],pis[3],pis[4]),".txt"),row.names=F,col.names=T,quote=F,sep="\t")
    
  }
  EGF.predict.IDs<-c()
  EGF.predict.IDs<-c(B2.N,B3.N,B4.N)
  tmp<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),EGF.predict.IDs),c(1,2)]

  tmp<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),B3.N),c(1,2)]
  FC.B3[rownames(tmp),]  
  tmp<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),B4.N),c(1,2)]
  FC.B4[rownames(tmp),]  
  tmp<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),B2.N),c(1,2)]
  FC.B2[rownames(tmp),]  
}

######## 
RawQPCR <- read.csv(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/qPCR/MyData.csv"))
qgenes<-c("TPR","CKAP2L","KIF5C","ROCK1","BPNT1","GALNS","GLIPR2","KEL","ALDH4A1","CDCP1","CLCA2")
#qgenes<-c("CKAP2L", "ROCK1", "TPR","ALDH4A1", "CLCA2","GALNS")

IDs<-cbind(rownames(x$E),x$genes);colnames(IDs)<-c("ILMN","SYMBOL")
#IDs[IDs["SYMBOL"]=="TPR",]

qgenesIDs<-list()
for(i in 1:length(qgenes)){
  qgenesIDs[[i]]<-as.character(IDs[which(qgenes[i]==IDs[,"SYMBOL"]),1])
}
names(qgenesIDs)<-qgenes
#intersect(c(unlist(qgenesIDs),"ILMN_1738468","ILMN_1712718","ILMN_1737949","ILMN_1737949"),EGF.predict.IDs)

qGenes<-c()
qGenes[[1]] <-    cbind(RawQPCR[1:6,c(3,4)],RawQPCR[13:18,c(3,4)],RawQPCR[25:30,c(3,4)])
qGenes[[2]]<-     cbind(RawQPCR[1:6,c(5,6)],RawQPCR[13:18,c(5,6)],RawQPCR[25:30,c(5,6)])
qGenes[[3]]<-     cbind(RawQPCR[1:6,c(7,8)],RawQPCR[13:18,c(7,8)],RawQPCR[25:30,c(7,8)])
qGenes[[4]]<-     cbind(RawQPCR[1:6,c(9,10)],RawQPCR[13:18,c(9,10)],RawQPCR[25:30,c(9,10)])
qGenes[[5]] <-    cbind(RawQPCR[1:6,c(11,12)],RawQPCR[13:18,c(11,12)],RawQPCR[25:30,c(11,12)])
qGenes[[6]] <-    cbind(RawQPCR[1:6,c(13,14)],RawQPCR[13:18,c(13,14)],RawQPCR[25:30,c(13,14)])
qGenes[[7]] <-    cbind(RawQPCR[1:6,c(15,16)],RawQPCR[13:18,c(15,16)],RawQPCR[25:30,c(15,16)])
qGenes[[8]] <-    cbind(RawQPCR[1:6,c(17,18)],RawQPCR[13:18,c(17,18)],RawQPCR[25:30,c(17,18)])
qGenes[[9]] <-    cbind(RawQPCR[1:6,c(19,20)],RawQPCR[13:18,c(19,20)],RawQPCR[25:30,c(19,20)])
qGenes[[10]] <-   cbind(RawQPCR[1:6,c(21,22)],RawQPCR[13:18,c(21,22)],RawQPCR[25:30,c(21,22)])
qGenes[[11]] <-   cbind(RawQPCR[1:6,c(23,24)],RawQPCR[13:18,c(23,24)],RawQPCR[25:30,c(22,23)])
names(qGenes)<-qgenes

qGenes.M<-c()
for(i in 1:length(qgenes)){
  qGenes.M[[i]]<-cbind(rowMeans(log2(qGenes[[i]][c(1,3,5)])),rowMeans(log2(qGenes[[i]][c(2,4,6)])),  log2(rowMeans(qGenes[[i]][c(2,4,6)])) )
  colnames(qGenes.M[[i]])<-c("log2geoMeanCT","log2geoMeanREL","log2MeanREL")
}
names(qGenes.M)<-qgenes

#qgenes<-c()
qGenes.FC<-list()
for(i in 1:length(qgenes)){
  qGenes.FC[[i]]<-colMeans(qGenes.M[[i]][c1,])- colMeans(qGenes.M[[i]][c0,])
}
names(qGenes.FC)<-qgenes
###
Paper.qgenes<-c("CKAP2L", "ROCK1", "TPR","ALDH4A1", "CLCA2","GALNS")
Paper.qgenesIDs<-c(qgenesIDs[[2]],qgenesIDs[[4]][2],qgenesIDs[[1]][2],qgenesIDs[[9]][3],qgenesIDs[[11]][1],qgenesIDs[[6]][1])
Paper.qgenes[!Paper.qgenesIDs %in% rownames(FC.B3)]
Paper.qgenes[Paper.qgenesIDs %in% rownames(FC.B3)]
Paper.qgenes.FC<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5)])
Paper.qgenes.FC2<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,3]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,3]))[,c(1,3,4,5)])

Paper.qgenes.FC.bb<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,4,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,5)])
Paper.qgenes.FC.cc<-rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2,6)]) - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5)])
Paper.qgenes.FC.dd<-t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(2)] - rowMeans(t(sapply(Paper.qgenes,function(x) qGenes.M[[x]][,2]))[,c(1,3,4,5,6)])


paper.1<-as.matrix(FC.B3[rownames(FC.B3) %in% Paper.qgenesIDs,])
paper.1<-paper.1[Paper.qgenesIDs,]
rownames(paper.1)<-Paper.qgenes
paper.1<-cbind(paper.1,Paper.qgenes.FC)
paper.1<-cbind(paper.1,Paper.qgenes.FC2)
P_m_x_c[rownames(P_m_x_c)%in% Paper.qgenesIDs,]


pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new1.pdf"),height=7,width=8)
barplot(log2(t(paper.1[,1:6])),ylim=c(0,10),beside=T,col=brewer.pal(11,"RdYlBu")[c(1:4,10,9)],legend.text=paste(sep=".","x",c(1,3,4,5,2,6)),ylab="logarithmic expression levels")#,ylab="avg microarray signal")
dev.off()

pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new2.pdf"),height=7,width=7)
barplot(t(paper.1[,c(1,5,2,3,4,6)]),beside=T,col=brewer.pal(11,"RdYlBu")[c(1,10,2:4,9)],legend.text=paste(sep=".","x",c(1:6)),ylab="avg microarray signal")
dev.off()


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
          axis.title.y = element_text(face="bold", colour=1,angle=90 ,vjust=0.5, size=base_size),
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

pdf(paste(sep="",out,"PAPER_qPCR_Array_expression_new4_limma.pdf"),height=7,width=10)
g1<-ggplot(df,aes(x=factor(X2),y=value,fill=C)) + geom_bar(position="dodge",stat="identity",width=.5)+facet_wrap(~ X1)+  #,title = paste(sep=" ","Reads after quality control"))+
  labs(x = "", y = "Logarithmic expression levels")+scale_fill_manual(values = cols,name=" Indicator variable",labels = list("g = 0","g = 1" )) + ylim(c(0,10))+ thememap(14,0.6) +
  scale_x_discrete(labels=c(expression(x[1]),expression(x[3]),expression(x[4]),expression(x[5]),expression(x[2]),expression(x[6])))+
  theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")

g1 #+ theme(legend.justification=c(.9,-1.7), legend.position=c(1,0)) #theme(legend.position="bottom")
dev.off()
###

qGenes.FC.comp<-c()
qGenes.FC.comp<-c(qGenes.FC[[1]][c(2,3)],FC.B3[c(qgenesIDs[[1]])[2],"log2FC"] )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[2]][c(2,3)],FC.B3[c(qgenesIDs[[2]])   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[4]][c(2,3)],FC.B3[c(qgenesIDs[[4]])[2],"log2FC"] ) )
#qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[3]][c(2,3)],log2(FC[c(qgenesIDs[[3]])[1] ])) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[5]][c(2,3)],FC.B3[c(qgenesIDs[[5]])   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[6]][c(2,3)],FC.B3[c(qgenesIDs[[6]])   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[7]][c(2,3)],FC.B3[c(qgenesIDs[[7]])   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[8]][c(2,3)],FC.B3[c(qgenesIDs[[8]])   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[9]][c(2,3)],FC.B3[intersect(B3.N,qgenesIDs[[9]])[1]   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[10]][c(2,3)],FC.B3[intersect(B3.N,qgenesIDs[[10]])[1]   ,"log2FC"] ) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,c(qGenes.FC[[11]][c(2,3)],FC.B3[intersect(B3.N,qgenesIDs[[11]])[1]   ,"log2FC"] ) )
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




grC.ID<-IDs[rownames(P_m_x_c),c(1,2)]
#grC.ID.qPCR<-c()
#for(g in qgenes[c(-3)]){
#  grC.ID.qPCR<-rbind(grC.ID.qPCR,grC.ID[which(grC.ID[,2]==g),])
#}
#grC.ID.qPCR<-grC.ID.qPCR[-c(6,9,11),]
grC.ID.qPCR<-grC.ID[Paper.qgenesIDs,]

grC.A<-X.TEST[as.character(rownames(grC.ID.qPCR)),]

grC.Q<-c()
for(g in as.character(grC.ID.qPCR[,2])){
  grC.Q<-rbind(grC.Q,qGenes.M[[g]][,2])
}
rownames(grC.Q)<-rownames(grC.A)

grC.Q.log2<-grC.Q #####log2(exp(grC.Q))
grC.A.log2<-grC.A #log2(exp(grC.A))

grC.Q.FC<-rowMeans(grC.Q.log2[,c1])- rowMeans(grC.Q.log2[,c0])
grC.A.FC<-rowMeans(grC.A.log2[,c1])- rowMeans(grC.A.log2[,c0])

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

pdf(paste(sep="",out,"PAPER_qPCR_FC_new_Limma.pdf"),height=6,width=8)
g1<-ggplot(df2,aes(x=factor(X1),y=value,fill=X2)) + geom_bar(position="dodge",stat="identity",width=.5)+labs(x = "", y = "Log2-fold change for group c") + ylim(c(-1.5,1.5))+ scale_fill_manual(values = cols,name="") +
  thememap(14,0.6) +theme(legend.background = element_rect(fill="grey90", size=.5, linetype="dotted"))+theme(legend.position="bottom")
g1  
#g1+ theme(legend.justification=c(.9,-1.7), legend.position=c(1,0)) #theme(legend.position="bottom")
dev.off()
cor(grC.FC)
plot(grC.FC[,1],grC.FC[,2],ylim=c(-1.5,1.5),xlim=c(-1.5,1.5))
abline(lm(grC.FC[,1]~grC.FC[,2]))

cor(grC.FC,method="spearman")

P_m_x_c[]
ooooo<-P_m_x_c[take,]



################# FoldChange Test
FC1_KO<-X.TEST[,2]-X.TEST[,1]
FC1_KO.M<-apply(X.TEST[,c(1,2)],MARGIN=1,FUN=mean)
FC1_ALL<-X.TEST[,4]-X.TEST[,3]
FC1_ALL.M<-apply(X.TEST[,c(3,4)],MARGIN=1,FUN=mean)
FC1_14<-X.TEST[,6]-X.TEST[,5]
FC1_14.M<-apply(X.TEST[,c(5,6)],MARGIN=1,FUN=mean)

FC1<-cbind(FC1_KO,FC1_ALL,FC1_14)
FC1.M<-cbind(FC1_KO.M,FC1_ALL.M,FC1_14.M)


ma.plot(FC1_KO.M,FC1_KO,cex=0.6,plot.method="smoothScatter",add.loess=T)
ma.plot(FC1_ALL.M,FC1_ALL,cex=0.6,plot.method="smoothScatter",add.loess=T)
ma.plot(FC1_14.M,FC1_14,cex=0.6,plot.method="smoothScatter",add.loess=T)


library(ggplot2)
MA.df <- data.frame(FC1_KO.M,FC1_KO);
colnames(MA.df) <- c("mean","diff");
dcolor<-"green"; lcolor<-"red";sp=0.2
MA.plot <- ggplot(data=MA.df,aes(x=mean,y=diff),na.rm=TRUE) +
  geom_point(alpha=.25, color="black", size=1.3) +  geom_point(alpha=.25, color=dcolor,size=.3) +  ylim(c(-5,5)) +
  #stat_smooth(method="loess",family="gaussian",span=sp,colour=lcolor, size=1,fill=lcolor) #+
  geom_abline(intercept=1,slope=0,colour="blue",size=.5) +   geom_abline(intercept=-1,slope=0,colour="blue",size=.5) 
MA.plot
hist(abs(MA.df$diff))



Cut.off.FC<-.5
#Cut.off.FC<-.8
#Cut.off.FC<-.75
Cut.off.FC<-1

GR.A_KO<- rownames(FC1)[abs(FC1[,1])  < Cut.off.FC ]
GR.A_ALL<-rownames(FC1)[abs(FC1[,2])  < Cut.off.FC ]
GR.A_14<- rownames(FC1)[abs(FC1[,3])  < Cut.off.FC ]
GR.A<-intersect(GR.A_ALL,intersect(GR.A_KO,GR.A_14))
intersect(Paper.qgenesIDs,GR.A)

GR.B_KO<- rownames(FC1)[abs(FC1[,1])  > Cut.off.FC ]
GR.B_ALL<-rownames(FC1)[abs(FC1[,2])  > Cut.off.FC ]
GR.B_14<- rownames(FC1)[abs(FC1[,3])  > Cut.off.FC ]
GR.B<-intersect(GR.B_ALL,intersect(GR.B_KO,GR.B_14))
intersect(Paper.qgenesIDs,GR.B)
  
GR.C_KO <-rownames(FC1)[abs(FC1[,1])  > Cut.off.FC ]
GR.C_ALL<-rownames(FC1)[abs(FC1[,2])  < Cut.off.FC ]
GR.C_14 <-rownames(FC1)[abs(FC1[,3])  > Cut.off.FC ]
GR.C<-intersect(GR.C_ALL,intersect(GR.C_KO,GR.C_14))
intersect(Paper.qgenesIDs,GR.C)

GR.D_KO <-rownames(FC1)[abs(FC1[,1]) > Cut.off.FC ]
GR.D_ALL<-rownames(FC1)[abs(FC1[,2]) < Cut.off.FC ]
GR.D_14 <-rownames(FC1)[abs(FC1[,3]) < Cut.off.FC ]
GR.D<-intersect(GR.D_ALL,intersect(GR.D_KO,GR.D_14))
intersect(Paper.qgenesIDs,GR.D)

length(unique(c(GR.A,GR.B,GR.C,GR.D))) / dim(X.TEST)[1]

pdf(paste(out,"Foldchange_Method_log2FC_",Cut.off.FC,".pdf"))
  plot(FC1.M[Paper.qgenesIDs,1],FC1[Paper.qgenesIDs,1],pch=21,xlim=c(4,10),ylim=c(-2,2),col=brewer.pal(8,"Dark2")[c(1:5,7)],ylab="M",xlab="A")
  points(FC1.M[Paper.qgenesIDs,2],FC1[Paper.qgenesIDs,2],pch=4,col=brewer.pal(8,"Dark2")[c(1:5,7)])
  points(FC1.M[Paper.qgenesIDs,3],FC1[Paper.qgenesIDs,3],pch=13,col=brewer.pal(8,"Dark2")[c(1:5,7)])
  legend("topleft",col=1,pch=c(21,4,13),c("KO","ALL","14"))
  legend("topright",col=brewer.pal(8,"Dark2")[c(1:5,7)],lty=1,Paper.qgenes)
  abline(h=Cut.off.FC);abline(h=(-1*Cut.off.FC))
  
  barplot(t(FC1[Paper.qgenesIDs,]),beside=T,ylim=c(-2,2),names.arg=Paper.qgenes,legend.text=c("KO","ALL","14"))
  abline(h=Cut.off.FC);abline(h=(-1*Cut.off.FC))
  
  f.input5(GR.A,GR.B,GR.C,GR.D,Paper.qgenesIDs,name=c("A","B","C","D","qPCR"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_a.85),name=c("A","B","C","D","MY.A.85"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_b.85),name=c("A","B","C","D","MY.B.85"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_c.85),name=c("A","B","C","D","MY.C.85"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_d.85),name=c("A","B","C","D","MY.D.85"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_a),name=c("A","B","C","D","MY.A"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_b),name=c("A","B","C","D","MY.B"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_c),name=c("A","B","C","D","MY.C"))
  f.input5(GR.A,GR.B,GR.C,GR.D,rownames(P_m_x_d),name=c("A","B","C","D","MY.D"))
dev.off()

dim(P_m_x_a);dim(P_m_x_b);dim(P_m_x_c);dim(P_m_x_d);
dim(P_m_x_a.85);dim(P_m_x_b.85);dim(P_m_x_c.85);dim(P_m_x_d.85);
dim(y2$E)


# GR.A.M<- rownames(FC1)[ abs(rowMeans(FC1) ) < Cut.off.FC]
# intersect(Paper.qgenesIDs,GR.A.M)
# GR.B.M<- rownames(FC1)[ abs(rowMeans(FC1) ) > Cut.off.FC]
# intersect(Paper.qgenesIDs,GR.B.M)
# GR.C_0 <-rownames(FC1)[abs(rowMeans(FC1[,c(1,3)]) )  > Cut.off.FC ]
# GR.C_1 <-rownames(FC1)[abs(FC1[,2])  < Cut.off.FC ]
# GR.C.M<-intersect(GR.C_0,GR.C_1)
# intersect(Paper.qgenesIDs,GR.C.M)
# GR.D_1 <-rownames(FC1)[abs(rowMeans(FC1[,c(2,3)]) )  < Cut.off.FC ]
# GR.D_0 <-rownames(FC1)[abs(FC1[,1])  > Cut.off.FC ]
# GR.D.M<-intersect(GR.D_0,GR.D_1)
# intersect(Paper.qgenesIDs,GR.D.M)
# f.input5(GR.A.M, GR.B.M, GR.C.M, GR.D.M,rownames(P_m_x_a.85),name=c("A","B","FC 1\nC","D","MY.A"))
# f.input5(GR.A.M, GR.B.M, GR.C.M, GR.D.M,rownames(P_m_x_b.85),name=c("A","B","FC 1\nC","D","MY.B"))
# f.input5(GR.A.M, GR.B.M, GR.C.M, GR.D.M,rownames(P_m_x_c.85),name=c("A","B","FC 1\nC","D","MY.C"))
# f.input5(GR.A.M, GR.B.M, GR.C.M, GR.D.M,rownames(P_m_x_d.85),name=c("A","B","FC 1\nC","D","MY.D"))
# f.input5(GR.A.M, GR.B.M, GR.C.M, GR.D.M,Paper.qgenesIDs,name=c("A","B","FC 1\nC","D","qPCR"))


################# t-Test
#Array<-X.TEST
#C1<-c(1,3,5);C2<-c(2,4,6);
t.Test.array<-function(Array,C1,C2,IDs="",What=""){
  
  pVal<-c()
  FC<-c()
  log2FC<-c()
  for(i in 1:length(Array[,1])){
    pVal[i]<-t.test(Array[i,C1],Array[i,C2],var.equal=TRUE,paired=FALSE)$p.value ## used
    
    #pVal[i]<-t.test(Array[i,C1],Array[i,C2],var.equal=F,paired=T)$p.value
    
    #FC[i]<-log2DES(c(mean(Array[i,1:2]),mean(Array[i,3:4])))
    FC[i]<- mean(Array[i,C2])-mean(Array[i,C1])
    
  }
  log2FC<-FC
  
  df<-t.test(Array[i,C1],Array[i,C2],var.equal=TRUE)$parameter
  print(df)
  #library(multtest)
  require(multtest)
  adjp<-mt.rawp2adjp(pVal, proc=c("BH"), alpha = 0.05, na.rm = FALSE) ## default alpha = 0.05 ! 
  or.adjp<-adjp$adjp[order(adjp$index),]
  
  res<-cbind(Array,log2FC,pVal,or.adjp[,2])
  colnames(res)[length(res[1,])]<-"BH"
  print(head(res))
  
  

  #resALLout<-cbind(as.character(IDs[rownames(resALL),2]),resALL)
  #colnames(resALLout)[1]<-"SYMBOL"
  #OUT=paste(sep="",folderOUT,"t.Test_",What)
  #print(OUT)
  # write.table( add_IDS(resALLout),file=paste(sep="",OUT,"_resALL.txt"),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE,sep = "\t")
  
  
  return(res)
}  

t.GR.B<-t.Test.array(X.TEST,C1<-c(1,3,5),C2<-c(2,4,6));sum(t.GR.B[,"BH"]<0.05)
t.GR.C<-t.Test.array(X.TEST,C1<-c(1,3,4,5),C2<-c(2,6));sum(t.GR.C[,"BH"]<0.05)
t.GR.D<-t.Test.array(X.TEST,C1<-c(1,3,4,5,6),C2<-c(2));sum(t.GR.D[,"BH"]<0.05)

t.GR.B[Paper.qgenesIDs,]
t.GR.C[Paper.qgenesIDs,]
t.GR.D[Paper.qgenesIDs,]


################# t-Test Limma
t.test.limma<-function(y2,targets,AW=1){
rna<-factor(targets$CellType)
design <- model.matrix (~0 + rna)
colnames ( design ) <- levels (rna)
if(AW==1){
  aw <- arrayWeights (y2 , design )
  barplot(aw, xlab="Array", ylab="Weight", col="white", las=2); abline(h=1, lwd=1, lty=2)
  print(aw)#
}else{  
  aw<-c(1,1,1,1,1,1)
}
fit <- lmFit (y2, design,weights=aw)
contrasts <- makeContrasts(mu1-mu0, levels=design)  ## important for logFC calculation 
contr.fit <- eBayes ( contrasts.fit(fit , contrasts ))
print(summary(decideTests(contr.fit, method="global")))
print(topTable ( contr.fit , coef = 1,))
topT.aw<-topTable(contr.fit, coef=1,number=length(y2$E[,1]))
return(topT.aw)
}
targets$CellType=c("mu0","mu1","mu0","mu1","mu0","mu1")
limma.GR.B<-t.test.limma(y2,targets) ; limma_aw1.GR.B<-t.test.limma(y2,targets,AW=0)
targets$CellType=c("mu0","mu1","mu0","mu0","mu0","mu1")
limma.GR.C<-t.test.limma(y2,targets) ; limma_aw1.GR.C<-t.test.limma(y2,targets,AW=0)
targets$CellType=c("mu0","mu1","mu0","mu0","mu0","mu0")
limma.GR.D<-t.test.limma(y2,targets) ; limma_aw1.GR.D<-t.test.limma(y2,targets,AW=0)

bb<-limma.GR.B[Paper.qgenesIDs,]
cc<-limma.GR.C[Paper.qgenesIDs,]
dd<-limma.GR.D[Paper.qgenesIDs,]
sum(cc$adj.P.Val<bb$adj.P.Val)


aPvalue<-0.05
#aPvalue<-0.01
limma.GR.B.SIG<-limma.GR.B[limma.GR.B$adj.P.Val<aPvalue,]
limma.GR.C.SIG<-limma.GR.C[limma.GR.C$adj.P.Val<aPvalue,]
limma.GR.D.SIG<-limma.GR.D[limma.GR.D$adj.P.Val<aPvalue,]

pdf(paste(out,"Limma_ttest_Method_qval_",aPvalue,".pdf"))
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),as.character(rownames(limma.GR.D.SIG)),name=c("B","C","D"))
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),Paper.qgenesIDs,name=c("B","C","qPCR"))
  barplot(t(cbind(bb$adj.P.Val,cc$adj.P.Val,dd$adj.P.Val)),beside=T,names.arg=Paper.qgenes,legend.text=c("B","C","D"),ylab="adj.P.Val")
  abline(h=aPvalue);
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_a.85)),name=c("B","C","MY.A.85"))
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_b.85)),name=c("B","C","MY.B.85"))
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_c.85)),name=c("B","C","MY.C.85"))
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_d.85)),name=c("B","C","MY.D.85"))
  f.input3(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_c)),name=c("B","C","MY.C"))
  
  f.input5(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(P_m_x_a)),as.character(rownames(P_m_x_b)),as.character(rownames(P_m_x_c)),as.character(rownames(P_m_x_d)),name=c("B","MY.A","MY.B","MY.C","MY.D"))
  f.input5(as.character(rownames(limma.GR.B.SIG)),as.character(rownames(P_m_x_a.85)),as.character(rownames(P_m_x_b.85)),as.character(rownames(P_m_x_c.85)),as.character(rownames(P_m_x_d.85)),name=c("B","MY.A.85","MY.B.85","MY.C.85","MY.D.85"))
  
  f.input5(as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_a)),as.character(rownames(P_m_x_b)),as.character(rownames(P_m_x_c)),as.character(rownames(P_m_x_d)),name=c("C","MY.A","MY.B","MY.C","MY.D"))
  f.input5(as.character(rownames(limma.GR.C.SIG)),as.character(rownames(P_m_x_a.85)),as.character(rownames(P_m_x_b.85)),as.character(rownames(P_m_x_c.85)),as.character(rownames(P_m_x_d.85)),name=c("C","MY.A.85","MY.B.85","MY.C.85","MY.D.85"))

dev.off()

limma.GR.B.SIG<-limma.GR.B[limma.GR.B$adj.P.Val<0.05,]
limma.GR.B.SIG[Paper.qgenesIDs,]
paper.1

bb.y2<-rowMeans(y2$E[Paper.qgenesIDs,c(2,4,6)]) - rowMeans(y2$E[Paper.qgenesIDs,c(1,3,5)])
cc.y2<-rowMeans(y2$E[Paper.qgenesIDs,c(2,6)])   - rowMeans(y2$E[Paper.qgenesIDs,c(1,3,4,5)])
dd.y2<-          y2$E[Paper.qgenesIDs,c(2)]     - rowMeans(y2$E[Paper.qgenesIDs,c(1,3,4,5,6)])

cor(cbind(Paper.qgenes.FC.bb,bb$logFC,bb.y2))
cor(cbind(Paper.qgenes.FC.cc,cc$logFC,cc.y2))
cor(cbind(Paper.qgenes.FC.dd,dd$logFC,dd.y2))

apply(cbind(Paper.qgenes.FC.bb,bb$logFC,bb.y2),1,function(x) abs(x[1])-abs(x[3]))
apply(cbind(Paper.qgenes.FC.cc,cc$logFC,cc.y2),1,function(x) abs(x[1])-abs(x[3]))
cal.euc.dist<-function(A){d=(A[1]-A[2])/sqrt(2);return(abs(d))}
ED.bb<-apply(cbind(Paper.qgenes.FC.bb,bb$logFC,bb.y2)[,c(1,3)],MARGIN=1,FUN=cal.euc.dist)
ED.cc<-apply(cbind(Paper.qgenes.FC.cc,cc$logFC,cc.y2)[,c(1,3)],MARGIN=1,FUN=cal.euc.dist)
mean(ED.cc)

FC.B2.85<-FC.B2[FC.B2[,9]>0.85,]
FC.B3.85<-FC.B3[FC.B3[,9]>0.85,]; FC.B3.85<-FC.B3.85[order(FC.B3.85[,9],decreasing=T),]
FC.B4.85<-FC.B4[FC.B4[,9]>0.85,]

venn(f.input3(
  rownames(limma.GR.C[limma.GR.C$adj.P.Val<0.05,]),
  rownames(FC.B3.85),
  Paper.qgenesIDs,c("limma","BIC","qPCR")))

venn(f.input3(
  rownames(limma.GR.B[limma.GR.B$adj.P.Val<0.05,]),
  rownames(FC.B2.85),
  Paper.qgenesIDs,c("limma","BIC","qPCR")))

venn(f.input3(
  rownames(limma.GR.D[limma.GR.D$adj.P.Val<0.05,]),
  rownames(FC.B4.85),
  Paper.qgenesIDs,c("limma","BIC","qPCR")))

venn(f.input3(
  rownames(limma.GR.B[limma.GR.B$adj.P.Val<0.05,]),
  rownames(limma.GR.C[limma.GR.C$adj.P.Val<0.05,]),
  rownames(limma.GR.D[limma.GR.D$adj.P.Val<0.05,]),
  c("limma.B","limma.C","limma.D")))

venn(f.input3(
  rownames(limma.GR.B[limma.GR.B$adj.P.Val<0.05,]),
  rownames(limma.GR.C[limma.GR.C$adj.P.Val<0.05,]),
  Paper.qgenesIDs,c("limma.B","limma.C","qPCR")))

f.input5(
  rownames(limma.GR.B[limma.GR.B$adj.P.Val<0.05,]),
  rownames(FC.B2.85),
  rownames(FC.B3.85),
  rownames(limma.GR.C[limma.GR.C$adj.P.Val<0.05,]),
  Paper.qgenesIDs,c("limma.B","BIC.B","BIC.C","limma.C","qPCR"))

qqqq<-setdiff(Paper.qgenesIDs,intersect(rownames(limma.GR.B[limma.GR.B$adj.P.Val<0.05,]),
          intersect(rownames(limma.GR.C[limma.GR.C$adj.P.Val<0.05,]),Paper.qgenesIDs)))






################# TEST for total design
targets$CellType
EGF<-factor(c("wE","E","wE","E","wE","E"))
RNAi<-factor(c("KO","KO","ALL","ALL","FZ","FZ"))
design <- model.matrix (~EGF*RNAi)
colnames(design)
if(AW==1){
  aw <- arrayWeights (y2 , design )
  barplot(aw, xlab="Array", ylab="Weight", col="white", las=2); abline(h=1, lwd=1, lty=2)
  print(aw)#
}else{  
  aw<-c(1,1,1,1,1,1)
}
fit <- lmFit (y2, design,weights=aw)
contrasts <- makeContrasts(RNAiFZ-EGFwE, levels=design)  ## important for logFC calculation 
contr.fit <- eBayes ( fit,fit=)
  contrasts.fit(fit , contrasts ))
print(summary(decideTests(contr.fit, method="global")))
print(topTable ( contr.fit , coef = 1,))
topT.aw<-topTable(contr.fit, coef=1,number=length(y2$E[,1]))

results <- classifyTestsF(contr.fit)
vennCounts(results)
################# TEST for total design
