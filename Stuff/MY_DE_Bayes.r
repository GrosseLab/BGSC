FOLDER="~/Promotion"
#FOLDER="~/SSHFS/SOL_home/weinhol/Promotion"
source(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/KW_functions.r"))
source(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/MY_DE_Bayes.functions.r"))

#### data
  #write.table(SF,"SF_matrix.txt")
  SF<-read.table(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/SF_matrix.txt"),T) ##
  FILE.N<-"qu_bkgd"
  #SF2<-SF
  #x<-get.log.data(SF[1:100,])

		ALT_15 <- read.delim(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/qu_nobkgd.FoldChange.T=1.5.txt"))
		ALT_125 <- read.delim(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/qu_nobkgd.FoldChange.T=1.25.txt"))
		ALT_15N<-as.character(ALT_15[,1])
		ALT_125N<-as.character(ALT_125[,1])
		
    ## VGL !!! DATEI !!! 
    FC_125_BG <- read.delim(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/NEW_BG/qu_bkgd.FoldChange.T=1.25.txt"))
    FC_125_BG.N<-as.character(FC_125_BG[,1])
    FC_125_BG.10 <- read.delim(paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/NEW_BG/qu_bkgd_corr10.FoldChange.T=1.25.txt"))
    FC_125_BG.N.10<-as.character(FC_125_BG.10[,1])


  	corrr=1 #nicht auf 0 sondern 1,2,3 setzen
  	correction=c(0,10,20)
  	cor<-correction[corrr]
		#x<-get.log2.data(SF-cor)
    x<-get.ln.data(SF-cor)
    ## AIC BIC mit VAR 
		x.L.var<-get.ML.logLs.var(x)
		x.Aka.var <- make.IC(x.L.var,c(2,3),2) # AIC
		x.SBay.var<- make.IC(x.L.var,c(2,3),6) # BIC with 6 = number of values 	
		
		Aka.V.c<-get.Sig.class(x.Aka.var,3,-5)
		Aka.V.d<-get.Sig.class(x.Aka.var,4,-5)
		Bay.V.c<-get.Sig.class(x.SBay.var,3,-5)    
		Bay.V.d<-get.Sig.class(x.SBay.var,4,-5)
		## AIC BIC mit VAR
 
TAUo<-seq(0.25,2,.25)
for(to in 1:length(TAUo)){
  
    tau<-TAUo[to]	
    #tau<-1
    print(tau)
    
		pis<-c(0.90,0.008,0.07,0.022)
		#pis<-c(0.25,0.25,0.25,0.25)
		mu0<-c(0,0,0,0,0,0,0)
		gam0<- 0.2
		btest <- do.bayestest(x,pis,tau,gam0,mu0)

	#### Auswertung
		#head(btest$Z)
		
		x.C2<-btest$C2
		x.C2.maxV<- apply( x.C2,MARGIN=1,FUN=maxIndex)
		table(x.C2.maxV)
		x.C<-btest$C
		x.C.maxV<-apply(x.C,1,get.C.maxV)
		table(x.C.maxV)

		#x.C2.sig<-btestC2[which(apply(btestC2[,c(2,3,4)], 1, max)>0.95),c(2,3,4)]

		TMP.c1<-x.C2[x.C2.maxV==3,]
		if( is.matrix(TMP.c1) ){
			TMP.or<-order(TMP.c1[,3],decreasing=T)
			TMP.c<-TMP.c1[TMP.or,]
			x.C2.c<-list(SIG=TMP.c,N=as.character(rownames(TMP.c))) ## RESULT for class c !!! 
		} 
		if( is.vector(TMP.c1) ){
			x.C2.c<-list(SIG=TMP.c1,N= as.character(names(x.C2.maxV[x.C2.maxV==3])))
		}
    input<-f.input3(ALT_15N  ,ALT_125N ,x.C2.c$N,c("FC_1.5","FC_1.25","C2"))  ;venn(input) ; 
    
	###### AIC , BIC

		x.L<-get.ML.logLs(x,tau)
			x.Aka <- make.IC(x.L,c(1,2),2) # AIC
			x.SBay<- make.IC(x.L,c(1,2),6) # BIC with 6 = number of values 

		Aka.c<-get.Sig.class(x.Aka,3)
		Bay.c<-get.Sig.class(x.SBay,3)
    
		ti<-get.Time()
		out=paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/List/")

		pdf(paste(sep="",out,"PLOT_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,".pdf"))
  		hist(x,200,freq=F)
  		lines(density(x),col=2)
  		input<-f.input4(Aka.V.c$N,Bay.V.c$N,Aka.c$N,Bay.c$N,c("AIC.Var","BIC.Var","AIC","BIC")) ;venn(input) ;
  		input<-f.input4(ALT_15N,ALT_125N,Aka.V.c$N,Bay.V.c$N,c("FC_1.5","FC_1.25","AIC.Var","BIC.Var"))  ;venn(input) ;
  		input<-f.input4(ALT_15N,ALT_125N,Aka.c$N,Bay.c$N,c("FC_1.5","FC_1.25","AIC","BIC"))  ;venn(input) ;
      input<-f.input3(ALT_15N  ,ALT_125N ,x.C2.c$N,c("FC_1.5","FC_1.25","C2"))  ;venn(input) ; 
  		input<-f.input3(Aka.V.c$N,Bay.V.c$N,x.C2.c$N,c("AIC.Var","BIC.Var","C2"))  ;venn(input) ;
  		input<-f.input3(Aka.c$N  ,Bay.c$N  ,x.C2.c$N,c("AIC","BIC","C2"))  ;venn(input) ;
    
      v1<-venn(input);v2<-v1[2:8,1];names(v2)<-c("A","B","B&C","C","A&C","A&B","A&B&C")
      plot(venneuler(v2))
		dev.off()

		write.table(fout(Aka.V.c),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_AIC.v_class_c.txt"),row.names=F,col.names=T,quote=F)
		write.table(fout(Aka.V.d),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_AIC.v_class_d.txt"),row.names=F,col.names=T,quote=F)
		write.table(fout(Aka.c),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_AIC_class_c.txt"),row.names=F,col.names=T,quote=F)
		write.table(fout(Bay.c),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_BIC_class_c.txt"),row.names=F,col.names=T,quote=F)
		write.table(fout(Bay.V.c),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_BIC.v_class_c.txt"),row.names=F,col.names=T,quote=F)
		write.table(fout(Bay.V.d),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_BIC.v_class_d.txt"),row.names=F,col.names=T,quote=F)
		write.table(fout(x.C2.c),paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_C2_class_c.txt"),row.names=F,col.names=T,quote=F)
    
    
    write.table(Bay.V.c$N,paste(sep="",out,"LIST_tau_",tau,"_cor_",cor,"_",ti$d,"_",ti$s,":",ti$m,"_BIC.v_class_c_Names.txt"),row.names=F,col.names=T,quote=F)
    
    
	}

  #TESTdata 
	TAU<-seq(0.25,2,.25)
 
  for(o in 1:length(TAU)){

		tau<-TAU[o]
    #tau<-1
		print(tau)	

    #nubV=500000
     nubV=10000
    k<-get.k(nubV)
#    k<-rep(4,10)
    d<-get.Data(k,tau,gam0,mu0)
    PIS=as.vector(table(k)/length(k))
#    PIS=c(0.001,0.001,0.001,0.997)
    test.d <- do.bayestest(d,PIS,tau,gam0,mu0)
  
    # evaluate
    d.C2<-test.d$C2
    d.C2.maxV<- apply(d.C2,MARGIN=1,FUN=maxIndex)
    CM.C2<-get.ConfMat(k,d.C2.maxV)
		
    d.C<-test.d$C    
		d.C.maxV<-apply(d.C,1,get.C.maxV)
		CM.C<-get.ConfMat(k,d.C.maxV)

    L<-get.ML.logLs(d,tau)
			Aka <- make.IC(L,c(1,2),2) # AIC
			SBay<- make.IC(L,c(1,2),6) # BIC with 6 = number of values 
			CM.A<-get.ev(k,Aka)
			CM.B<-get.ev(k,SBay)				
    
		L.var<-get.ML.logLs.var(d)
		Aka <- make.IC(L.var,c(2,3),2) # AIC
		SBay<- make.IC(L.var,c(2,3),6) # BIC with 6 = number of values 		
		CM.A.v<-get.ev(k,Aka)
		CM.B.v<-get.ev(k,SBay)

		CM<-list(C2<-CM.C2,C<-CM.C,A<-CM.A,B<-CM.B,Av<-CM.A.v,Bv<-CM.B.v)	
		vgl<-matrix(0,6,5,dimnames=list(c("C2","C","AIC","BIC","AIC.VAR","BIC.VAR"),c("a","b","c","d","all")))
		for(j in 1:6){
			vgl[j,1:4]<-diag(CM[[j]])/rowSums(CM[[j]])
			vgl[j,5]<-sum(diag(CM[[j]]))/length(k)
		}
		print(vgl)
  	saveSet(d,k,vgl,tau)    
}




#################FOR HENRI BIC results ###########################
###
out=paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/Auswertung/List/")
cor2=0



X.TEST<-get.ln.data(SF-cor2)
X.TEST<-get.ln.data(SF[goodPV,]-cor2)
#X.TEST<-get.log2.data(SF-cor2)

X.TEST_L<-get.ML.logLs.var(X.TEST)
A <- make.IC(X.TEST_L,c(2,3),2) # AIC
B <- make.IC(X.TEST_L,c(2,3),6) # BIC with 6 = number of values   

A.o<-A$IC[A$minV==3,]
B.o<-B$IC[B$minV==3,]
or<-order(A.o[,3],decreasing=F)
A.o2<-A.o[or,]
A.N<-as.character(rownames(A.o2))
length(A.N)
or2<-order(B.o[,3],decreasing=F)
B.o2<-B.o[or2,]
B.N<-as.character(rownames(B.o2))
length(B.N)    

hist(A.o2[,3],1000)
hist(B.o2[,3],1000)


input<-f.input2(A.N,B.N,c("AIC.all","BIC.all"))  ;venn(input) ;
input<-f.input3(ALT_125N,FC_125_BG.N,A.N,c("FC_1.25","FC_1.25_BG","AIC.all"))  ;venn(input) ;
input<-f.input3(ALT_125N,FC_125_BG.N,B.N,c("FC_1.25","FC_1.25_BG","BIC.all"))  ;venn(input) ;

SIGV<- (-5)
A.c.SIG<-get.Sig.class(A,3,SIGV)
length(A.c.SIG$N)

B.a.SIG<-get.Sig.class(B,1,SIGV)
length(B.a.SIG$N)
B.b.SIG<-get.Sig.class(B,2,SIGV)
length(B.b.SIG$N)
B.c.SIG<-get.Sig.class(B,3,SIGV)
length(B.c.SIG$N)
B.d.SIG<-get.Sig.class(B,4,SIGV)
length(B.d.SIG$N)

input<-f.input2(ALT_125N,B.c.SIG$N,c("FC_1.25","BIC.var"))  ;venn(input) ;
input<-f.input2(FC_125_BG.N,B.c.SIG$N,c("FC_1.25_BG","BIC.var"))  ;venn(input) ;    
input<-f.input3(ALT_125N,FC_125_BG.N,B.c.SIG$N,c("FC_1.25","FC_1.25_BG","BIC.var"))  ;venn(input) ;

input<-f.input2(ALT_125N,A.c.SIG$N,c("FC_1.25","AIC.var"))  ;venn(input) ;
input<-f.input2(FC_125_BG.N,A.c.SIG$N,c("FC_1.25_BG","AIC.var"))  ;venn(input) ;    
input<-f.input3(ALT_125N,FC_125_BG.N,A.c.SIG$N,c("FC_1.25","FC_1.25_BG","AIC.var"))  ;venn(input) ;
input<-f.input2(A.c.SIG$N,B.c.SIG$N,c("AIC.var","BIC.var"))  ;venn(input) ;

input<-f.input4(B.a.SIG$N,B.b.SIG$N,B.c.SIG$N,B.d.SIG$N,c("BIC.a","BIC.b","BIC.c","BIC.d"))  ;venn(input) ;

## FOR HENRI !!!!

write.table(fout(B.c.SIG),paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_c.txt"),row.names=F,col.names=T,quote=F,sep="\t")    
write.table(B.c.SIG$N,paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_c_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(fout(B.d.SIG),paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_d.txt"),row.names=F,col.names=T,quote=F,sep="\t")    
write.table(B.c.SIG$N,paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_d_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")

write.table(fout_2(B.a.SIG),paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_a.txt"),row.names=F,col.names=T,quote=F,sep="\t")    
write.table(B.a.SIG$N,paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_a_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(fout_2(B.b.SIG),paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_b.txt"),row.names=F,col.names=T,quote=F,sep="\t")    
write.table(B.b.SIG$N,paste(sep="",out,"qu_bkgd_cor_",cor2,"_BIC.V_sig_",SIGV,"_class_b_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")



pdf(paste(sep="",out,"VENN_foldChange1.25_vs_BIC.V.pdf"))
input<-f.input2(FC_125_BG.N,B.c.SIG$N,c("FC_1.25_BG","BIC.var(-5)_corr10"))  ;venn(input) ;    
dev.off()
write.table(setdiff(FC_125_BG.N,B.c.SIG$N),paste(sep="",out,"VENN_foldChange1.25_vs_BIC.V__FC_only_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(setdiff(B.c.SIG$N,FC_125_BG.N),paste(sep="",out,"VENN_foldChange1.25_vs_BIC.V__BIC_only_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(intersect(B.c.SIG$N,FC_125_BG.N),paste(sep="",out,"VENN_foldChange1.25_vs_BIC.V__BIC_and_FC_NAMES.txt"),row.names=F,col.names=T,quote=F,sep="\t")
## FOR HENRI !!!!


#B.c.SIG.10<-B.c.SIG

pdf(paste(sep="",out,"VENN_foldChange1.25_vs_BIC.V_4er.pdf"),11,11)
input<-f.input4(FC_125_BG.N,FC_125_BG.N.10,B.c.SIG.10$N,B.c.SIG$N,c("FC_1.25_BG_corr0","FC_1.25_BG_corr10","BIC.var(-5)_corr10","BIC.var(-5)_corr0"))  ;venn(input) ;venn(input,simplify=T);    

input<-f.input2(FC_125_BG.N,FC_125_BG.N.10 ,c("FC_1.25_BG_corr0","FC_1.25_BG_corr10"))  ;venn(input) ;    
input<-f.input2(B.c.SIG.10$N,B.c.SIG$N     ,c("BIC.var(-5)_corr10","BIC.var(-5)_corr0"))  ;venn(input) ;    
input<-f.input2(FC_125_BG.N,B.c.SIG$N      ,c("FC_1.25_BG_corr0","BIC.var(-5)_corr0"))  ;venn(input) ;    
input<-f.input2(FC_125_BG.N.10,B.c.SIG.10$N,c("FC_1.25_BG_corr10","BIC.var(-5)_corr10"))  ;venn(input) ;    

dev.off()

 #########################
# Posterior P(m|x ) aus BIC berechnen 
#-> BIC=-2 ln(P(x|m)) =>  P(x|m)= exp(-BIC/2)       => function q
# P(m|x) = (P(x|m) * P(m) [bzw pi]) / sum_m Zähler  => function w
  qB<-function(B){
    return(exp(-B/2))
  }
  
  P_x_m<-apply(B$IC,MARGIN=1,FUN=qB)
  P_x_m<-t(log(P_x_m))
  
  #
  pis<-c(0.90,0.008,0.07,0.022)
  #pis<-c(0.25,0.25,0.25,0.25)
  
  wB<-function(B,pis){
    t<-c()
    for(i in 1:4){
      t[i]<-B[i]*pis[i]
    }
    return(t/sum(t))
  }
  
  P_m_x<-apply(exp(P_x_m),MARGIN=1,FUN=wB,pis   )
  P_m_x<-t(P_m_x)
  
  #m1<-apply((P_x_m),MARGIN=1,FUN=maxIndex)
  m2<-apply((P_m_x),MARGIN=1,FUN=maxIndex)
  #m3<-apply(X.TEST_L,MARGIN=1,FUN=maxIndex)
  #m4<-apply(B$IC,MARGIN=1,FUN=minIndex)
  
  #print(table(m1))
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
  
zz<-P_m_x_c
input=f.input2( rownames(P_m_x_c),rownames(zz),c("Paper","equal"))
venn(input,simplify=F)
input=f.input2( rownames(head(P_m_x_c,2000)),rownames(zz),c("Paper","equal"))
venn(input,simplify=F)

tID<-intersect(rownames(B$IC),EGFR_genes.IDs[,1])
head(B$IC[tID,])
head(X.TEST_L[tID,])
table(apply((X.TEST_L),MARGIN=1,FUN=maxIndex))
table(apply((B$IC),MARGIN=1,FUN=minIndex))

  write.table(fout_3(P_m_x),paste(sep="",out,"qu_bkgd_BIC.V_class_abcd_mit_Gewichten.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  write.table(fout_3(P_m_x_a),paste(sep="",out,"qu_bkgd_BIC.V_class_a_mit_Gewichten.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  write.table(fout_3(P_m_x_b),paste(sep="",out,"qu_bkgd_BIC.V_class_b_mit_Gewichten.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  write.table(fout_3(P_m_x_c),paste(sep="",out,"qu_bkgd_BIC.V_class_c_mit_Gewichten.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  write.table(fout_3(P_m_x_d),paste(sep="",out,"qu_bkgd_BIC.V_class_d_mit_Gewichten.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  write.table(rownames(P_m_x_a),paste(sep="",out,"qu_bkgd_BIC.V_class_a_mit_Gewichten_N.txt"),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(rownames(P_m_x_b),paste(sep="",out,"qu_bkgd_BIC.V_class_b_mit_Gewichten_N.txt"),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(rownames(P_m_x_c),paste(sep="",out,"qu_bkgd_BIC.V_class_c_mit_Gewichten_N.txt"),row.names=F,col.names=F,quote=F,sep="\t")
  write.table(rownames(P_m_x_d),paste(sep="",out,"qu_bkgd_BIC.V_class_d_mit_Gewichten_N.txt"),row.names=F,col.names=F,quote=F,sep="\t")

qu.bg <- read.delim("~/Promotion/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_quantile_bkgd_SF767-1-799-6.txt",sep="\t")
FILE<-qu.bg
FILE.N<-"qu_bkgd_corr10"

rownames(FILE)<-FILE[,1] 


write.table(FILE[rownames(P_m_x_a),c(1,2,52,53,54,55)],paste(sep="",out,"qu_bkgd_BIC.V_class_a_mit_Gewichten_Gen_Informations.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(FILE[rownames(P_m_x_b),c(1,2,52,53,54,55)],paste(sep="",out,"qu_bkgd_BIC.V_class_b_mit_Gewichten_Gen_Informations.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(FILE[rownames(P_m_x_c),c(1,2,52,53,54,55)],paste(sep="",out,"qu_bkgd_BIC.V_class_c_mit_Gewichten_Gen_Informations.txt"),row.names=F,col.names=T,quote=F,sep="\t")
write.table(FILE[rownames(P_m_x_d),c(1,2,52,53,54,55)],paste(sep="",out,"qu_bkgd_BIC.V_class_d_mit_Gewichten_Gen_Informations.txt"),row.names=F,col.names=T,quote=F,sep="\t")


  ## pVal Filter
  S<-0.01
  goodPV<-(rowSums(cbind(SF.pval[,1]<S,SF.pval[,2]<S,SF.pval[,3]<S,SF.pval[,4]<S,SF.pval[,5]<S,SF.pval[,6]<S))==6)
  RNsPV<-rownames(SF.pval)
  length(RNsPV[goodPV])
  RNs.sig<-intersect(RNs,RNsPV[goodPV])
  
  SF.pval["ILMN_2212999",] 

  ## ADD FC
  a0<-c(1,2,3,4,5,6);   b0<-c(1,3,5);       b1<-c(2,4,6);         c0<-c(1,3,5,4);      c1<-c(2,6);         d0<-c(1,3,5,4,6);    d1<-c(2);            
  qu.bg <- read.delim( paste(sep="",FOLDER,"/Kappler_Wichmann_Medizin/KW_120813_1343/Analysen_SF767-1-799-6/Einzelanalyse_quantile_bkgd_SF767-1-799-6.txt"),sep="\t")
  FILE<-qu.bg
  FILE.N<-"qu_bkgd"
  IDs<-FILE[,c(1:2,54,55)]
  rownames(IDs)<-as.character(FILE[,1])
  head(IDs)
  EGFR_genes<-function(){
    EGFR_genes<-c("ABL1","ABL2","ACTA1","ADRBK1","ADRM1","AGTR1","AHNAK","AHSA1","AIP","AKAP12","ALCAM","ALDOA","AMH","ANKS1A","ANKS1B","ANXA1","APBA3","APBB1","APBB2","APBB3","APPL1","AR","AREG","ARF4","ARF6","ATIC","ATP1A1","ATP1B1","ATP5C1","BLK","BMX","BTC","CALM1","CALM2","CALM3","CAMK2A","CAMK2G","CAMLG","CASP1","CAV1","CAV2","CAV3","CBL","CBLB","CBLC","CBR1","CD44","CD59","CD82","CDC25A","CDH1","CEACAM1","CEBPB","CEND1","CFL1","CISH","CLNK","CLTA","CLTCL1","CMTM8","CNTN2","COL9A3","COX2","CRK","CSRP1","CTNNB1","CTNND1","DCN","DCTN2","DEGS1","DNAJC4","DOK2","DOK4","DOK5","DOK6","DUSP3","DYNC1H1","EEF1G","EGF","ELF3","EPB41","EPHA2","EPPK1","EPS15","EPS8","ERBB2","ERBB3","ERBB4","EREG","ERRFI1","ESR1","EZR","FAH","FAS","FBXO25","FER","FES","FKBP4","FKBP8","FN1","GAB1","GABARAPL2","GAPDH","GNAI2","GOT1","GPM6B","GRB10","GRB14","GRB2","GRB7","HBEGF","HEXIM1","HGS","HIST3H3","HOXC10","HSP90B1","HSPA1A","HSPA4","HSPA8","HTT","ICAM1","INPPL1","IRS1","IRS4","ITGA5","JAK2","JUP","KCTD9P2","KRT17","KRT18","KRT7","KRT8","LYN","MAP2K1","MAP3K12","MAP3K14","MAP4K1","MAPK1","MAPK14","MAPK8IP1","MAPK8IP2","MAST1","MATR3","MDH1","MEGF6","MET","MGARP","MOB4","MUC1","NCAM1","NCK1","NCK2","NDN","NEDD4","NUMB","NUMBL","OAZ1","OLFM1","PAK1","PCNA","PDCD6IP","PDGFRB","PDK3","PIK3C2A","PIK3C2B","PIK3R1","PIK3R2","PIK3R3","PIN4","PITPNA","PKIA","PLCG1","PLCG2","PLD1","PLD2","PLEC","PLSCR1","PPP5C","PRCC","PRDX1","PRKACA","PRKAR1A","PRKAR1B","PRKCA","PRKCB","PRKCD","PRKD1","PSMA7","PSMD4","PTGDS","PTK2","PTK2B","PTK6","PTPN1","PTPN11","PTPN12","PTPN2","PTPN6","PTPRJ","PTPRS","RAB3A","RAP1GDS1","RASA1","RGS16","RIN2","RIPK1","RNF115","RNF126","RQCD1","RUSC2","S100A7","S100A9","SCAMP1","SCAMP3","SEC13","SEPP1","SFN","SGSM2","SH2B1","SH2B3","SH2D1A","SH2D2A","SH2D3A","SH3BGRL3","SHC1","SHC2","SHC3","SLA","SLC3A2","SLC9A3R1","SMURF2","SNRPD2","SNX1","SNX2","SNX4","SNX6","SOCS1","SOCS3","SOCS4","SOCS5","SOCS7","SORBS2","SOS1","SOS2","SPARCL1","SPCS2","SRC","STAM2","STAT1","STAT3","STAT5A","STAT5B","STIP1","STUB1","SYK","TGFA","TJP1","TLN1","TMCO3","TNC","TNK2","TNS4","TPI1","TPM1","TRPV1","TUBA1A","TUBA4A","UBB","UBE2V2","UCHL1","UROD","VAPA","VAV1","VAV2","VAV3","XRCC6","YWHAZ","ZAP70","ZNF259")
    EGFR_genes.IDs<-IDs[IDs[,2]==EGFR_genes[1],]
    for(i in 2:length(EGFR_genes)) EGFR_genes.IDs<-rbind(EGFR_genes.IDs,IDs[IDs[,2]==EGFR_genes[i],])
    P_m_x[intersect(rownames(P_m_x),EGFR_genes.IDs[,1]),]
    EGFR_genes.IDs_P_m_x<- P_m_x[intersect(rownames(P_m_x),EGFR_genes.IDs[,1]),]
    EGFR_genes.IDs_P_m_x[EGFR_genes.IDs_P_m_x[,3]>0.90,]
    EGFR_genes.IDs[rownames(EGFR_genes.IDs_P_m_x[EGFR_genes.IDs_P_m_x[,3]>0.95,]),]
    IDs[RNs,][202,]
    P_m_x_c2<-P_m_x_c[P_m_x_c[,3]>0.95,]
    RNs<-rownames(P_m_x_c2)
  }
  #plot for number of genes at which prob
  prob_plot<-function(){
    sum(P_m_x_c[,3]>0.999999)
    qq<-c()
    for(i in seq(from=0.3,to=1,by=.05))qq<-c(qq,sum(P_m_x_c[,3]>i) )
    names(qq)<-as.character(seq(from=0.3,to=1,by=.05) ) 
    barplot(sort(qq,decreasing=T),las=2)
    plot(seq(from=0.3,to=1,by=.05)*100,qq)
  
    sum(P_m_x_d[,4]>0.999)
    qq<-c()
    for(i in seq(from=0.3,to=1,by=.05))qq<-c(qq,sum(P_m_x_d[,4]>i) )
    names(qq)<-as.character(seq(from=0.3,to=1,by=.05) ) 
    barplot(sort(qq,decreasing=T),las=2)
    plot(seq(from=0.3,to=1,by=.05)*100,qq)
    
    P_m_x_d2<-P_m_x_d[P_m_x_d[,4]>0.95,]
    RNs<-rownames(P_m_x_d2)
    FC<-exp((X.TEST[RNs,d1])-rowMeans(X.TEST[RNs,d0]))
    P_m_x_d_FC<-cbind(P_m_x_d2,FC,log2(FC),exp(rowMeans(X.TEST[RNs,])))
    P_m_x_d_FC<-P_m_x_d_FC[order(abs(P_m_x_d_FC[,6]),decreasing=T),]
    
    RNsFC<-rownames(P_m_x_d_FC)
    GrD<-cbind(P_m_x_d_FC[RNsFC,c(3)],SF[RNsFC,],P_m_x_d_FC[RNsFC,c(7,5,6)])
    colnames(GrD)<-c("prob_for_D",colnames(SF),"geoMeanExp","FC","log2FC")  
    B.tmp<-as.matrix(signif(GrD,3))
    B<-as.matrix(cbind(rownames(GrD), as.character(IDs[rownames(GrD),2]), B.tmp))
    colnames(B)[c(1,2)]<-c("ID","Sym")
    head(B)
    write.table(B,paste(sep="",out,"NEW_qu_bkgd_BIC_grD_AllInfo_arrayPval_0.01.txt"),row.names=F,col.names=T,quote=F,sep="\t")
    
  }
  ##log
  FC<-exp(rowMeans(X.TEST[RNs,c1])-rowMeans(X.TEST[RNs,c0]))
  #  FC<-exp(rowMeans(X.TEST[RNs,c1])-rowMeans(X.TEST[RNs,c0]))
  P_m_x_c_FC<-cbind(P_m_x_c2,FC,log2(FC),exp(rowMeans(X.TEST[RNs,])))
  P_m_x_c_FC<-P_m_x_c_FC[order(abs(P_m_x_c_FC[,6]),decreasing=T),]
  #log2
  #FC<-2^(rowMeans(X.TEST[RNs,c1])-rowMeans(X.TEST[RNs,c0]))
  #P_m_x_c_FC<-cbind(P_m_x_c2,FC,log2(FC),2^(rowMeans(X.TEST[RNs,])))
  #P_m_x_c_FC<-P_m_x_c_FC[order(abs(P_m_x_c_FC[,6]),decreasing=T),]
  
  ## RNs.sig pVal-Filter !!!!
    FC<-exp(rowMeans(X.TEST[RNs.sig,c1])-rowMeans(X.TEST[RNs.sig,c0]))
    P_m_x_c_FC<-cbind(P_m_x_c2[RNs.sig,],FC,log2(FC),exp(rowMeans(X.TEST[RNs.sig,])))
    P_m_x_c_FC<-P_m_x_c_FC[order(abs(P_m_x_c_FC[,6]),decreasing=T),]
  ####
  
  RNsFC<-rownames(P_m_x_c_FC)
  GrC<-cbind(P_m_x_c_FC[RNsFC,c(3)],SF[RNsFC,],P_m_x_c_FC[RNsFC,c(7,5,6)])
  colnames(GrC)<-c("prob_for_C",colnames(SF),"geoMeanExp","FC","log2FC")  
  B.tmp<-as.matrix(signif(GrC,3))
  B<-as.matrix(cbind(rownames(GrC), as.character(IDs[rownames(GrC),2]), B.tmp))
  colnames(B)[c(1,2)]<-c("ID","Sym")
  head(B)
  write.table(B,paste(sep="",out,"NEW_qu_bkgd_BIC_grC_AllInfo.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  #with RNs.sig
  write.table(B,paste(sep="",out,"NEW_qu_bkgd_BIC_grC_AllInfo_arrayPval_0.01.txt"),row.names=F,col.names=T,quote=F,sep="\t")


  write.table(fout_4(P_m_x_c_FC,ID),paste(sep="",out,"NEW_qu_bkgd_BIC_grC.txt"),row.names=F,col.names=T,quote=F,sep="\t")
  pdf(paste(sep="",out,"NEW_qu_bkgd_BIC_grC_expressions.pdf"),15,10)
  RNsFC<-rownames(P_m_x_c_FC)
  for(i in 1:length(RNsFC)){
  if(length(setdiff(RNsFC[i],rownames(X))==RNsFC[i])==1){
    Merge<-as.matrix(rbind(SF[RNsFC[i],],c(0,0,0,0,0,0)))
  }else{
    Merge<-as.matrix(rbind(SF[RNsFC[i],],X[RNsFC[i],]))
  }
   
  both<-"b";
    if(both=="a"){
    COLs<-c(brewer.pal(12,"Paired")[c(2,1)],brewer.pal(12,"Paired")[c(4,3)],brewer.pal(12,"Paired")[c(2,1)],
            brewer.pal(12,"Paired")[c(2,1)],brewer.pal(12,"Paired")[c(2,1)],brewer.pal(12,"Paired")[c(4,3)])
    barplot(Merge,beside=T,main=IDs[rownames(Merge)[1],2],names.arg=c("KO","KO+EGF","ALL+EGF","ALL","14","14+EGF"),legend.text=c("SF (C-)","799(C-)","SF (C+)","799 (C+)"),las=2,col=COLs)
    }else{
    COLs<-c(brewer.pal(12,"Paired")[c(2)],brewer.pal(12,"Paired")[c(4)],brewer.pal(12,"Paired")[c(2)],
            brewer.pal(12,"Paired")[c(2)],brewer.pal(12,"Paired")[c(2)],brewer.pal(12,"Paired")[c(4)])
    barplot(Merge[1,],beside=T,main=paste(sep=" -> ",rownames(Merge)[1],IDs[rownames(Merge)[1],2]),names.arg=c("KO","KO+EGF","ALL+EGF","ALL","14","14+EGF"),legend.text=c("SF (C-)","SF (C+)"),las=2,col=COLs)
    }
  }
  dev.off()
  
  #save.image(paste(sep="",out,"AddFC.RData")) # 07.03.2014 # ACHTUNG Mit log2 gerechnete Werte !!!! 
                    

#################QPCR ################################
RawQPCR <- read.csv("~/Promotion/Kappler_Wichmann_Medizin/Auswertung/qPCR/MyData.csv")
qgenes<-c("TPR","CKAP2L","KIF5C","ROCK1")
qgenesIDs<-c()
for(i in 1:4){
qgenesIDs[[i]]<-as.character(IDs[which(qgenes[i]==IDs[,"SYMBOL"]),1])
}
names(qgenesIDs)<-qgenes

qGenes<-c()
qGenes[[1]] <-    cbind(RawQPCR[1:6,c(3,4)],RawQPCR[13:18,c(3,4)],RawQPCR[25:30,c(3,4)])
qGenes[[2]]<-  cbind(RawQPCR[1:6,c(5,6)],RawQPCR[13:18,c(5,6)],RawQPCR[25:30,c(5,6)])
qGenes[[3]]<-   cbind(RawQPCR[1:6,c(7,8)],RawQPCR[13:18,c(7,8)],RawQPCR[25:30,c(7,8)])
qGenes[[4]]<-   cbind(RawQPCR[1:6,c(9,10)],RawQPCR[13:18,c(9,10)],RawQPCR[25:30,c(9,10)])
names(qGenes)<-qgenes

qGenes.M<-c()
for(i in 1:4){
qGenes.M[[i]]<-cbind(rowMeans(log2(qGenes[[i]][c(1,3,5)])),rowMeans(log2(qGenes[[i]][c(2,4,6)])),  log2(rowMeans(qGenes[[i]][c(2,4,6)])) )
colnames(qGenes.M[[i]])<-c("log2geoMeanCT","log2geoMeanREL","log2MeanREL")
}
names(qGenes.M)<-qgenes

qGenes.FC<-c()
for(i in 1:4){
  qGenes.FC[[i]]<-colMeans(qGenes.M[[i]][c1,])- colMeans(qGenes.M[[i]][c0,])
}

log2(FC[c(qgenesIDs[[1]])])
log2(FC[c(qgenesIDs[[2]])])
log2(FC[c(qgenesIDs[[3]])])
log2(FC[c(qgenesIDs[[4]])])

qGenes.FC.comp<-c()
qGenes.FC.comp<-c(qGenes.FC[[1]][c(2,3)],log2(FC[c(qgenesIDs[[1]])[2] ]))
qGenes.FC.comp<-rbind(qGenes.FC.comp,
                      c(qGenes.FC[[2]][c(2,3)],log2(FC[c(qgenesIDs[[2]]) ])) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,
                      c(qGenes.FC[[4]][c(2,3)],log2(FC[c(qgenesIDs[[4]])[2] ])) )
qGenes.FC.comp<-rbind(qGenes.FC.comp,
                      c(qGenes.FC[[3]][c(2,3)],log2(FC[c(qgenesIDs[[3]])[1] ])) )
rownames(qGenes.FC.comp)<-qgenes[c(1,2,4,3)]

barplot(t(qGenes.FC.comp[,c(1,3)]),beside=T,las=2,legend.text=c("qPCR","Array"),main="log2Foldchange for class C ")

################################################




if(0){
  
  ### IVO test fr.22.2
  hist(log(SF2-80),200)
  y<-get.log.data( (SF-10))
  qqnorm(y)
  l<-length(y[,1])
  x6<-apply(x,1,sd)    
  mean(x6)
  y6<-apply(y,1,sd)    
  mean(y6)
  1/(median(y6))^2  
  y1<-5/6*apply(y,1,var)
  hist(y1)
  y2<-rowMeans(y)
  y3<-(l-1)/l*var(y2) 
  y4<-(6*l-1)/6/l*var(as.vector(y))  
  y4-y3-mean(y1)
  hist(y1,220)
  ### IVO test fr.22.2
  
  #TEST1
  x.L.var<-get.ML.logLs.var(x)
  x.SBay.var<- make.IC(x.L.var,c(2,3),6) # BIC with 6 = number of values
  Bay.V.c<-get.Sig.class(x.SBay.var,3,-5)
  
  xl2.L.var<-get.ML.logLs.var(xl2)
  xl2.SBay.var<- make.IC(xl2.L.var,c(2,3),6) # BIC with 6 = number of values
  l2.Bay.V.c<-get.Sig.class(xl2.SBay.var,3,-5)
  
  hist(xl2.SBay.var$IC[xl2.SBay.var$minV==3,3])
  hist(x.SBay.var$IC[x.SBay.var$minV==3,3])
  plot(density(xl2.SBay.var$IC[xl2.SBay.var$minV==3,3]))
  lines(density(x.SBay.var$IC[x.SBay.var$minV==3,3]))
  input<-f.input2(Bay.V.c$N ,l2.Bay.V.c$N,c("BIC.ln","BIC.log2"))  ;venn(input) ;
  
  #TEST2
  x<-get.ln.data(SF-0)
  ## AIC BIC mit VAR 
  x.L.var<-get.ML.logLs.var(x)
  x.SBay.var<- make.IC(x.L.var,c(2,3),6) # BIC with 6 = number of values
  Bay.V.c<-get.Sig.class(x.SBay.var,3,-5)
  
  x1<-get.ln.data(SF-10)
  ## AIC BIC mit VAR 
  x1.L.var<-get.ML.logLs.var(x1)
  x1.SBay.var<- make.IC(x1.L.var,c(2,3),6) # BIC with 6 = number of values
  Bay.V.c1<-get.Sig.class(x1.SBay.var,3,-5)
  
  x2<-get.ln.data(SF-20)
  ## AIC BIC mit VAR 
  x2.L.var<-get.ML.logLs.var(x2)
  x2.SBay.var<- make.IC(x2.L.var,c(2,3),6) # BIC with 6 = number of values
  Bay.V.c2<-get.Sig.class(x2.SBay.var,3,-5)    
  input<-f.input3(Bay.V.c$N ,Bay.V.c1$N,Bay.V.c2$N,c("BIC.0","BIC.1","BIC.2"))  ;venn(input) ;
  
  t<-as.character(names(x.SBay.var$IC[x.SBay.var$minV==3,3]))
  t1<-as.character(names(x1.SBay.var$IC[x1.SBay.var$minV==3,3]))
  t2<-as.character(names(x2.SBay.var$IC[x2.SBay.var$minV==3,3]))
  input<-f.input3(t ,t1,t2,c("BIC.0","BIC.1","BIC.2"))  ;venn(input) ;
  
  x.SBay.var$IC[setdiff(t1, t),]
  x.SBay.var$minV[setdiff(t1, t)]   
  x[setdiff(t1, t)[3],]
  x1[setdiff(t1, t)[3],]
  x2[setdiff(t1, t)[3],]
}


