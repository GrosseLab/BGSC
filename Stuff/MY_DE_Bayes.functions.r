library("limma")
library("VennDiagram")
#library("venneuler")
library(hexbin)
#library(ggplot2)
library(gplots)
#install.packages("ggdendro")
#library(ggdendro)
#install.packages("rpart")
library(rpart)
library(hexbin)
#library(calibrate)
library(MASS)
library(RColorBrewer)

#source("http://bioconductor.org/biocLite.R")  
#biocLite("DESeq")
#biocLite("edgeR")
#biocLite("baySeq")
#biocLite("DEXSeq")
#biocLite("vsn")
#install.packages("VennDiagram")
#install.packages("venneuler")
#install.packages("hexbin")
#install.packages("ggplot2")
#install.packages("gplots")
#install.packages("ggdendro")
#install.packages("rpart")
#install.packages("hexbin")
#install.packages("calibrate")
#install.packages("MASS")

#### functions
  do.bayestest <- function(d,pis,tau,gam0,mu0){
    # bstat <- matrix(NA,nrow=nrow(d),ncol=2)
    
    #####   functions for , taux, mux , Bx 
    taux <- function(N,tau,gam0){
      return(tau*(N+gam0))
    }
    mux <- function(N,x,gam0,mu0){
      return((N*mean(x)+gam0*mu0)/(N+gam0))
    }
    Ax <- function(N,x,gam0,mu0){
      return((sum(x^2)+gam0*mu0^2)/(N+gam0))
    }
    Bx <- function(N,x,tau,gam0,mu0){
      return(sqrt(gam0*tau/taux(N,tau,gam0))*(tau/(2*pi))^(N/2)*exp(-0.5*taux(N,tau,gam0)*(Ax(N,x,gam0,mu0)-mux(N,x,gam0,mu0)^2)))
    }
    
    #####   calc Bx
    f.B<-function(x,sub,N,tau,gam0,mu0,a){
      #a == vector mit Auswahl der Gene
      return( Bx(N[a],x[sub[[a]]],tau,gam0,mu0[a])) 
    }
    
    N<-c()
    a0<-c(1,2,3,4,5,6);   N[1] <- length(a0);
    b0<-c(1,3,5);         N[2] <- length(b0);
    b1<-c(2,4,6);         N[3] <- length(b1);
    c0<-c(1,3,5,4);       N[4] <- length(c0);
    c1<-c(2,6);           N[5] <- length(c1);
    d0<-c(1,3,5,4,6);     N[6] <- length(d0);
    d1<-c(2);             N[7] <- length(d1);
    sub<-list(a0=a0,b0=b0,b1=b1,c0=c0,c1=c1,d0=d0,d1=d1)
    
    B.ALL<-apply(d,MARGIN=1,FUN=f.B ,sub,N,tau,gam0,mu0,1)
    for(p in 2:length(sub)){
      TMP<-apply(d,MARGIN=1,FUN=f.B ,sub,N,tau,gam0,mu0,p)
      B.ALL<-cbind(B.ALL,TMP)    
    }
    colnames(B.ALL)<-c( "B","B.b0","B.b1","B.c0","B.c1","B.d0","B.d1")
    #print(head(B.ALL))
    
    #####   calc C 
    # aus pdf
    # C := ln (pi2/pi1)+ ln( (Bx*By) / B)  
   
    f.C<-function(B,a,pi2,pi1){
      return(log(pi2/pi1)+ log( (B[a[2]]*B[a[3]])/B[a[1]]) )
    }
    C.b<-apply(B.ALL,MARGIN=1,FUN=f.C,c(1,2,3),pis[2],pis[1])
    C.c<-apply(B.ALL,MARGIN=1,FUN=f.C,c(1,4,5),pis[3],pis[1])
    C.d<-apply(B.ALL,MARGIN=1,FUN=f.C,c(1,6,7),pis[4],pis[1])
    C<-cbind(C.b,C.c,C.d)
    #print(head(C))
   
    ## Z-test
    f.Z<-function(x,tau,N,sub,a0,a1){
      Mean.a0<-mean(x[sub[[a0]]])
      Mean.a1<-mean(x[sub[[a1]]])
      z<-sqrt( (N[a0]*N[a1]) /(N[a0]+N[a1]) ) * sqrt(tau) * (Mean.a0-Mean.a1)
      return(z)  
    }
    
    Z.b<-apply(d,1,f.Z,tau,N,sub,2,3)
    Z.c<-apply(d,1,f.Z,tau,N,sub,4,5)
    Z.d<-apply(d,1,f.Z,tau,N,sub,6,7)
    Z<-cbind(Z.b,Z.c,Z.d)
    #print(head(Z))
    
    ### calc "C2"               numerator   / denominator
    # P(e \in {a,b,c,d } | x ) = B_e * pi_e / sum_(k \in {a,b,c,d})  B_k * pi_k  
    f.C2<-function(B.ALL,pis){
      numerators = c(B.ALL[1]         *pis[1],
                     B.ALL[2]*B.ALL[3]*pis[2],
                     B.ALL[4]*B.ALL[5]*pis[3], 
                     B.ALL[6]*B.ALL[7]*pis[4])
      denominator=sum(numerators)
      return(numerators/denominator)
    }
    C2<-apply(B.ALL,MARGIN=1,FUN=f.C2,pis)
    C2<-restructur(C2)
    #print(head(C2))            
    
    ##### OUT  
    return( list(B.ALL=B.ALL,C=C,C2=C2,Z=Z) )
  } 
  
  maxIndex<-function(x){
    return(which(x==max(x)))
  }
  minIndex<-function(x){
    return(which(x==min(x)))
  }


remove.neg.val<-function(x){
  if(min(x)<=0){
    return(FALSE)
  }else{
    return(TRUE)
  }
}

get.pos.data<-function(x){
  
  xr<-apply(x,1,remove.neg.val)
  print(sum(xr))
  x<-x[xr,]
  x<-as.matrix(x)
  return(x)
}


  get.log2.data<-function(x){
    remove.neg.val<-function(x){
      if(min(x)<=0){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
    
    xr<-apply(x,1,remove.neg.val)
    print(sum(xr))
    x<-x[xr,]
    x<-log2(x)
    x<-as.matrix(x)
    return(x)
  }

  get.ln.data<-function(x){
    remove.neg.val<-function(x){
      if(min(x)<=0){
        return(FALSE)
      }else{
        return(TRUE)
      }
    }
    
    xr<-apply(x,1,remove.neg.val)
    print(sum(xr))
    x<-x[xr,]
    x<-log(x)
    x<-as.matrix(x)
    return(x)
  }

###### AIC , BIC

calc.MUs<-function(d,sub){  
  MUs<-c()
  for(i in 1:length(sub)){
    MUs[i]<-mean(d[sub[[i]]])
  }
  return(MUs)
}

get.ML.logLs<-function(x,tau){
  
  NV<-function(d,mu,tau){
    #tau is precision
    return(dnorm(d,mu,sqrt(1/tau)))
    #return((tau/(2*pi))^(1/2)*exp(-0.5*tau*(x-mu)^2)) 
  }
  
  logLs<-function(d,sub,tau){
    MUs<-calc.MUs(d,sub)
    
    logL<-c()   
    logL[1]<-sum( log(NV(d[sub[[1]]],MUs[1],tau)) )
    logL[2]<-sum( log(NV(d[sub[[2]]],MUs[2],tau)) , log(NV(d[sub[[3]]],MUs[3],tau)) ) 
    logL[3]<-sum( log(NV(d[sub[[4]]],MUs[4],tau)) , log(NV(d[sub[[5]]],MUs[5],tau)) )  
    logL[4]<-sum( log(NV(d[sub[[6]]],MUs[6],tau)) , log(NV(d[sub[[7]]],MUs[7],tau)) )
    
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
  
  logL<-apply(x,1,logLs,sub,tau)
  logL<-restructur(logL)
  #print(head(logL))
  
  return(logL)
}
  
get.ML.logLs.var<-function(x){  # with var-estimator
       
  NV<-function(d,mu,tau){
      #tau is precision
      return(dnorm(d,mu,sqrt(1/tau))) # f(x) = 1/(sqrt(2 pi) sigma) e^-((x - mu)^2/(2 sigma^2))  
      #                                     return((tau/(2*pi))^(1/2)*exp(-0.5*tau*(x-mu)^2)) 
  }
    
  get.var<-function(d,sub,MUs,N,i0,i1){
	  return(	( sum((d[sub[[i0]]]-MUs[i0])^2) + sum((d[sub[[i1]]]-MUs[i1])^2)) / ( N[i0]+N[i1]-2 ) )
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
	
    logL<-c()   
    logL[1]<-sum( log(NV(d[sub[[1]]],MUs[1],TAUs[1])) )
    logL[2]<-sum( log(NV(d[sub[[2]]],MUs[2],TAUs[2])) , log(NV(d[sub[[3]]],MUs[3],TAUs[2])) ) 
    logL[3]<-sum( log(NV(d[sub[[4]]],MUs[4],TAUs[3])) , log(NV(d[sub[[5]]],MUs[5],TAUs[3])) )  
    logL[4]<-sum( log(NV(d[sub[[6]]],MUs[6],TAUs[4])) , log(NV(d[sub[[7]]],MUs[7],TAUs[4])) )  
    
    #logL2<-c()
    #logL2[1]<-sum( log(dnorm(d[sub[[1]]],MUs[1],SDs[1])) )
    #logL2[2]<-sum( log(dnorm(d[sub[[2]]],MUs[2],SDs[2])) , log(dnorm(d[sub[[3]]],MUs[3],SDs[2])) ) 
    #logL2[3]<-sum( log(dnorm(d[sub[[4]]],MUs[4],SDs[3])) , log(dnorm(d[sub[[5]]],MUs[5],SDs[3])) )  
    #logL2[4]<-sum( log(dnorm(d[sub[[6]]],MUs[6],SDs[4])) , log(dnorm(d[sub[[7]]],MUs[7],SDs[4])) )
    #logL2 == logL
    
    #print(MUs)
    #print(TAUs)
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


#xr<-x[c("ILMN_1713088","ILMN_1661717"),]
#LLL<-get.ML.logLs.var(xr)
#LLL
#mean(xr[1,sub[[4]]])
#mean(xr[1,sub[[5]]])
#c
#hist(rnorm(1:10000,7.420398,sqrt(1/3954.1225)),100,F)
#hist(rnorm(1:10000,7.254977,sqrt(1/3954.1225)),100,F)
#plot( density(rnorm(1:10000,7.420398,sqrt(1/3954.1225))) , typ="l",xlim=c(7,7.5))
#lines( density(rnorm(1:10000,7.254977,sqrt(1/3954.1225))) , typ="l",col=2)
#dnorm(7.419872,7.420398,sqrt(1/3954.1225))
##b
#plot( density(rnorm(1:10000,7.420573,sqrt(1/209.0072))) , typ="l",xlim=c(7,7.8))
#lines( density(rnorm(1:10000, 7.309942,sqrt(1/209.0072))) , typ="l",col=2)
##d
#plot( density(rnorm(1:10000,7.390271,sqrt(1/213.5900))) , typ="l",xlim=c(7,7.8))
#lines( density(rnorm(1:10000, 7.240192,sqrt(1/213.5900))) , typ="l",col=2)


restructur<-function(A){
  A<-t(A)
  colnames(A)<-c("a","b","c","d")
  return(A)
}
  
myIC<-function(logL,npar,k=2){
  if(k==2){
    #AIC== -2*loglikelihood + k*npar, npar represents the number of parameters in the fitted model, and k = 2 
    aic<- (-2)*logL + k*npar
    return(aic)
  }else if(k>2){
    #BIC == AIC(object, ..., k = log(nobs(object)))
    
    k=log(k) # k -> number of ovservations (hier 6)
    bic<- (-2)*logL + k*npar
    return(bic)
  }else{
    print(c(k,"<2 -- ERROR"))
  }
}
  
cal.myIC<-function(L,npar,k){
  IC<-c()
  IC[1]<-myIC(L[1],npar[1],k)
  IC[2]<-myIC(L[2],npar[2],k)
  IC[3]<-myIC(L[3],npar[2],k)
  IC[4]<-myIC(L[4],npar[2],k)
  return(IC)
}
 
make.IC<-function(L,npar,k){
    IC<-apply(L,1,cal.myIC,npar,k)
    IC<-restructur(IC)
    Ind<- apply(IC,MARGIN=1,FUN=minIndex)
    Names<-as.character(rownames(IC))
    print(table(Ind))
    return(list(IC=IC,minV=Ind,N=Names))
} 

fout<-function(A){
	if( is.matrix(A$SIG)){
		B<-cbind(A$N,A$SIG)
		B<-as.matrix(B)
		for(j in 2:5){
		  B[,j]<-signif(as.double(B[,j]),4) ## abschneiden der letzten Positionen 
		}
		colnames(B)<-c("ID","a","b","c","d")
	}
	else if( is.vector(A$SIG)){
		B<-c(A$N,A$SIG)
		names(B)<-c("ID","a","b","c","d")
	}
	return(B)
}
    
fout_2<-function(A){
  #if fout not work :(
  B.tmp<-matrix(signif(as.numeric(unlist(A$SIG)),4),length(A$N),4)
  B<-as.matrix(cbind(A$N, B.tmp))
  colnames(B)<-c("ID","a","b","c","d")
  return(B)
}
fout_3<-function(A){
  #if fout not work :(
  B.tmp<-matrix(signif(A,3),length(rownames(A)),4)
  B<-as.matrix(cbind(rownames(A), B.tmp))
  colnames(B)<-c("ID","a","b","c","d")
  return(B)
}
fout_4<-function(A,ID){
  #if fout not work :(
  B.tmp<-matrix(signif(A,3),length(rownames(A)),7)
  B<-as.matrix(cbind(rownames(A), as.character(IDs[rownames(A),2]), B.tmp))
  

  colnames(B)<-c("ID","Sym","a","b","c","d","FC","log2FC","geoMeanExp")
  return(B)
}

  #IN<-x.SBay.var
  #class<-3
	get.Sig.class<-function(IN,class,sig=0){
				#print(table(IN$minV))
			o<-IN$IC[IN$minV==class,]
				#print(range(o[,class]))
        #hist(o,200)
      print(length(o[,1]))
			
      if(sig==0){
					o.s<-o
			}else{
				 	#o.s<-o[o[,class]<log(sig),]
			    o.s<-o[o[,class]<sig,]
			}
			print(head(o.s))
			or<-order(o.s[,class],decreasing=F)
			o.sig<-o.s[or,]
			print(head(o.sig))
			Names<-as.character(rownames(o.sig))
				#print(length(Names))
			return(list(SIG=o.sig,N=Names))
	}
		

########### make test data 
#  get.probs<-function(a,range){
#
#    for(i in 1:a){
#    samp <- sample(1:a,range)
#    samp <- samp/sum(samp)
#    q<-rbind(q,samp)
#    }
#    return(q)
#  }

get.k<-function(nubV){
  q<-c()
  for(i in 1:nubV){
    q[i]<-sample(1:4,1)
  }
  return(q)
}
get.v<-function(tau,MUs){
  return(rnorm(1,MUs,1/sqrt(tau)))
}

get.Data<-function(k,tau,gam0,mu0){
  get.e<-function(k,tau,gam0,mu0){
  
  		#MUs<-rnorm(1:7,mu0,1/sqrt(tau*gam0))
  		MUs<-c() # for different mu0 
  		for(i in 1:length(mu0)){
  			MUs[i]<-rnorm(1,mu0[i],1/sqrt(tau*gam0))
  		}
  
      if(k==1){
        v<-rnorm(1:6,MUs[1],1/sqrt(tau))
      }else if(k==2){
        a<-c(2,3,2,3,2,3)  
        v<-c(get.v(tau,MUs[a[1]]),get.v(tau,MUs[a[2]]),get.v(tau,MUs[a[3]]),get.v(tau,MUs[a[4]]),get.v(tau,MUs[a[5]]),get.v(tau,MUs[a[6]])) 
      }else if(k==3){
        a<-c(4,5,4,4,4,5)  
        v<-c(get.v(tau,MUs[a[1]]),get.v(tau,MUs[a[2]]),get.v(tau,MUs[a[3]]),get.v(tau,MUs[a[4]]),get.v(tau,MUs[a[5]]),get.v(tau,MUs[a[6]])) 
      }else if(k==4){
        a<-c(6,7,6,6,6,6)  
        v<-c(get.v(tau,MUs[a[1]]),get.v(tau,MUs[a[2]]),get.v(tau,MUs[a[3]]),get.v(tau,MUs[a[4]]),get.v(tau,MUs[a[5]]),get.v(tau,MUs[a[6]])) 
      }
      return(v)
  }
  e<-apply(as.matrix(k),1,get.e,tau,gam0,mu0)
  d<-t(e)
  rownames(d)<-k
  
  return(d)
}  

get.MUS<-function(k,mu0,tau,gam0){
  MUs<-c() # for different mu0 
  for(i in 1:length(mu0)){
    MUs[i]<-rnorm(1,mu0[i],1/sqrt(tau*gam0))
  }
  return(MUs)
}
#MUS<-t(apply(as.matrix(k),1,get.MUS,mu0,tau,gam0))
#plot(density(MUS[,1]),typ="l")
#for(l in 2:length(mu0)){
#lines(density(MUS[,l]),col=l)
#}

 get.C.maxV<-function(cc){
	if(max(cc)>0){
			return( (which(cc==max(cc))+1) )
	}else{return(1)}
}				
get.ConfMat<-function(actual,predict){
	CofMa<-matrix(0,4,4,dimnames=list(c("a","b","c","d"),c("a","b","c","d")))
	for(i in 1:length(actual)){
		CofMa[actual[i],predict[i]]<-CofMa[actual[i],predict[i]]+1
	}
	return(CofMa)
}

get.ev<-function(k,AB){
			a<-AB$IC
			a.maxV<- apply(a,MARGIN=1,FUN=minIndex)
			CM<-get.ConfMat(k,a.maxV)
			return(CM)
}

#save_data
get.Time<-function(){

  mydate<-strsplit(date()," ")
  mytime<-strsplit(as.character(mydate[[1]][4]),":")
  std<-as.integer(mytime[[1]][1])
  min<-as.integer(mytime[[1]][2])
	day<-as.integer(mydate[[1]][3])

return(list(d=day,s=std,m=min))
	
}

saveSet<-function(d,k,vgl,tau){ 
  mydate<-strsplit(date()," ")
  mytime<-strsplit(as.character(mydate[[1]][4]),":")
  std<-as.integer(mytime[[1]][1])
  min<-as.integer(mytime[[1]][2])
  
  write.table(d,paste(sep="","~/Promotion/Kappler_Wichmann_Medizin/Auswertung/TESTDATA/Set_tau_",tau,"_",as.integer(mydate[[1]][3]),"_",std,":",min,"_d.txt"),row.names=F,col.names=F)
  write.table(k,paste(sep="","~/Promotion/Kappler_Wichmann_Medizin/Auswertung/TESTDATA/Set_tau_",tau,"_",as.integer(mydate[[1]][3]),"_",std,":",min,"_k.txt"))
  write.table(vgl,paste(sep="","~/Promotion/Kappler_Wichmann_Medizin/Auswertung/TESTDATA/Set_tau_",tau,"_",as.integer(mydate[[1]][3]),"_",std,":",min,"_vgl.txt"))
  pdf(paste(sep="","~/Promotion/Kappler_Wichmann_Medizin/Auswertung/TESTDATA/Set_tau_",tau,"_",as.integer(mydate[[1]][3]),"_",std,":",min,"_d.pdf"))
  hist(d)
  dev.off()
  
  #load_dat
  if(0){
    k<-read.table("~/Promotion/Kappler_Wichmann_Medizin/Auswertung/TESTDATA/k_21_15:57.txt",T)
    d<-as.matrix(read.table("~/Promotion/Kappler_Wichmann_Medizin/Auswertung/TESTDATA/d_21_15:57.txt",F))
    k<-k[,1]
    rownames(d)<-k
  }
  
}


f.input2 = function (p,q,name=c("A","B"),plotVENN=TRUE){
  input  <-list(A=p,B=q)
  names(input)<-name
  #print(input)
  if(plotVENN){ venn(input) }
  
  i<-intersect(input[[1]],input[[2]])
  s1<-setdiff(input[[1]],input[[2]])
  s2<-setdiff(input[[2]],input[[1]])
  
  return(list(inter=i,diffAB=s1,diffBa=s2))
}

f.input4 = function (p,q,r,s,t,name=c("A","B","C","D"),plotVENN=TRUE){
  input  <-list(A=p,B=q,C=r,D=s)
  names(input)<-name
  
  #if(plotVENN){ venn(input,simplify=TRUE) }
  if(plotVENN){ venn(input,simplify=FALSE) }
  #print(input)
  return(intersect4( p,q,r,s))
  
}

f.input5 = function (p,q,r,s,t,name=c("A","B","C","D","E"),plotVENN=TRUE){
  input  <-list(A=p,B=q,C=r,D=s,E=t)
  names(input)<-name
  
  #if(plotVENN){ venn(input,simplify=TRUE) }
  if(plotVENN){ venn(input,simplify=FALSE) }
  #print(input)
  return(intersect5( p,q,r,s,t))
  
}

f.input5.pretty = function (p,q,r,s,e,name,VennName=""){
  input  <-list(A=p,B=q,C=r,D=s,E=e)
  names(input)<-name
  
  require("VennDiagram")
  venn.plot <- venn.diagram(
    x = input,
    filename = paste(sep="_",VennName,"Venn_5set_pretty.tiff"),
    col = "black",
    fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    alpha = 0.50,
    cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,
            1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
    cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.05
  );
  return(input)  
}


f.input3.pretty = function (p,q,r,name,VennName=""){
  input  <-list(A=p,B=q,C=r)
  names(input)<-name
  venn.plot <- venn.diagram(
    x = input,
    euler.d = TRUE,
    filename = "Euler_3set_scaled.tiff",
    cex = 2.5,
    cat.cex = 2.5,
    cat.pos = 0
  );
}


f.input4.pretty = function (p,q,r,s,name,VennName=""){
  input  <-list(A=p,B=q,C=r,D=s)
  names(input)<-name
  
  require("VennDiagram")
  venn.plot <- venn.diagram(
    x = input,
    filename = paste(sep="_",VennName,"Venn_4set_pretty.tiff"),
    col = "transparent",
    fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
    alpha = 0.50,
    label.col = c("orange", "white", "darkorchid4", "white", 
                  "white", "white", "white", "white", "darkblue", "white", 
                  "white", "white", "white", "darkgreen", "white"),
    cex = 1.5,
    fontfamily = "serif",
    fontface = "bold",
    cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
    cat.cex = 1.5,
    cat.pos = 0,
    cat.dist = 0.07,
    cat.fontfamily = "serif",
    rotation.degree = 270,
    margin = 0.2
  );
}

f.input3 = function (p,q,r,name=c("A","B","C"),plotVENN=TRUE){
  input  <-list(A=p,B=q,C=r)
  names(input)<-name
  #print(input)
  
  if(plotVENN){ venn(input) }
  return(intersect3(p,q,r))
}

intersect5<-function(A,B,C,D,E){
  O<-intersect(A,intersect(B,intersect(C,intersect(D,E))))
  return(O)
}

intersect4<-function(A,B,C,D){
  O<-intersect(A,intersect(B,intersect(C,D)))
  return(O)
}

intersect3<-function(A,B,C){
  O<-intersect(A,intersect(B,C))
  return(O)
}

