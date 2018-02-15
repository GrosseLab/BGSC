head(X.TEST)

mfc<-cbind(X.TEST[,2]-X.TEST[,1],X.TEST[,4]-X.TEST[,3],X.TEST[,6]-X.TEST[,5])

head(mfc)
########################### c #############################
t1<-0.5;t2<-0.5;t3<-0.5;
t1<-0.5;t2<-1;t3<-0.5;
t1<-1;t2<-1;t3<-1;

mfc.c<-mfc[abs(mfc[,1])>t1,]
mfc.c<-mfc.c[abs(mfc.c[,2])<t2,]
mfc.c<-mfc.c[abs(mfc.c[,3])>t3,]

venn(f.input2(rownames(mfc.c),rownames(P_m_x_c),c("FC","BIC")))
venn(f.input2(rownames(mfc.c),rownames(P_m_x_c[P_m_x_c[,3]>0.9,]),c("FC","BIC")))

mfc.c.inter1<-intersect(rownames(mfc.c),rownames(P_m_x_c))
mfc.c.gene1<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.c.inter1),c(1,2)]

mfc.c.inter2<-intersect(rownames(mfc.c),rownames(P_m_x_c[P_m_x_c[,3]>0.9,]))
mfc.c.gene2<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.c.inter2),c(1,2)]

mfc.c.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.c)),c(1,2)]
mfc.c.gene3.b<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_c)),c(1,2)]
mfc.c.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_c[P_m_x_c[,3]>0.5,])),c(1,2)]
length(mfc.c.gene3[,1]);length(mfc.c.gene3.b[,1]);length(mfc.c.gene3.c[,1]);

################# d ###################
t1<-0.5;t2<-0.5;t3<-0.5;
t1<-0.75;t2<-0.5;t3<-0.5;
t1<-1;t2<-1;t3<-1;

mfc.d<-mfc[abs(mfc[,1])>t1,]
mfc.d<-mfc.d[abs(mfc.d[,2])<t2,]
mfc.d<-mfc.d[abs(mfc.d[,3])<t3,]

venn(f.input2(rownames(mfc.d),rownames(P_m_x_d),c("FC","BIC")))
venn(f.input2(rownames(mfc.d),rownames(P_m_x_d[P_m_x_d[,4]>0.9,]),c("FC","BIC")))

mfc.d.inter1<-intersect(rownames(mfc.d),rownames(P_m_x_d))
mfc.d.gene1<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.d.inter1),c(1,2)]

mfc.d.inter2<-intersect(rownames(mfc.d),rownames(P_m_x_d[P_m_x_d[,4]>0.9,]))
mfc.d.gene2<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.d.inter2),c(1,2)]

mfc.d.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.d)),c(1,2)]
mfc.d.gene3.b<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_d)),c(1,2)]
mfc.d.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_d[P_m_x_d[,4]>0.5,])),c(1,2)]
length(mfc.d.gene3[,1]);length(mfc.d.gene3.b[,1]);length(mfc.d.gene3.c[,1]);


################# b ###################
t1<-0.5;t2<-0.5;t3<-0.5;
t1<-0.75;t2<-0.5;t3<-0.5;
t1<-1;t2<-1;t3<-1;

mfc.b<-mfc[abs(mfc[,1])>t1,]
mfc.b<-mfc.b[abs(mfc.b[,2])>t2,]
mfc.b<-mfc.b[abs(mfc.b[,3])>t3,]

venn(f.input2(rownames(mfc.b),rownames(P_m_x_b),c("FC","BIC")))
venn(f.input2(rownames(mfc.b),rownames(P_m_x_b[P_m_x_b[,2]>0.9,]),c("FC","BIC")))

mfc.b.inter1<-intersect(rownames(mfc.b),rownames(P_m_x_b))
mfc.b.gene1<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.d.inter1),c(1,2)]

mfc.b.inter2<-intersect(rownames(mfc.d),rownames(P_m_x_b[P_m_x_b[,2]>0.9,]))
mfc.b.gene2<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.d.inter2),c(1,2)]

mfc.b.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.b)),c(1,2)]
mfc.b.gene3.b<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_b)),c(1,2)]
mfc.b.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_b[P_m_x_b[,2]>0.5,])),c(1,2)]
length(mfc.b.gene3[,1]);length(mfc.b.gene3.b[,1]);length(mfc.b.gene3.c[,1]);


################### a 

t1<-0.5;t2<-0.5;t3<-0.5;
t1<-0.25;t2<-0.25;t3<-0.25;
t1<-1;t2<-1;t3<-1;

mfc.a<-mfc[abs(mfc[,1])<t1,]
mfc.a<-mfc.a[abs(mfc.a[,2])<t2,]
mfc.a<-mfc.a[abs(mfc.a[,3])<t3,]

venn(f.input2(rownames(mfc.a),rownames(P_m_x_a),c("FC","BIC")))
venn(f.input2(rownames(mfc.a),rownames(P_m_x_a[P_m_x_a[,1]>0.9,]),c("FC","BIC")))

mfc.a.inter1<-intersect(rownames(mfc.a),rownames(P_m_x_a))
mfc.a.gene1<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.d.inter1),c(1,2)]

mfc.a.inter2<-intersect(rownames(mfc.d),rownames(P_m_x_a[P_m_x_a[,1]>0.9,]))
mfc.a.gene2<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),mfc.d.inter2),c(1,2)]

mfc.a.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.a)),c(1,2)]
mfc.a.gene3.b<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_a)),c(1,2)]
mfc.a.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_a[P_m_x_a[,1]>0.5,])),c(1,2)]
length(mfc.a.gene3[,1]);length(mfc.a.gene3.b[,1]);length(mfc.a.gene3.c[,1]);


#
total<-length(intersect(as.character(EGFR_genes.IDs[,1]),rownames(X.TEST)))
mfc.RES<-c()
mfc.RES<-rbind(mfc.RES,c(length(mfc.a.gene3[,1]),length(mfc.a.gene3.b[,1]),length(mfc.a.gene3.c[,1])))
mfc.RES<-rbind(mfc.RES,c(length(mfc.b.gene3[,1]),length(mfc.b.gene3.b[,1]),length(mfc.b.gene3.c[,1])))
mfc.RES<-rbind(mfc.RES,c(length(mfc.c.gene3[,1]),length(mfc.c.gene3.b[,1]),length(mfc.c.gene3.c[,1])))
mfc.RES<-rbind(mfc.RES,c(length(mfc.d.gene3[,1]),length(mfc.d.gene3.b[,1]),length(mfc.d.gene3.c[,1])))

mfc.RES[1,]/total
colSums(mfc.RES[-1,])/total
colSums(mfc.RES)/total


### case
total<-length(intersect(as.character(EGFR_genes.IDs[,1]),rownames(X.TEST)))
#total<-length(unique(EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(X.TEST)),c(2)]))
#t1<-0.3;t2<-0.3;t3<-0.3;

mfc.gene3<-c()
mfc.gene3.ALL<-c()
tseq<-seq(0.1,1,by=.05)
L.mfc<-list()
for(i in tseq){
  
  t1<-i;t2<-i;t3<-i;
  
  mfc.a<-mfc[abs(mfc[,1])<t1,]
  mfc.a<-mfc.a[abs(mfc.a[,2])<t2,]
  mfc.a<-mfc.a[abs(mfc.a[,3])<t3,]
  mfc.b<-mfc[abs(mfc[,1])>t1,]
  mfc.b<-mfc.b[abs(mfc.b[,2])>t2,]
  mfc.b<-mfc.b[abs(mfc.b[,3])>t3,]
  mfc.c<-mfc[abs(mfc[,1])>t1,]
  mfc.c<-mfc.c[abs(mfc.c[,2])<t2,]
  mfc.c<-mfc.c[abs(mfc.c[,3])>t3,]
  mfc.d<-mfc[abs(mfc[,1])>t1,]
  mfc.d<-mfc.d[abs(mfc.d[,2])<t2,]
  mfc.d<-mfc.d[abs(mfc.d[,3])<t3,]
  
  mfc.a.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.a)),c(1,2)]
  mfc.b.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.b)),c(1,2)]
  mfc.c.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.c)),c(1,2)]
  mfc.d.gene3<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.d)),c(1,2)]
  
  tmp<-c(length(mfc.a.gene3[,1]),length(mfc.b.gene3[,1]),length(mfc.c.gene3[,1]),length(mfc.d.gene3[,1]))
  #tmp<-c(length(unique(mfc.a.gene3[,2])),length(unique(mfc.b.gene3[,2])),length(unique(mfc.c.gene3[,2])),length(unique(mfc.d.gene3[,2])))
  
  #venn(f.input4(rownames(mfc.a),rownames(mfc.b),rownames(mfc.c),rownames(mfc.d),c("a","b","c","d")))
  #venn(f.input4(intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.a)),intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.b)),intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.c)),intersect(as.character(EGFR_genes.IDs[,1]),rownames(mfc.d)),c("a","b","c","d")))
  
  mfc.gene3<-rbind(mfc.gene3,tmp)
  
  tmp<-c(length(mfc.a[,1]),length(mfc.b[,1]),length(mfc.c[,1]),length(mfc.d[,1]))
  mfc.gene3.ALL<-rbind(mfc.gene3.ALL,tmp)
  
  L.mfc[[as.character(i)]][[1]]<-mfc.a.gene3
  L.mfc[[as.character(i)]][[2]]<-mfc.b.gene3
  L.mfc[[as.character(i)]][[3]]<-mfc.c.gene3
  L.mfc[[as.character(i)]][[4]]<-mfc.d.gene3
  names(L.mfc[[as.character(i)]])<-c("a","b","c","d")
}
colnames(mfc.gene3)<-c("a","b","c","d")
rownames(mfc.gene3)<-tseq
colnames(mfc.gene3.ALL)<-c("a","b","c","d")
rownames(mfc.gene3.ALL)<-tseq
mfc.gene3[,3]/mfc.gene3.ALL[,3]


mfc.gene3.c<-c()
mfc.gene3.c.ALL<-c()
L.mfc.c<-list()
tseq2<-seq(0.1,1,by=.05)
for(t in tseq2){
  mfc.a.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_a[P_m_x_a[,1]>t,])),c(1,2)]
  mfc.b.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_b[P_m_x_b[,2]>t,])),c(1,2)]
  mfc.c.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_c[P_m_x_c[,3]>t,])),c(1,2)]
  mfc.d.gene3.c<-EGFR_genes.IDs[intersect(as.character(EGFR_genes.IDs[,1]),rownames(P_m_x_d[P_m_x_d[,4]>t,])),c(1,2)]
  
  tmp<-c(length(mfc.a.gene3.c[,1]),length(mfc.b.gene3.c[,1]),length(mfc.c.gene3.c[,1]),length(mfc.d.gene3.c[,1]))
  #tmp<-c(length(unique(mfc.a.gene3.c[,2])),length(unique(mfc.b.gene3.c[,2])),length(unique(mfc.c.gene3.c[,2])),length(unique(mfc.d.gene3.c[,2])))
  mfc.gene3.c<-rbind(mfc.gene3.c,tmp)

  tmp<-c(length(rownames(P_m_x_a[P_m_x_a[,1]>t,])),length(rownames(P_m_x_b[P_m_x_b[,2]>t,])),length(rownames(P_m_x_c[P_m_x_c[,3]>t,])),length(rownames(P_m_x_d[P_m_x_d[,4]>t,])))
  mfc.gene3.c.ALL<-rbind(mfc.gene3.c.ALL,tmp)
  
  
  L.mfc.c[[as.character(i)]][[1]]<-mfc.a.gene3.c
  L.mfc.c[[as.character(i)]][[2]]<-mfc.b.gene3.c
  L.mfc.c[[as.character(i)]][[3]]<-mfc.c.gene3.c
  L.mfc.c[[as.character(i)]][[4]]<-mfc.d.gene3.c
  names(L.mfc.c[[as.character(i)]])<-c("a","b","c","d")
  
}  
colnames(mfc.gene3.c)<-c("a","b","c","d")
rownames(mfc.gene3.c)<-tseq2
colnames(mfc.gene3.c.ALL)<-c("a","b","c","d")
rownames(mfc.gene3.c.ALL)<-tseq2
mfc.gene3.c[,3]/mfc.gene3.c.ALL[,3]


rowSums(mfc.gene3.c[,-1])
rowSums(mfc.gene3[,-1])
rowSums(mfc.gene3.c[,c(3,4)])
rowSums(mfc.gene3[,c(3,4)])
mfc.gene3.c[,c(3)]
mfc.gene3[,c(3)]

tseq3<-seq(0.4,1,by=.05)
plot(tseq3,mfc.gene3.c[as.character(tseq3),3],type="l",col=2,ylim=c(0,40))
lines(tseq3,mfc.gene3[as.character(tseq3),3])
plot(tseq3,mfc.gene3.c[as.character(tseq3),4],type="l",col=2,ylim=c(0,40))
lines(tseq3,mfc.gene3[as.character(tseq3),4])
plot(tseq3,rowSums(mfc.gene3.c[as.character(tseq3),c(3,4)]),type="l",col=2,ylim=c(0,45))
lines(tseq3,rowSums(mfc.gene3[as.character(tseq3),c(3,4)]))
plot(tseq3,rowSums(mfc.gene3.c[as.character(tseq3),c(1,2)]),type="l",col=2,ylim=c(0,250))
lines(tseq3,rowSums(mfc.gene3[as.character(tseq3),c(1,2)]))

Sen.mfc.gene3.c<-rowSums(mfc.gene3.c[,c(2,3,4)])/(mfc.gene3.c[,1]+rowSums(mfc.gene3.c[,c(2,3,4)]))
Sen.mfc.gene3<-rowSums(mfc.gene3[,c(2,3,4)])/(mfc.gene3[,1]+rowSums(mfc.gene3[,c(2,3,4)]))

pdf(paste(sep="",out,"Compare.pdf"),14,8)
barplot(t(mfc.gene3.c),beside=T,legend.text=T,col=brewer.pal(9,"Set1")[1:4])
barplot(t(mfc.gene3),beside=T,legend.text=T,col=brewer.pal(9,"Set1")[1:4])
plot(tseq2,Sen.mfc.gene3.c,type="l",col=2,lty=2,ylim=c(0,1),xlab="t is: 1. FoldChange threshold \nt is: 2. prob for class threshold",ylab="SEN")
lines(tseq2,Sen.mfc.gene3,type="l",col=1)
legend("top",c("Prob Method","FoldChange Method"),lty=c(2,1),col=c(2,1))
dev.off()


MM<-mfc.gene3.c[,3]/mfc.gene3.c.ALL[,3]
NN<-mfc.gene3[,3]/mfc.gene3.ALL[,3]

MM<-rowSums(mfc.gene3.c[,c(3,4)])/rowSums(mfc.gene3.c.ALL[,c(3,4)])
NN<-rowSums(mfc.gene3[,c(3,4)])/rowSums(mfc.gene3.ALL[,c(3,4)])

plot(tseq2,MM,type="l",col=2,lty=2,ylim=c(0,.1),xlab="t is: 1. FoldChange threshold \nt is: 2. prob for class threshold",ylab="SEN")
lines(tseq2,NN,type="l",col=1)
legend("top",c("Prob Method","FoldChange Method"),lty=c(2,1),col=c(2,1))
