main <- function(){
  
  ### normalize Data ----------------------------------------------------------------------
    normData <- normalizeExpData()
  
  ### calulate logLik for class  Data ----------------------------------------------------------------------
    Lsets <- get.Lset()
    normDataLogLik <- logLikelihoodOfnormData(normData$E)
  
    rownames(ALL.MUs) <- rownames(normData$E)
    colnames(ALL.MUs) <- c('a0','a1','b0','b1','c0','c1','d0','d1')
    rownames(ALL.VARs) <- rownames(normData$E)
    colnames(ALL.VARs) <- c('a','b','c','d')
    
  ### calulate BIC from logLik  ----------------------------------------------------------------------
    npar <- sapply(Lsets, function(x) sum(!sapply(x,is.null ) )) + 1  ## number parapeters for LogLilk -> mean + var 
    k <-  sapply(Lsets, function(x) sum( sapply(x,length) ))
    normDataBIC <- get.IC(normDataLogLik , npar, k , IC = 'BIC')
    BICminInd <- apply( normDataBIC,MARGIN = 1, FUN = minIndex)
    print(table(BICminInd))
  
  ### calulate Posterior from BIC ----------------------------------------------------------------------
    normDataPosterior <- get.Posterior( normDataBIC ,Pis = c(0.7,0.1,0.1,0.1))
    POSTmaxInd <- apply(normDataPosterior, MARGIN = 1 ,FUN = maxIndex)
    print(table(POSTmaxInd))
    
    PostClass <- get.gene.classes(data = normDataPosterior,indexing = "max",filter = 0.75, DoPlot = TRUE)
  
  ### compare to qPCR ----------------------------------------------------------------------
    MeanFoldChangeClass <- get.log2Mean.and.log2FC(normData = normData)
    qCPRdataC <- getQPCR()
    
    IDs.dt <- data.table::data.table(normData$genes,keep.rownames = T,key = 'rn')
    IDs.dt.c <- IDs.dt[rownames(PostClass$resFilter$c),]
    data.table::setkey(IDs.dt.c,'SYMBOL')
    qgenesIDs <- lapply(qCPRdataC$rn, function(qg) as.character(IDs.dt.c[qg,][['rn']]) )
    names(qgenesIDs) <- qCPRdataC$rn
  
    make.plot.data.FC.Ill.qPCR <- function(qCPRdata,MeanFoldChangeClass, class="c"){
      TMP <- data.frame()
      for (tmpG in qCPRdataC$rn  ) {
        expC <- MeanFoldChangeClass[[class]][qgenesIDs[[tmpG]], ]
        TMPexp <- data.frame( "Gene" = tmpG,
                              "ExpID" = expC$rn,
                              "Set"="Illumina",
                              "FC" =expC[["s1s0FC"]] ,
                              "stderr" = expC[["s0s1stderr"]]
        )
        
        TMPqpc <- data.frame( "Gene" = tmpG,
                              "Set" = "RT-qPCR",
                              "FC" = qCPRdataC[tmpG,][["relative.Werte.geoMean.C1C0.logFC"]] ,
                              #"stderr" = qCPRdataC[tmpG,][["relative.Werte.pooeld.var.C0C1_SatterthwaiteApproximation"]],
                              "stderr" = qCPRdataC[tmpG,][["s0s1stderr"]]
        )
        
        for (i in 1:nrow(TMPexp)) {
          tmp <- TMPqpc 
          tmp$ExpID <-  as.character(TMPexp[i,]$ExpID)
          TMP <- rbind(TMP , rbind(TMPexp[i, ],tmp[,  colnames(TMPexp) ]))
        }
      }  
      # TMP <- TMP[ TMP$Gene != 'KIF5C',]
      TMP$pid <- paste0(TMP$Gene,'::',TMP$ExpID)
      TMP$pid <- factor(TMP$pid)
      return(TMP)
      
    }
    
    make.plot.data.exp.Ill.qPCR <- function(plotGens,MeanFoldChangeClass, class="c"){
      TMP <- data.frame()
      for (tmpG in qCPRdataC$rn  ) {
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
    
    
  ### ggplot2 theme  ----------------------------------------------------------------------
    require(ggplot2)
    cols <- c(RColorBrewer::brewer.pal(9,"RdGy")[8],RColorBrewer::brewer.pal(9,"Reds")[6])
    base_size <- 12
    .thememap <- function(base_size = 12, legend_key_size = 0.4, base_family = "", col = "grey70") {
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
                       panel.grid.major =  ggplot2::element_line(colour = col, size = 0.2),
                       panel.grid.minor =  ggplot2::element_blank(),
                       panel.background = ggplot2::element_rect(fill="white",size = 0.2,),
                       # #panel.grid.minor.y = element_line(size=3),
                       # panel.grid.major = ggplot2::element_line(colour = "white"),
                       
                       # Force the plot into a square aspect ratio
                       # aspect.ratio = 1,
                       aspect.ratio = 9 / 16,
                       
                       legend.key.size  = ggplot2::unit(legend_key_size, "cm"),
                       plot.title = ggplot2::element_text(hjust = 0.5 , vjust= 1)          
        )
    }
  
  ### bar log2 FC qPCR Illumina  ----------------------------------------------------------------------
    PlotDataFC <- make.plot.data.FC.Ill.qPCR(qCPRdata,MeanFoldChangeClass)
    
    rns <- levels(PlotDataFC$Gene)[c(2,5,6,1,3,4)] ; PlotDataFC$Gene <- factor(PlotDataFC$Gene, levels = rns)
    rns <- levels(PlotDataFC$pid)[c(2,5,6,1,3,4)]  ; PlotDataFC$pid <- factor(PlotDataFC$pid, levels = rns)
    
    # cols<-c(RColorBrewer::brewer.pal(9,"RdGy")[8],RColorBrewer::brewer.pal(9,"Reds")[6])
    cols <- RColorBrewer::brewer.pal(11,"PRGn")[c(2,10)]
    limits <- aes(ymax = PlotDataFC$FC + PlotDataFC$stderr , ymin=PlotDataFC$FC - PlotDataFC$stderr)
    dodge <- position_dodge(width=0.9)
    
    # g1 <- ggplot(PlotDataFC,aes(x=factor(pid),y=FC,fill=factor(Set))) +
    g1 <- ggplot(PlotDataFC,aes(x=factor(Gene),y=FC,fill=factor(Set))) +  
      geom_bar(position="dodge",stat="identity") + 
      geom_errorbar(limits, position=dodge, width=0.25) +
      ylim(c(-2 ,2)) +
      scale_fill_manual(values = cols,name="") + 
      labs(x = "", y = "Log2-fold change") +  
      .thememap(base_size = base_size,legend_key_size = 0.6) + 
      theme(legend.position="bottom") 
    g1    
  
  ### bar log2 FC qPCR Illumina  ----------------------------------------------------------------------
    PlotDataEXP <- make.plot.data.exp.Ill.qPCR(qCPRdata,MeanFoldChangeClass)
    
    rns <- levels(PlotDataEXP$Gene)[c(2,5,6,1,3,4)] ; PlotDataEXP$Gene <- factor(PlotDataEXP$Gene, levels = rns)
    rns <- levels(PlotDataEXP$pid)[c(2,5,6,1,3,4)]  ; PlotDataEXP$pid <- factor(PlotDataEXP$pid, levels = rns)
    
    # cols <- RColorBrewer::brewer.pal(11,"PRGn")[c(2,10)]
    cols <-  RColorBrewer::brewer.pal(8,'Paired')[c(6,2)] 
    
    limits <- aes(ymax = PlotDataEXP$Mean + PlotDataEXP$stderr , ymin=PlotDataEXP$Mean - PlotDataEXP$stderr)
    dodge <- position_dodge(width=0.9)
    
    g2 <- ggplot(PlotDataEXP,aes(x=factor(Gene),y=Mean,fill=factor(Set))) +  
      geom_bar(position="dodge",stat="identity") + 
      geom_errorbar(limits, position=dodge, width=0.25) +
      ylim(c(0 ,10)) +
      scale_fill_manual(values = cols,name="") + 
      labs(x = "", y = "Log2 mean expression") +  
      .thememap(base_size = base_size,legend_key_size = 0.6) + 
      theme(legend.position="bottom") 
    g2    

    gridExtra::grid.arrange(g2  , g1, ncol=1)#,top=tair)
    
  
  ### pheatmap of normData  ----------------------------------------------------------------------
      Leset <- get.Lset()
      pmat <- t(matrix(c(0,1,0,0,0,1),2,3))
      dimnames(pmat) <- list( c('no RNAi','RNAi by siRNA ALL','RNAi by siRNA I'),c('no EGF','EGF'))  
      pmat <- pmat[c(1,3,2),] ## to get same oder as in paper
      bk = seq(0,max(1),by = 1)
      hmcols <- c(RColorBrewer::brewer.pal(9,"Reds")[6],RColorBrewer::brewer.pal(9,"Blues")[7])# colorRampPalette(c("blue", "red"))(length(bk))
      pheatmap::pheatmap(pmat,color = hmcols,breaks = seq(-1,1,by = 1) ,display_numbers = T,fontsize_number = 20,number_format = '%.0f',cluster_rows = FALSE,cluster_cols = FALSE,main='expression pattern group c',fontsize = 12,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black',legend = FALSE)
      
      # grid.C <- pheatmap::pheatmap(pmat,color = hmcols,breaks = seq(-1,1,by = 1) ,display_numbers = T,fontsize_number = 20,number_format = '%.0f',cluster_rows = FALSE,cluster_cols = FALSE,main='expression pattern group c',fontsize = 12,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black',legend = FALSE ,silent = T)
      
      library(gridExtra)
      library(grid)
      library(gtable)
    
      Gene.gtable <- list()
      for( tmoG in names(qgenesIDs)){
        pmat <- t(matrix(normData$E[qgenesIDs[[tmoG]],],2,3))
        dimnames(pmat) <- list( c('no RNAi','RNAi by siRNA ALL','RNAi by siRNA I'),c('no EGF','EGF'))  
        pmat <- pmat[c(1,3,2),] ## to get same oder as in paper
        
        pmat <- 2^pmat
        tmpnScore <-  ceiling(pmat/100)*100
        # tmpnScore <-  ceiling(pmat/10)*10
        
        bk = seq(0,max(tmpnScore),by = 1)
        hmcols <- colorRampPalette(c("white", "darkred"))(length(bk)-1)
        pheatmap::pheatmap(pmat,color = hmcols,breaks = seq(0,max(tmpnScore),by = 1) ,display_numbers = T,fontsize_number = 12,cluster_rows = FALSE,cluster_cols = FALSE,main=tmoG,fontsize = 12,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black', )
                           # filename = paste0(plotDir,'/Heatmap_Score_',tair,'.pdf') ,width = 10 ,height = 10)
        Gene.gtable[[tmoG]] <- pheatmap::pheatmap(pmat,silent = T,color = hmcols,breaks = seq(0,max(tmpnScore),by = 1) ,display_numbers = T,fontsize_number = 15,cluster_rows = FALSE,cluster_cols = FALSE,main=tmoG,fontsize = 12,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black')
      }
  
      # grid.arrange(
      #   Gene.gtable$CKAP2L$gtable,
      #   Gene.gtable$ROCK1$gtable,
      #   Gene.gtable$TPR$gtable,
      #   Gene.gtable$ALDH4A1$gtable,
      #   Gene.gtable$CLCA2$gtable,
      #   Gene.gtable$GALNS$gtable,
      #   nrow=2,ncol=3)
      
      addBorder_gtable <- function(g){ g <- gtable_add_grob(g,
                                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                          t= 1 , b = nrow(g), l = 1, r = ncol(g))
                    # g <- gtable_add_grob(g,
                    #                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                    #                      t = 1, l = 1, r = ncol(g))
                    
                    return(g)
      }
      grid.arrange(
        addBorder_gtable(Gene.gtable$CKAP2L$gtable),
        addBorder_gtable(Gene.gtable$ROCK1$gtable),
        addBorder_gtable(Gene.gtable$TPR$gtable),
        addBorder_gtable(Gene.gtable$ALDH4A1$gtable),
        addBorder_gtable(Gene.gtable$CLCA2$gtable),
        addBorder_gtable(Gene.gtable$GALNS$gtable),
        nrow=2,ncol=3)              

  ### DES   ----------------------------------------------------------------------
      
      IDs <- do.call(c,qgenesIDs)
      
      i <- 6
      id <- "ILMN_1730999"
      
      
      tmpDensity <- lapply(qgenesIDs, function(id) Density.NV.fit.plot(id = id,normData,DOplot = TRUE) )
      
      aaa <- Density.NV.fit.plot(id = 'ILMN_1687840',normData,useGroup = "a",DOplot = TRUE)
      bbb <- Density.NV.fit.plot(id = 'ILMN_1684585',normData,useGroup = "b",DOplot = TRUE)
      ccc <- Density.NV.fit.plot(id = 'ILMN_1730999',normData,useGroup = "c",DOplot = TRUE)
      ddd <- Density.NV.fit.plot(id = 'ILMN_2320964',normData,useGroup = "d",DOplot = TRUE)
      
      grid.arrange(aaa + theme(legend.position="none")  , bbb + theme(legend.position="none") , ccc + theme(legend.position="none") , ddd + theme(legend.position="none") ,ncol=2,nrow=2)
      
      head(PostClass$resFilter$d,20)
      
  ### TEST   ----------------------------------------------------------------------
      

  
  ### qPCR
  RawQPCR
  
  
  
  
  ###   enrichment analysis   ----------------------------------------------------------------------

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
  
  
  
  # library("org.Hs.eg.db")
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
  f.input4(EGFR_overlapping_GO[['a']]$Term,
           EGFR_overlapping_GO[['b']]$Term,
           EGFR_overlapping_GO[['c']]$Term,
           EGFR_overlapping_GO[['d']]$Term,vennOut = T)
  
  tmp <- REFSEQ_MRNA_to_Illm(EGFR_overlapping_GO[['c']][EGFR_overlapping_GO[['c']]$Term == "GO:0005524~ATP binding",]$Genes)
  
  normData$E[intersect(rownames(normData$E),tmp),]
  
  table(POSTmaxInd[intersect(rownames(normData$E),tmp)])
  
  sapply(ChartReportCategoryList[[ChartReportCategory]][[i]]$Genes, function(x) REFSEQ_MRNA_to_Illm(x,asString = TRUE))
  
}