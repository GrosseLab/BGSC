main <- function(){
  
  require(ggplot2)
  library(gridExtra)
  library(grid)
  library(gtable)
  # library(BGSC)
  
  runPredection <- function(set){
    ### normalize Data ----------------------------------------------------------------------
    # normData <- normalizeExpData()
    # normData <- normalizeExpData(DetectionPvalNumber = 1)
    normData <- normalizeExpData(set = set)
    
    ### calulate logLik for class  Data ----------------------------------------------------------------------
    Lsets <- get.Lset()
    IndicatorVar <- lapply(Lsets,function(x){
      vec <- rep(0,6);names(vec) <- c(1:6)
      if (!is.null(x$s1)) { vec[x$s1] <- 1 }
      return(vec)
    })
    # pheatmap::pheatmap(t(do.call(rbind, IndicatorVar)),cluster_rows = F,cluster_cols = F)
    
    normDataLogLikData <- logLikelihoodOfnormData(normData$E)
    normDataLogLik <- normDataLogLikData[['logL']]
    ALL.MUs  <- normDataLogLikData[['ALL.MUs']]
    ALL.VARs <- normDataLogLikData[['ALL.VARs']]
    
    # TPRid = 'ILMN_1730999'
    # exp(normDataLogLik[TPRid,]) 
    
    ### calulate BIC from logLik  ----------------------------------------------------------------------
    npar <- sapply(Lsets, function(x) sum(!sapply(x,is.null ) )) + 1  ## number parapeters for LogLilk -> mean + var 
    k <-  sapply(Lsets, function(x) sum( sapply(x,length) ))
    normDataBIC <- get.IC(normDataLogLik , npar, k , IC = 'BIC')
    BICminInd <- apply( normDataBIC,MARGIN = 1, FUN = minIndex)
    print(table(BICminInd))
    
    # normDataBIC[TPRid,]
    # exp(normDataLogLik[TPRid,]) / ( sqrt(6)^npar )
    
    ### calulate Posterior from BIC ----------------------------------------------------------------------
    # normDataPosterior <- get.Posterior( normDataBIC ,Pis = c(0.7,0.1,0.1,0.1))
    pis = c(0.7,0.1,0.1,0.1)
    normDataPosterior <- get.Posterior( normDataBIC ,Pis = pis)
    
    POSTmaxInd <- apply(normDataPosterior, MARGIN = 1 ,FUN = maxIndex)
    print(table(POSTmaxInd))
    
    PostClass <- get.gene.group(data = normDataPosterior,indexing = "maximal",filter = 0.75, DoPlot = TRUE)
    
    return(list("PostClass" = PostClass,'normData'=normData,'ALL.MUs' = ALL.MUs,'ALL.VARs' = ALL.VARs,'normDataPosterior' = normDataPosterior,'normDataBIC' = normDataBIC))
  }
  
  set <- 'SF767'
  predSF767 <- runPredection(set)
  normData <- predSF767$normData
  ALL.MUs <- predSF767$ALL.MUs
  ALL.VARs <- predSF767$ALL.VARs
  normDataBIC <- predSF767$normDataBIC
  normDataPosterior <- predSF767$normDataPosterior
  PostClass <- predSF767$PostClass
  
  writePostClass <- function(PredectionRes,set){
    MeanFoldChangeClass <- get.log2Mean.and.log2FC(normData = PredectionRes$normData)

    outRes <- list()
    for(class in names(PredectionRes$PostClass$res)){
      tmp <- round(PredectionRes$PostClass$res[[class]],4)
      colnames(tmp) <- paste0("Posterior_class_",colnames(tmp))
      tmpkey <- paste0("Posterior_class_",class)
      OutOrderID <- sort(tmp[,tmpkey],decreasing = T)
      tmp <- as.data.table(tmp,keep.rownames = T)
      setkey(tmp,'rn')
    
      # tmpMu <- PredectionRes$ALL.MUs
      # tmpMu <- tmpMu[,colnames(tmpMu)[stringr::str_detect(string = colnames(tmpMu),pattern = class)  ]]
      # colnames(tmpMu) <- paste0("MeanExp_",colnames(tmpMu))
      # removeCol <- colnames(tmpMu)[is.nan(tmpMu[1,])]
      # tmpMu <- as.data.table(tmpMu,keep.rownames = T)
      # if(length(removeCol)>0){
      #   for(i in removeCol){
      #     tmpMu <- tmpMu[, !i, with=FALSE] 
      #   }
      # }
      # setkey(tmpMu,'rn')
      # tmp <-  merge(tmpMu,tmp,by.x = 'rn',by.y = 'rn')
      # setkey(tmp,'rn')
      
      tmpFC <- MeanFoldChangeClass[[class]]
      colnames(tmpFC) <- stringr::str_replace_all(colnames(tmpFC),pattern = 's',replacement = class)
      colnames(tmpFC) <- stringr::str_replace_all(colnames(tmpFC),pattern = 'M',replacement = "Mean")
      colnames(tmpFC) <- stringr::str_replace_all(colnames(tmpFC),pattern = 'FC',replacement = "log2FC")
      colnames(tmpFC) <- stringr::str_replace_all(colnames(tmpFC),pattern = paste0(class,'tderr'),replacement = 'stderr')
      
      setkey(tmpFC,'rn')
      tmp <-  merge(tmpFC,tmp,by.x = 'rn',by.y = 'rn')
      
      tmpSYMBOL <- as.data.table(PredectionRes$normData$genes[tmp$rn,],keep.rownames = T)
      setkey(tmpSYMBOL,'rn')
      tmp <-  merge(tmpSYMBOL[,.(rn,SYMBOL)],tmp,by.x = 'rn',by.y = 'rn')
      setkey(tmp,'rn')
      
      tmp <- tmp[names(OutOrderID),]
      outRes[[class]] <- tmp
      
      colnames(tmp)[1] <- 'ID'
      # write.csv2(tmp,paste0('./PredectionResult/PredectionResultClass_',class,'.csv') )    
      write.csv2( format(tmp, digits=3),paste0('./PredectionResult/PredectionResultClass_',class,'.csv') )    
    }
    
    return(outRes)
  }
  predSF767outRes <- writePostClass(predSF767,set)
  
  ### compare to qPCR ----------------------------------------------------------------------
    MeanFoldChangeClass <- get.log2Mean.and.log2FC(normData = normData)
    
    data(qPCR_SF767, envir = environment()) 
    qgenes <-  c('CKAP2L','ROCK1','TPR','ALDH4A1','CLCA2','GALNS')
    comparative_qPCR <- list()
    qPCR_Mean <- qPCR_log2FC <- data.frame()
    for(n in qgenes){
      Ref_GAPDH <- as.double(qPCR_SF767[,'GAPDH'])
      
      tmp <-  as.double(qPCR_SF767[,n])
      comparative_qPCR[[n]] <- comparativeMethod_qPCR.analysis(Gene = tmp,Ref = Ref_GAPDH)
      comparative_qPCR[[n]]$CellLine <- qPCR_SF767$CellLine
      comparative_qPCR[[n]]$Treatment <- qPCR_SF767$Treatment
    } 
    SF.log2FC <- comparativeMethod_qPCR.RNAi.log2FC(comparative_qPCR)
    qCPRdataC <- SF.log2FC$ddCT.C1C0
    data.table::setkey(qCPRdataC,Gene)
    # qCPRdataC <- getQPCR()

    IDs.dt <- data.table::data.table(normData$genes,keep.rownames = T,key = 'rn')
    IDs.dt.c <- IDs.dt[rownames(PostClass$resFilter$c),]
    data.table::setkey(IDs.dt.c,'SYMBOL')
    qgenesIDs <- lapply(as.character(qCPRdata$Gene), function(qg) as.character(IDs.dt.c[qg,][['rn']]) )
    names(qgenesIDs) <- as.character(qCPRdata$Gene)
  
    hist(MeanFoldChangeClass$c$s1s0FC,50);abline(v=c(0.5,-.5))
    
    MeanFoldChangeClass$c[abs(MeanFoldChangeClass$c$s1s0FC) > 0.5,]
    
  ### bar log2 FC qPCR Illumina  ----------------------------------------------------------------------
    PlotDataFC <- make.plot.data.FC.Ill.qPCR(qCPRdata = qCPRdataC, qgenesIDs = qgenesIDs, MeanFoldChangeClass, class="c" )
    
    cor.test( PlotDataFC[PlotDataFC$Set=='Microarray','FC'],PlotDataFC[PlotDataFC$Set=='qPCR','FC'] ) 
    
    # print(PlotDataFC)
    
    rns <- qgenes #levels(PlotDataFC$Gene)[c(2,5,6,1,3,4)] ; 
    PlotDataFC$Gene <- factor(PlotDataFC$Gene, levels = rns)
    # rns <- levels(PlotDataFC$pid)[c(2,5,6,1,3,4)]  ; PlotDataFC$pid <- factor(PlotDataFC$pid, levels = rns)
    
    # cols <- RColorBrewer::brewer.pal(11,"PRGn")[c(2,10)]
    cols <- RColorBrewer::brewer.pal(11,"PRGn")[c(3,9)]
    limits <- aes(ymax = PlotDataFC$FC + PlotDataFC$stderr , ymin=PlotDataFC$FC - PlotDataFC$stderr)
    dodge <- position_dodge(width=0.9)
    
    base_size <- 20
    # g1 <- ggplot(PlotDataFC,aes(x=factor(pid),y=FC,fill=factor(Set))) +
    PlotFC <- ggplot(PlotDataFC,aes(x=factor(Gene),y=FC,fill=factor(Set))) +  
      geom_bar(position="dodge",stat="identity") + 
      geom_errorbar(limits, position=dodge, width=0.4) +
      ylim(c(-2 ,2)) +
      # labs(x = "", y = "Log2-fold change",title="Comparison of Illumina data and RT-qPCR") +
      labs(x = "", y = "Log2-fold change") +
      scale_fill_manual(values = cols,name="") + 
      thememapBarplot(base_size = base_size,legend_key_size = 0.6) + 
      # theme( aspect.ratio = 9 / 16) + 
      theme(legend.position="bottom") 
    print(PlotFC)        
    ggsave("./PaperPlot/Microary_qPCR_barplot.pdf",device = 'pdf',width = 10,height = 6)
    
  ### bar Illumina expression ----------------------------------------------------------------------
    PlotDataEXP <- make.plot.data.exp.Ill.qPCR(qCPRdata = qCPRdataC, MeanFoldChangeClass = MeanFoldChangeClass,class="c")
    
    # print(PlotDataEXP)
    
    rns <- qgenes #levels(PlotDataFC$Gene)[c(2,5,6,1,3,4)] ; 
    PlotDataFC$Gene <- factor(PlotDataFC$Gene, levels = rns)
    # rns <- levels(PlotDataEXP$pid)[c(2,5,6,1,3,4)]  ; PlotDataEXP$pid <- factor(PlotDataEXP$pid, levels = rns)
    
    cols <-  RColorBrewer::brewer.pal(8,'Paired')[c(6,2)] 
    limits <- aes(ymax = PlotDataEXP$Mean + PlotDataEXP$stderr , ymin=PlotDataEXP$Mean - PlotDataEXP$stderr)
    dodge <- position_dodge(width = 0.9)
    
    base_size <- 20
    PlotExp <- ggplot(PlotDataEXP,aes(x=factor(Gene),y=Mean,fill=factor(Set))) +  
      geom_bar(position="dodge",stat="identity") + 
      geom_errorbar(limits, position=dodge, width=0.4) +
      ylim(c(0 ,10)) +
      scale_fill_manual(values = cols,name="") + 
      # labs(x = "", y = "Log2 mean expression",title="mean expression data of Illumina") +  
      labs(x = "", y = "Log2 mean expression") +  
      thememapBarplot(base_size = base_size,legend_key_size = 0.6) + 
      # theme( aspect.ratio = 9 / 16) + 
      theme(legend.position="bottom") 
    print(PlotExp) 
    ggsave("/Users/weinhol/GitHub/BGSC/PaperPlot/Microary_mean_barplot.pdf",device = 'pdf',width = 10,height = 6)
    
    gridExtra::grid.arrange(PlotExp  , PlotFC, ncol=1)
    ggsave("/Users/weinhol/GitHub/BGSC/PaperPlot/Microary_mean_Microary_qPCR__barplot.pdf",plot =  gridExtra::grid.arrange(PlotExp  , PlotFC, ncol=1) ,device = 'pdf',width = 10,height = 12)
    
  
  ### pheatmap of normData  ----------------------------------------------------------------------
    addBorder_gtable <- function(g){ g <- gtable_add_grob(g,
                                                          grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
                                                          t= 1 , b = nrow(g), l = 1, r = ncol(g))
      # g <- gtable_add_grob(g,
      #                      grobs = rectGrob(gp = gpar(fill = NA, lwd = 1)),
      #                      t = 1, l = 1, r = ncol(g))
      
      return(g)
    }
    
    Leset <- get.Lset()
    pmat <- t(matrix(c(0,1,0,0,0,1),2,3))
    dimnames(pmat) <- list( c('no RNAi','RNAi by siRNA ALL','RNAi by siRNA I'),c('no EGF','EGF'))  
    pmat <- pmat[c(1,3,2),] ## to get same oder as in paper
    bk = seq(0,max(1),by = 1)
    #hmcols <- c(RColorBrewer::brewer.pal(9,"Reds")[6],RColorBrewer::brewer.pal(9,"Blues")[7])# colorRampPalette(c("blue", "red"))(length(bk))
    hmcols <- RColorBrewer::brewer.pal(8,'Paired')[c(6,2)] 
    pheatmap::pheatmap(pmat,color = hmcols,breaks = seq(-1,1,by = 1) ,display_numbers = T,fontsize_number = 20,number_format = '%.0f',cluster_rows = FALSE,cluster_cols = FALSE,main='expression pattern group c',fontsize = 10,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black',legend = FALSE) 
    
    Lsets <- get.Lset()
    IndicatorVar <- lapply(Lsets,function(x){
      vec <- rep(0,6);names(vec) <- c(1:6)
      if (!is.null(x$s1)) { vec[x$s1] <- 1 }
      return(vec)
    })
    IndicatorVar.gtable <- purrr::map2( IndicatorVar , names(IndicatorVar), function(vec,.y){
      
      pmat <- t(matrix(vec,2,3))
      dimnames(pmat) <- list( c('no RNAi','RNAi by siRNA ALL','RNAi by siRNA I'),c('no EGF','EGF'))  
      pmat <- pmat[c(1,3,2),] ## to get same oder as in paper
      bk = seq(0,max(1),by = 1)
      #hmcols <- c(RColorBrewer::brewer.pal(9,"Reds")[6],RColorBrewer::brewer.pal(9,"Blues")[7])# colorRampPalette(c("blue", "red"))(length(bk))
      hmcols <- RColorBrewer::brewer.pal(8,'Paired')[c(6,2)] 
      
      if(.y == 'a') hmcols <- 'slategray'
      
      pheatmap::pheatmap(pmat,silent = T,color = hmcols,breaks = seq(-1,1,by = 1) ,display_numbers = T,fontsize_number = 16,number_format = '%.0f',cluster_rows = FALSE,cluster_cols = FALSE,main=paste0('Expression pattern of for group ', .y),fontsize = 12 ,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black',legend = FALSE) 
    } )
    
    grid.arrange(
      addBorder_gtable(IndicatorVar.gtable$a$gtable),
      addBorder_gtable(IndicatorVar.gtable$b$gtable),
      addBorder_gtable(IndicatorVar.gtable$c$gtable),
      addBorder_gtable(IndicatorVar.gtable$d$gtable),
      nrow=2,ncol=2
    )   
    ggsave("/Users/weinhol/GitHub/BGSC/PaperPlot/Schematic_ExpressionPattern.pdf",
           plot= grid.arrange(
                              addBorder_gtable(IndicatorVar.gtable$a$gtable),
                              addBorder_gtable(IndicatorVar.gtable$b$gtable),
                              addBorder_gtable(IndicatorVar.gtable$c$gtable),
                              addBorder_gtable(IndicatorVar.gtable$d$gtable),
                              nrow=2,ncol=2) 
           ,device = 'pdf',width = 11,height = 7)
    
    Gene.gtable <- list()
    for(tmoG in names(qgenesIDs)){
      pmat <- t(matrix(normData$E[qgenesIDs[[tmoG]],],2,3))
      dimnames(pmat) <- list( c('no RNAi','RNAi by siRNA ALL','RNAi by siRNA I'),c('no EGF','EGF'))  
      pmat <- pmat[c(1,3,2),] ## to get same oder as in paper
      
      pmat <- 2^pmat
      tmpnScore <-  ceiling(pmat/100)*100
      # tmpnScore <-  ceiling(pmat/10)*10
      
      bk = seq(0,max(tmpnScore),by = 1)
      hmcols <- colorRampPalette(c("white", "darkred"))(length(bk)-1)
      #pheatmap::pheatmap(pmat,color = hmcols,breaks = seq(0,max(tmpnScore),by = 1) ,display_numbers = T,fontsize_number = 12,cluster_rows = FALSE,cluster_cols = FALSE,main=tmoG,fontsize = 12,gaps_row = c(1,2),gaps_col = 1,border_color = 'black',number_color = 'black', )
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
      
   
    grid.arrange(
      addBorder_gtable(Gene.gtable$CKAP2L$gtable),
      addBorder_gtable(Gene.gtable$ROCK1$gtable),
      addBorder_gtable(Gene.gtable$TPR$gtable),
      addBorder_gtable(Gene.gtable$ALDH4A1$gtable),
      addBorder_gtable(Gene.gtable$CLCA2$gtable),
      addBorder_gtable(Gene.gtable$GALNS$gtable),
      nrow=2,ncol=3)                  

    
    
    
 
  ### bar log2 FC qPCR EGF  ----------------------------------------------------------------------
    qgenes <-  c('CKAP2L','ROCK1','TPR','ALDH4A1','CLCA2','GALNS')
    qgenes2 <- qgenes[qgenes!='CLCA2']
    
    data(qPCR_SF767_LNZ308, envir = environment()) 

    comparative_qPCR_CL_SFLNZ <- list()
    for(n in qgenes2){
      Ref_GAPDH <- as.double(qPCR_SF767_LNZ308[,'GAPDH'])
      tmp <-  as.double(qPCR_SF767_LNZ308[,n])
      comparative_qPCR_CL_SFLNZ[['GAPDH']][[n]] <- comparativeMethod_qPCR.analysis(Gene = tmp,Ref = Ref_GAPDH)
    } 
    CL2 <- c('SF767','SF767','LNZ308','LNZ308')
    Treat2 <- c("Co" ,"EGF","Co" ,"EGF")  
    CL_SFLNZ.log2FC <- comparativeMethod_qPCR.EGF.log2FC(comparative_qPCR_CL_SFLNZ[['GAPDH']],CL = CL2 ,Treat=Treat2,nREP=3)
    qCPRdataEGF <- CL_SFLNZ.log2FC$ddCT.EGF
    data.table::setkey(qCPRdataC,Gene)
    
    pdata <- data.frame(qCPRdataEGF)
    pdata$Gene <- as.factor(pdata$Gene)
    pdata$CellLine <- as.factor(pdata$CellLine)
    pdata$CellLine <- droplevels(pdata$CellLine)
    colnames(pdata)[c(7,8)] <- c("log2FC",'stderr')
    rns <- c('CKAP2L','ROCK1','TPR','ALDH4A1','GALNS'); pdata$Gene <- factor(pdata$Gene, levels = rns)
    rns <- c("SF767","LNZ308"); pdata$CellLine <- factor(pdata$CellLine, levels = rns)
    cols <- RColorBrewer::brewer.pal(11,"PRGn")[c(3,9)]
    SetCor <- cor(pdata$log2FC[c(1,3,5,7,9)],pdata$log2FC[c(2,4,6,8,10)])
    
    plotOUT <- '/Volumes/Macintosh HD/Users/weinhol/Promotion/ServerTMP/KapplerWichmann/PLOTS_SF767-1-799-6/'
    dodge <- position_dodge(width=0.9)
    base_size <- 20
    
    PlotFC <- ggplot(pdata,aes(x=factor(Gene),y=log2FC,fill=factor(CellLine),
                               ymin= log2FC - stderr, ymax= log2FC + stderr)) + 
      geom_bar(position="dodge",stat="identity") + 
      geom_errorbar(position=dodge, width=0.4)  +
      ylim(c(-2 ,2)) +
      # labs(x = "", y = "Log2-fold change",title=set) +
      labs(x = "", y = "Log2-fold change") +
      scale_fill_manual(values = cols,name="") + 
      thememapBarplot(base_size = base_size,legend_key_size = 0.6) + 
      theme(legend.position="bottom") 
    print(PlotFC)    
    ggsave(plot = PlotFC,filename ="./PaperPlot/qPCR_ddCT_log2FC_SF767_LNZ308_barplot.pdf",device = 'pdf',width = 10,height = 6)
    
    
    
    
  ### DES   ----------------------------------------------------------------------

      GeneExample <- c('ILMN_1687840','ILMN_1684585','ILMN_1730999','ILMN_2320964') 
      names(GeneExample) <- c('a','b','c','d')
      tmpPlot <- purrr::map2(GeneExample,names(GeneExample),function(.x,.y) Density.NV.fit.plot(id = .x ,normData, ALL.MUs, ALL.VARs, useGroup = .y ,DOplot = FALSE) )
      grid.arrange(tmpPlot$a + theme(legend.position = "none"),
                   tmpPlot$b + theme(legend.position = "none"),
                   tmpPlot$c + theme(legend.position = "none"),
                   tmpPlot$d + theme(legend.position = "none") ,ncol=2,nrow=2)
      
      
      Density.NV.fit.plot(id = 'ILMN_1730999' ,normData, ALL.MUs, ALL.VARs,useGroup = 'c' ,DOplot = FALSE,basesize = 20,GrAblack = TRUE,onlySYMBOL = TRUE )
      ggsave("/Users/weinhol/GitHub/BGSC/PaperPlot/DensityPlot_TPR.pdf",device = 'pdf',width = 10,height = 6)
      
      
      # id <- "ILMN_1730999"
      # tmpDensity <- lapply(qgenesIDs, function(id) Density.NV.fit.plot(id = id,normData,DOplot = TRUE) )
      # 
      # aaa <- Density.NV.fit.plot(id = 'ILMN_1687840',normData, ALL.MUs, ALL.VARs,useGroup = "a",DOplot = TRUE)
      # bbb <- Density.NV.fit.plot(id = 'ILMN_1684585',normData, ALL.MUs, ALL.VARs, useGroup = "b",DOplot = TRUE)
      # ccc <- Density.NV.fit.plot(id = 'ILMN_1730999',normData, ALL.MUs, ALL.VARs, useGroup = "c",DOplot = TRUE)
      # ddd <- Density.NV.fit.plot(id = 'ILMN_2320964',normData, ALL.MUs, ALL.VARs, useGroup = "d",DOplot = TRUE)
      # 
      
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
  
  ### DAVID -> /Users/weinhol/GitHub/BGSC/perl/CallDAVID.txt
  #
  #/Users/weinhol/miniconda2/bin/perl /Users/weinhol/GitHub/BGSC/perl/chartReport_modify2.pl /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_a.txt 'REFSEQ_MRNA' /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_a_chartReport.txt
  #/Users/weinhol/miniconda2/bin/perl /Users/weinhol/GitHub/BGSC/perl/chartReport_modify2.pl /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_b.txt 'REFSEQ_MRNA' /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_b_chartReport.txt
  #/Users/weinhol/miniconda2/bin/perl /Users/weinhol/GitHub/BGSC/perl/chartReport_modify2.pl /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_c.txt 'REFSEQ_MRNA' /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_c_chartReport.txt
  #/Users/weinhol/miniconda2/bin/perl /Users/weinhol/GitHub/BGSC/perl/chartReport_modify2.pl /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_d.txt 'REFSEQ_MRNA' /Users/weinhol/GitHub/BGSC/inst/extdata/Filtering_REFSEQ_MRNA_d_chartReport.txt
  ##
  
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
      
      # ChartReportCategorySigList[[ChartReportCategory]][[i]] <-ChartReportCategoryList[[ChartReportCategory]][[i]][which(Benjamini < 0.1),]
      ChartReportCategorySigList[[ChartReportCategory]][[i]] <-ChartReportCategoryList[[ChartReportCategory]][[i]][which(FDR < (0.05 * 100)),] ### FDR is multiplied with 100
      
      setkey(ChartReportCategorySigList[[ChartReportCategory]][[i]],'Term')
      
      ChartReportCategoryTermsList[[ChartReportCategory]] <- unique(c(ChartReportCategoryTermsList[[ChartReportCategory]], ChartReportCategorySigList[[ChartReportCategory]][[i]]$Term) )
      
    }
  }
  
  do.call(rbind, purrr::map2(ChartReportCategorySigList$KEGG_PATHWAY,names(ChartReportCategorySigList$KEGG_PATHWAY), 
                            function(x,group) {
                                       tmp <- x[,c('Category','Term', 'Count','%','Pvalue','FDR')]; 
                                       if (nrow(tmp) > 1) { tmp <- cbind(tmp,group) } else {tmp$group <- NA}
                                       tmp
                            })
         )

 do.call(rbind, purrr::map2(ChartReportCategorySigList$GOTERM_BP_FAT,names(ChartReportCategorySigList$KEGG_PATHWAY), 
                            function(x,group) {
                              tmp <- x[,c('Category','Term', 'Count','%','Pvalue','FDR')]; 
                              if (nrow(tmp) > 1) { tmp <- cbind(tmp,group) } else {tmp$group <- NA}
                              tmp
                            })
 )
   
  
    ChartReportCategorySigMatList <- list()
    # ChartReportCategory <- 'KEGG_PATHWAY' #ChartReportCategory <- 'GOTERM_MF_FAT'
    for(ChartReportCategory in names(ChartReportCategorySigList)){
      print(ChartReportCategory)
      ChartReportCategorySigMat <- matrix(0, length(names(ChartReportCategorySigList[[ChartReportCategory]])), length(ChartReportCategoryTermsList[[ChartReportCategory]]),
                                          dimnames = list(names(ChartReportCategorySigList[[ChartReportCategory]]),ChartReportCategoryTermsList[[ChartReportCategory]]) )
      print(dim(ChartReportCategorySigMat))
      if(dim(ChartReportCategorySigMat)[2] > 1){
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
        
        tmp <- do.call(rbind, purrr::map2(ChartReportCategorySigList[[ChartReportCategory]],names(ChartReportCategorySigList$KEGG_PATHWAY), 
                                          function(x,group) {
                                            tmp <- x[,c('Category','Term', 'Count','%','Pvalue','FDR', "List Total","Pop Hits","Pop Total","Fold Enrichment")]; 
                                            if (nrow(tmp) > 0) { tmp <- cbind(tmp,group) } else {tmp$group <- NA}
                                            tmp
                                          })
        )
        print(tmp)
        write.csv2(tmp,paste0(enrichmentFolder,'/ChartReportSig_',ChartReportCategory,'.csv'))
        
        
        
      } 
      
      ChartReportCategorySigMatList[[ChartReportCategory]] <- ChartReportCategorySigMat
      
  
      
    }
    

  
  #######################################
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