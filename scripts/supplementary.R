plotOUT <- './PaperPlot/'
dodge <- position_dodge(width=0.9)
base_size <- 20

# PlotsPaper mRNA G0G1 Proliferation  -------------------------------------
ExcelFiles <- c("./rawdata/supplementary/RNA_Daten_zur_SF767_CCND1.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_EGFRall.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_EGFRv1.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_EGFRv4.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_G0G1.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_GAPDH.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_MMP2.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_proliferation48H.csv",
                "./rawdata/supplementary/RNA_Daten_zur_SF767_SVV.csv")
ExcelFilesNames <- c("CCND1","EGFRall","EGFRv1","EGFRv4","G0G1","GAPDH","MMP2","proliferation48H","SVV")
names(ExcelFiles) <- ExcelFilesNames
ExcelFilesList <- lapply(ExcelFilesNames,function(x) read.csv(ExcelFiles[x],sep = ';',dec = ',') )
names(ExcelFilesList) <-  ExcelFilesNames
OutNames <- c("CCND1","EGFR variants I-IV ","EGFR variant I","EGFR variant IV","G0G1","GAPDH","MMP2","proliferation48H","BIRC5")
names(OutNames) <- ExcelFilesNames

# PlotsPaper mRNA -------------------------------------
  mRNANames <- c("EGFRall","EGFRv1","EGFRv4","GAPDH","MMP2","SVV","CCND1")
  PData <- data.frame()
  for(mRNA in mRNANames){
    mat <- matrix(as.double(unlist(ExcelFilesList[[mRNA]][,c(5:8)])),4,4)
    dimnames(mat) <- list(as.character(ExcelFilesList[[mRNA]][,'Behandlung']),c('rep1','rep2','rep3','rep4'))
    # rownames(mat)[2] <- 'siRNA-luci'
    rownames(mat)[2] <- 'siRNA-nonsense'
    print(mat)
    pdata <- data.frame( "Treat" = rownames(mat),
                         "Gene" = mRNA,
                         "Out" = OutNames[mRNA],
                         "GeoMean" = (2^rowMeans(log2(mat),na.rm = T )),
                         "Mean" = rowMeans(mat,na.rm = T),
                         "stderr" = apply(mat,1,function(row) plotrix::std.error((row),na.rm = T)),
                         "sd" = apply(mat,1,function(row) sd((row),na.rm = T)))
    
    PData <-  rbind(PData,pdata)
    
    rns <- levels(pdata$Treat)[c(1,4,2,3)] ; pdata$Treat <- factor(pdata$Treat, levels = rns)
    
    limitsMeanErr <- aes(ymax = pdata$Mean + pdata$stderr , ymin=pdata$Mean - pdata$stderr) ; limitsMeanErrMAX <- round(max(pdata$Mean + pdata$stderr)+10)
    limitsGeoMeanErr <- aes(ymax = pdata$GeoMean + pdata$stderr , ymin=pdata$GeoMean - pdata$stderr); limitsGeoMeanErrMAX <- round(max(pdata$GeoMean + pdata$stderr)+10)
    limitsMeanSd <- aes(ymax = pdata$Mean + pdata$sd , ymin=pdata$Mean - pdata$sd); limitsMeanSdMAX <- round(max(pdata$Mean + pdata$sd)+10)
    limitsGeoMeanSd <- aes(ymax = pdata$GeoMean + pdata$sd , ymin=pdata$GeoMean - pdata$sd); limitsGeoMeanSdMAX <- round(max(pdata$GeoMean + pdata$sd)+10)
    
    PlotMeanSd <- ggplot(pdata,aes(x=factor(Treat),y=Mean)) +  
      geom_bar(position="dodge",stat="identity") + 
      geom_errorbar(limitsMeanSd, position=dodge, width=0.4) +
      #ylim(c(0 ,limitsMeanSdMAX)) +
      labs(x = "", y = paste0("relative mRNA level of ",OutNames[mRNA]," copies / \n HPRT copies"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
      thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
      scale_y_continuous(breaks=seq(0,limitsMeanSdMAX,20))
    #theme(legend.position="bottom") 
    # scale_fill_manual(values = cols,name="") + 
    print(PlotMeanSd)
    ggsave(plot = PlotMeanSd,filename = paste0(plotOUT,'mRNA_',mRNA,'_MeanSd.pdf'),width = 12,height = 8)
    
    # PlotMeanErr <- ggplot(pdata,aes(x=factor(Treat),y=Mean)) +  
    #   geom_bar(position="dodge",stat="identity") + 
    #   geom_errorbar(limitsMeanErr, position=dodge, width=0.4) +
    #   # ylim(c(0 ,limitsMeanErrMAX)) +
    #   labs(x = "", y = paste0("relative mRNA level of ",OutNames[mRNA]," copies / \n HPRT copies"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
    #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
    #   scale_y_continuous(breaks=seq(0,limitsMeanErrMAX,20))
    # #theme(legend.position="bottom") 
    # # scale_fill_manual(values = cols,name="") + 
    # print(PlotMeanErr)
    # ggsave(plot = PlotMeanErr,filename = paste0(plotOUT,'mRNA_',mRNA,'_MeanErr.pdf'),width = 12,height = 8)
    # 
    # PlotGeoMeanSd <- ggplot(pdata,aes(x=factor(Treat),y=GeoMean)) +  
    #   geom_bar(position="dodge",stat="identity") + 
    #   geom_errorbar(limitsGeoMeanSd, position=dodge, width=0.4) +
    #   #ylim(c(0 ,limitsGeoMeanSdMAX)) +
    #   labs(x = "", y = paste0("relative mRNA level of ",OutNames[mRNA]," copies / \n HPRT copies"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
    #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
    #   scale_y_continuous(breaks=seq(0,limitsGeoMeanSdMAX,20))
    # #theme(legend.position="bottom") 
    # # scale_fill_manual(values = cols,name="") + 
    # print(PlotGeoMeanSd)
    # ggsave(plot = PlotGeoMeanSd,filename = paste0(plotOUT,'mRNA_',mRNA,'_GeoMeanSd.pdf'),width = 12,height = 8)
    # 
    # PlotGeoMeanErr <- ggplot(pdata,aes(x=factor(Treat),y=GeoMean)) +  
    #   geom_bar(position="dodge",stat="identity") + 
    #   geom_errorbar(limitsGeoMeanErr, position=dodge, width=0.4) +
    #   #ylim(c(0 ,limitsGeoMeanErrMAX)) +
    #   labs(x = "", y = paste0("relative mRNA level of ",OutNames[mRNA]," copies / \n HPRT copies"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
    #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
    #   scale_y_continuous(breaks=seq(0,limitsGeoMeanErrMAX,20))
    # #theme(legend.position="bottom") 
    # # scale_fill_manual(values = cols,name="") + 
    # print(PlotGeoMeanErr)
    # ggsave(plot = PlotGeoMeanErr,filename = paste0(plotOUT,'mRNA_',mRNA,'_GeoMeanErr.pdf'),width = 12,height = 8)
  }  
  
  rns <- levels(PData$Treat)[c(1,4,2,3)] ; PData$Treat <- factor(PData$Treat, levels = rns)
  limitsMeanErr <- aes(ymax = PData$Mean + PData$stderr , ymin=PData$Mean - PData$stderr) ; limitsMeanErrMAX <- round(max(PData$Mean + PData$stderr)+10)
  limitsGeoMeanErr <- aes(ymax = PData$GeoMean + PData$stderr , ymin=PData$GeoMean - PData$stderr); limitsGeoMeanErrMAX <- round(max(PData$GeoMean + PData$stderr)+10)
  limitsMeanSd <- aes(ymax = PData$Mean + PData$sd , ymin=PData$Mean - PData$sd); limitsMeanSdMAX <- round(max(PData$Mean + PData$sd)+10)
  limitsGeoMeanSd <- aes(ymax = PData$GeoMean + PData$sd , ymin=PData$GeoMean - PData$sd); limitsGeoMeanSdMAX <- round(max(PData$GeoMean + PData$sd)+10)
  
  PlotMeanSd <- ggplot(PData,aes(x=factor(Treat),y=Mean)) +  
    geom_bar(position="dodge",stat="identity") + 
    geom_errorbar(limitsMeanSd, position=dodge, width=0.4) +
    #ylim(c(0 ,limitsMeanSdMAX)) +
    labs(x = "", y = paste0("relative mRNA level of gene copies / HPRT copies")) +
    thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
    scale_y_continuous(breaks=seq(0,limitsMeanSdMAX,20)) +
    facet_grid(Out ~ .)
  #theme(legend.position="bottom") 
  # scale_fill_manual(values = cols,name="") + 
  print(PlotMeanSd)
  ggsave(plot = PlotMeanSd,filename = paste0(plotOUT,'mRNA_ALL_MeanSd.pdf'),width = 12,height = 20)
  
  
  barCol <- RColorBrewer::brewer.pal(11,'Spectral')[6]
  
  for(i in c("EGFR",'Control','Proli','EGFRControl')){
    set <- switch(i,
           "EGFR" = as.character( levels(PData$Gene)[1:3] ),
           'Control' = as.character( levels(PData$Gene)[4:5] ),
           'Proli' = as.character( levels(PData$Gene)[6:7] ),
           'EGFRControl'=  as.character( levels(PData$Gene)[1:5] )
    )
    tmp <- PData[PData$Gene %in% set,]
    tmp_limitsMeanSd <- aes(ymax = tmp$Mean + tmp$sd , ymin=tmp$Mean - tmp$sd); limitsMeanSdMAX <- round(max(tmp$Mean + tmp$sd)+10)
    
    PlotMeanSd <- ggplot(tmp,aes(x=factor(Treat),y=Mean)) +  
      geom_bar(position="dodge",stat="identity",fill=barCol,color='black') +
      geom_errorbar(tmp_limitsMeanSd, position=dodge, width=0.4) +
      labs(x = "", y = paste0("relative mRNA level of gene copies / HPRT copies")) +
      thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
      scale_y_continuous(breaks=seq(0,limitsMeanSdMAX,20)) +
      facet_grid(Out ~ .)
    print(PlotMeanSd)
    ggsave(plot = PlotMeanSd,filename = paste0(plotOUT,'mRNA_ALL_MeanSd_',i,'.pdf'),width = 8,height = 15)
  }
  
  # PlotMeanErr <- ggplot(PData,aes(x=factor(Treat),y=Mean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsMeanErr, position=dodge, width=0.4) +
  #   # ylim(c(0 ,limitsMeanErrMAX)) +
  #   labs(x = "", y = paste0("relative mRNA level of gene copies / HPRT copies")) +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
  #   scale_y_continuous(breaks=seq(0,limitsMeanErrMAX,20)) +
  #   facet_grid(Out ~ .)
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotMeanErr)
  # ggsave(plot = PlotMeanErr,filename = paste0(plotOUT,'mRNA_ALL_MeanErr.pdf'),width = 12,height = 20)
  # 
  # PlotGeoMeanSd <- ggplot(PData,aes(x=factor(Treat),y=GeoMean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsGeoMeanSd, position=dodge, width=0.4) +
  #   #ylim(c(0 ,limitsGeoMeanSdMAX)) +
  #   labs(x = "", y = paste0("relative mRNA level of gene copies / HPRT copies")) +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
  #   scale_y_continuous(breaks=seq(0,limitsGeoMeanSdMAX,20)) +
  #   facet_grid(Out ~ .)
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotGeoMeanSd)
  # ggsave(plot = PlotGeoMeanSd,filename = paste0(plotOUT,'mRNA_ALL_GeoMeanSd.pdf'),width = 12,height = 20)
  # 
  # PlotGeoMeanErr <- ggplot(PData,aes(x=factor(Treat),y=GeoMean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsGeoMeanErr, position=dodge, width=0.4) +
  #   #ylim(c(0 ,limitsGeoMeanErrMAX)) +
  #   labs(x = "", y = paste0("relative mRNA level of gene copies / HPRT copies")) +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) +
  #   scale_y_continuous(breaks=seq(0,limitsGeoMeanErrMAX,20)) +
  #   facet_grid(Out ~ .)
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotGeoMeanErr)
  # ggsave(plot = PlotGeoMeanErr,filename = paste0(plotOUT,'mRNA_ALL_GeoMeanErr.pdf'),width = 12,height = 20)

# PlotsPaper proliferation48H -------------------------------------
  set <- "proliferation48H"
  mat <- matrix(as.double(unlist(ExcelFilesList[[set]][,c(3:7)])),10,5)
  dimnames(mat) <- list(as.character(ExcelFilesList[[set]][,'Behandlung']),c('rep1','rep2','rep3','rep4','rep5'))
  
  mat <-  mat[c("K","Lu","All",'Full'),]
  # rownames(mat) <- c('control','siRNA-luci','siRNA-ALL','siRNA-I')
  rownames(mat) <- c('control','siRNA-nonsense','siRNA-ALL','siRNA-I')
  
  pdata <- data.frame( "Treat" = rownames(mat),
                       # "Out" = OutNames[mRNA],
                       "GeoMean" = (2^rowMeans(log2(mat),na.rm = T )),
                       "Mean" = rowMeans(mat,na.rm = T),
                       "stderr" = apply(mat,1,function(row) plotrix::std.error((row),na.rm = T)),
                       "sd" = apply(mat,1,function(row) sd((row),na.rm = T)))
  rns <- levels(pdata$Treat)[c(1,4,2,3)] ; pdata$Treat <- factor(pdata$Treat, levels = rns)
  
  limitsMeanErr <- aes(ymax = pdata$Mean + pdata$stderr , ymin=pdata$Mean - pdata$stderr) ; limitsMeanErrMAX <- round(max(pdata$Mean + pdata$stderr)+10)
  limitsGeoMeanErr <- aes(ymax = pdata$GeoMean + pdata$stderr , ymin=pdata$GeoMean - pdata$stderr); limitsGeoMeanErrMAX <- round(max(pdata$GeoMean + pdata$stderr)+10)
  limitsMeanSd <- aes(ymax = pdata$Mean + pdata$sd , ymin=pdata$Mean - pdata$sd); limitsMeanSdMAX <- round(max(pdata$Mean + pdata$sd)+10)
  limitsGeoMeanSd <- aes(ymax = pdata$GeoMean + pdata$sd , ymin=pdata$GeoMean - pdata$sd); limitsGeoMeanSdMAX <- round(max(pdata$GeoMean + pdata$sd)+10)
  
  PlotMeanSd <- ggplot(pdata,aes(x=factor(Treat),y=Mean)) +  
    geom_bar(position="dodge",stat="identity",fill=barCol,color='black') + 
    geom_errorbar(limitsMeanSd, position=dodge, width=0.4) +
    ylim(c(0 ,limitsMeanSdMAX)) +
    # scale_y_continuous(breaks=seq(0,limitsMeanSdMAX,20)) +
    labs(x = "", y = paste0("cell number after 48 hours \n of proliferation"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
    thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  #theme(legend.position="bottom") 
  # scale_fill_manual(values = cols,name="") + 
  print(PlotMeanSd)
  ggsave(plot = PlotMeanSd,filename = paste0(plotOUT,'proliferation48H_MeanSd.pdf'),width = 12,height = 8)
  
  # PlotMeanErr <- ggplot(pdata,aes(x=factor(Treat),y=Mean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsMeanErr, position=dodge, width=0.4) +
  #   ylim(c(0 ,limitsMeanErrMAX)) +
  #   labs(x = "", y = paste0("cell number after 48 hours of proliferation"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotMeanErr)
  # ggsave(plot = PlotMeanErr,filename = paste0(plotOUT,'proliferation48H_MeanErr.pdf'),width = 12,height = 8)
  # 
  # PlotGeoMeanSd <- ggplot(pdata,aes(x=factor(Treat),y=GeoMean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsGeoMeanSd, position=dodge, width=0.4) +
  #   ylim(c(0 ,limitsGeoMeanSdMAX)) +
  #   labs(x = "", y = paste0("cell number after 48 hours of proliferation"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotGeoMeanSd)
  # ggsave(plot = PlotGeoMeanSd,filename = paste0(plotOUT,'proliferation48H_GeoMeanSd.pdf'),width = 12,height = 8)
  # 
  # PlotGeoMeanErr <- ggplot(pdata,aes(x=factor(Treat),y=GeoMean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsGeoMeanErr, position=dodge, width=0.4) +
  #   ylim(c(0 ,limitsGeoMeanErrMAX)) +
  #   labs(x = "", y = paste0("cell number after 48 hours of proliferation"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotGeoMeanErr)
  # ggsave(plot = PlotGeoMeanErr,filename = paste0(plotOUT,'proliferation48H_GeoMeanErr.pdf'),width = 12,height = 8)

# PlotsPaper G0G1 -------------------------------------
  set <- "G0G1"
  mat <- matrix(as.double(unlist(ExcelFilesList[[set]][,c(3:8)])),10,6)
  dimnames(mat) <- list(as.character(ExcelFilesList[[set]][,'Bezeichnung']),c('rep1','rep2','rep3','rep4','rep5','rep6'))
  
  mat <-  mat[c("Ko","LU","All",'Full'),]
  rownames(mat) <- c('control','siRNA-luci','siRNA-ALL','siRNA-I')
  rownames(mat) <- c('control','siRNA-nonsense','siRNA-ALL','siRNA-I')
  pdata <- data.frame( "Treat" = rownames(mat),
                       # "Out" = OutNames[mRNA],
                       "GeoMean" = (2^rowMeans(log2(mat),na.rm = T )),
                       "Mean" = rowMeans(mat,na.rm = T),
                       "stderr" = apply(mat,1,function(row) plotrix::std.error((row),na.rm = T)),
                       "sd" = apply(mat,1,function(row) sd((row),na.rm = T)))
  rns <- levels(pdata$Treat)[c(1,4,2,3)] ; pdata$Treat <- factor(pdata$Treat, levels = rns)
  
  limitsMeanErr <- aes(ymax = pdata$Mean + pdata$stderr , ymin=pdata$Mean - pdata$stderr) ; limitsMeanErrMAX <- round(max(pdata$Mean + pdata$stderr)+10)
  limitsGeoMeanErr <- aes(ymax = pdata$GeoMean + pdata$stderr , ymin=pdata$GeoMean - pdata$stderr); limitsGeoMeanErrMAX <- round(max(pdata$GeoMean + pdata$stderr)+10)
  limitsMeanSd <- aes(ymax = pdata$Mean + pdata$sd , ymin=pdata$Mean - pdata$sd); limitsMeanSdMAX <- round(max(pdata$Mean + pdata$sd)+10)
  limitsGeoMeanSd <- aes(ymax = pdata$GeoMean + pdata$sd , ymin=pdata$GeoMean - pdata$sd); limitsGeoMeanSdMAX <- round(max(pdata$GeoMean + pdata$sd)+10)
  
  PlotMeanSd <- ggplot(pdata,aes(x=factor(Treat),y=Mean)) +  
    geom_bar(position="dodge",stat="identity",fill=barCol,color='black') + 
    geom_errorbar(limitsMeanSd, position=dodge, width=0.4) +
    # ylim(c(0 ,limitsMeanSdMAX)) +
    scale_y_continuous(breaks=seq(0,limitsMeanSdMAX,10)) +
    labs(x = "", y = paste0("relative level of G0/G1 cells"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
    thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  print(PlotMeanSd)
  ggsave(plot = PlotMeanSd,filename = paste0(plotOUT,'G0G1_MeanSd.pdf'),width = 12,height = 8)
  
  # PlotMeanErr <- ggplot(pdata,aes(x=factor(Treat),y=Mean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsMeanErr, position=dodge, width=0.4) +
  #   scale_y_continuous(breaks=seq(0,limitsMeanErrMAX,10)) +
  #   labs(x = "", y = paste0("relative level of G0/G1 cells"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotMeanErr)
  # ggsave(plot = PlotMeanErr,filename = paste0(plotOUT,'G0G1_MeanErr.pdf'),width = 12,height = 8)
  # 
  # PlotGeoMeanSd <- ggplot(pdata,aes(x=factor(Treat),y=GeoMean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsGeoMeanSd, position=dodge, width=0.4) +
  #   scale_y_continuous(breaks=seq(0,limitsGeoMeanSdMAX,10)) +
  #   labs(x = "", y = paste0("relative level of G0/G1 cells"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotGeoMeanSd)
  # ggsave(plot = PlotGeoMeanSd,filename = paste0(plotOUT,'G0G1_GeoMeanSd.pdf'),width = 12,height = 8)
  # 
  # PlotGeoMeanErr <- ggplot(pdata,aes(x=factor(Treat),y=GeoMean)) +  
  #   geom_bar(position="dodge",stat="identity") + 
  #   geom_errorbar(limitsGeoMeanErr, position=dodge, width=0.4) +
  #   scale_y_continuous(breaks=seq(0,limitsGeoMeanErrMAX,10)) +
  #   labs(x = "", y = paste0("relative level of G0/G1 cells"))+ # ,title="PAPER PLOT - Comparison of SF767 - Illumina data and RT-qPCR") +
  #   thememapBarplot(base_size = base_size,legend_key_size = 0.6) 
  # #theme(legend.position="bottom") 
  # # scale_fill_manual(values = cols,name="") + 
  # print(PlotGeoMeanErr)
  # ggsave(plot = PlotGeoMeanErr,filename = paste0(plotOUT,'G0G1_GeoMeanErr.pdf'),width = 12,height = 8)
