main <- function(){
  
  # normalize Data
  normData <- normalizeExpData()
  
  # calulate logLik for class  Data
  Lsets <- get.Lset()
  normDataLogLik <- logLikelihoodOfnormData(normData$E)
  
  # calulate BIC from logLik 
  npar <- sapply(Lsets, function(x) sum(!sapply(x,is.null ) )) + 1  ## number parapeters for LogLilk -> mean + var 
  k <-  sapply(Lsets, function(x) sum( sapply(x,length) ))
  normDataBIC <- get.IC(normDataLogLik , npar, k , IC = 'BIC')
  BICminInd <- apply( normDataBIC,MARGIN = 1, FUN = minIndex)
  print(table(BICminInd))
  
  # calulate Posterior from BIC
  normDataPosterior <- get.Posterior( normDataBIC ,Pis = c(0.7,0.1,0.1,0.1))
  POSTmaxInd <- apply(normDataPosterior, MARGIN = 1 ,FUN = maxIndex)
  print(table(POSTmaxInd))
  
  PostClass <- get.gene.classes(data = normDataPosterior,indexing = "max",filter = 0.75, DoPlot = TRUE)
  
  
  ### enrichment analysis
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