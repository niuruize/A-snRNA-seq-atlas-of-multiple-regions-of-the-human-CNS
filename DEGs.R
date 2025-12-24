
library(MAST);library(Seurat);library(dplyr)
##整理表达矩阵-构建SingleCellAssay对象
#Idents(scRNA)="celltype"
#table(Idents(scRNA))
#scRNA = subset(scRNA,downsample=20000)
#scRNA$GenesPerUMI <- log10(scRNA$nFeature_RNA)/log10(scRNA$nCount_RNA)

#setwd("/home/niuruize/Aging/scRNA/MAST")
load("/home/niuruize2/aging/DS/PFC_STP/DS_PFC_STP_seurat.RData")
table(scRNA$celltype);table(scRNA$Group);table(scRNA$celltype,scRNA$Group)
table(scRNA$Region)
scRNA <- scRNA[,scRNA$Region %in% c("PFC")]
table(scRNA$celltype,scRNA$Group)
scRNA1 = scRNA
celltypes = c("ExN","Astro","Micro","MOL","InN","OPC","Epend","Endo")

for (celltype in celltypes)
{
  Idents(scRNA1) <- 'Group'; table(Idents(scRNA1))
  scRNA2  <- scRNA1[,scRNA1$celltype %in% celltypes]
  scRNA2 <- subset(scRNA2, downsample=10000)
  fData = data.frame(symbolid=rownames(scRNA2),primerid=rownames(scRNA2))
  rownames(fData)=fData$symbolid
  cData = scRNA2@meta.data
  cData$wellKey <- rownames(cData)
  sca = FromMatrix(as.matrix(scRNA2@assays$RNA@data), cData = cData,fData = fData,check_sanity = FALSE)
  rm(scRNA2)
  gc()
  dim(sca)
  table(colData(sca)$Group)
  cond<-factor(colData(sca)$Group)
  cond<-relevel(cond,"Adult")
  colData(sca)$condition<-cond
  
  colData(sca)$Sex=factor(colData(sca)$Sex)
  #colData(sca)$Age=as.numeric(colData(sca)$Age)
  colData(sca)$percent.mt=as.numeric(colData(sca)$percent.mt)
  #colData(sca)$percent.rb=as.numeric(colData(sca)$percent.rb)
  colData(sca)$nCount_RNA=as.numeric(colData(sca)$nCount_RNA)
  #colData(sca)$datasets=factor(colData(sca)$datasets)
  
  ##-----
  ## (1) 校正cngeneson协变量:默认参数
  zlmCond <- zlm(~condition + nCount_RNA + percent.mt + Sex, sca, method="bayesglm", ebayes=TRUE)
  summaryCond <- summary(zlmCond,doLRT='conditionAging')
  summaryDt <- summaryCond$datatable
  levels(summaryDt$contrast)
  
  #整理结果
  df_pval = summaryDt %>% 
    dplyr::filter(contrast=='conditionAging') %>% 
    dplyr::filter(component=='H') %>% 
    dplyr::select(primerid, `Pr(>Chisq)`)
  
  df_logfc = summaryDt %>% 
    dplyr::filter(contrast=='conditionAging') %>% 
    dplyr::filter(component=='logFC') %>% 
    dplyr::select(primerid, coef, ci.hi, ci.lo)
  
  df_stat = dplyr::inner_join(df_logfc, df_pval) %>% 
    dplyr::rename("symbol"="primerid") %>% 
    dplyr::rename("pval"="Pr(>Chisq)","logFC"="coef") %>% 
    dplyr::mutate("fdr" = p.adjust(pval)) %>% 
    dplyr::arrange(fdr)
  head(df_stat)
  
  df_stat$FC<-10^(abs(df_stat$logFC))
  df_stat$FC<-ifelse(df_stat$logFC>0,df_stat$FC*(1),df_stat$FC*-1) 
  df_stat$log2FC <- log2(abs(df_stat$FC))
  df_stat$log2FC <- ifelse(df_stat$FC>0,df_stat$log2FC*(1),df_stat$log2FC*-1)
  
  write.csv(df_stat, file = paste0('PFC_',celltype,"_Aging_Adult.csv"))
}

