#加载R包
library(Seurat);library(tidyverse);library(cowplot);library(patchwork);library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db);library(MySeuratWrappers);library(scRNAtoolVis);library(readxl);library(harmony)
#remotes::install_github("lyc-1995/MySeuratWrappers")
#devtools::install_github('junjunlab/scRNAtoolVis')
#‘geomtextpath’, ‘ggh4x’, ‘ggunchull’, ‘jjAnno’

#load data
load("~/Aging/OPC/ACC_OPC.RData")
load("~/Aging/OPC/AMY_OPC.RData")
load("~/Aging/OPC/EC_OPC.RData")
load("~/Aging/OPC/HIP_OPC.RData")
load("~/Aging/OPC/MB_OPC.RData")
load("~/Aging/OPC/PFC_OPC.RData")
load("~/Aging/OPC/SC_OPC.RData")

#合并数据
scRNA = merge(ACC_OPC, y=c(AMY_OPC,EC_OPC,HIP_OPC,MB_OPC,PFC_OPC,SC_OPC))
rm(ACC_OPC,AMY_OPC,EC_OPC,HIP_OPC,MB_OPC,PFC_OPC,SC_OPC)
save(scRNA,file="OPC_merge.RData")

#整合数据
###--数据标准化
scRNA <- SCTransform(scRNA)
###--PCA
scRNA <- RunPCA(scRNA, npcs= 30, verbose=FALSE)
#ElbowPlot(scRNA, ndims = 30)
###--数据整合-group.by.vars参数是设置按哪个分组整合
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=30)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:30) %>% RunUMAP(reduction="harmony", dims = 1:30)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:30) %>% FindClusters(dims = 1:30, resolution = 0.8)
save(scRNA,file="OPC.RData")



