#加载R包
library(Seurat);library(tidyverse);library(cowplot);library(patchwork);library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db);library(MySeuratWrappers);library(scRNAtoolVis);library(readxl);library(harmony)
#remotes::install_github("lyc-1995/MySeuratWrappers")
#devtools::install_github('junjunlab/scRNAtoolVis')
#‘geomtextpath’, ‘ggh4x’, ‘ggunchull’, ‘jjAnno’

#load data
load("~/Aging/Endo/ACC_Endo.RData")
load("~/Aging/Endo/AMY_Endo.RData")
load("~/Aging/Endo/EC_Endo.RData")
load("~/Aging/Endo/HIP_Endo.RData")
load("~/Aging/Endo/MB_Endo.RData")
load("~/Aging/Endo/PFC_Endo.RData")
load("~/Aging/Endo/Retina_Endo.RData")
load("~/Aging/Endo/SC_Endo.RData")

#合并数据
scRNA = merge(ACC_Endo, y=c(AMY_Endo,EC_Endo,HIP_Endo,MB_Endo,PFC_Endo,Retina_Endo,SC_Endo))
rm(ACC_Endo,AMY_Endo,EC_Endo,HIP_Endo,MB_Endo,PFC_Endo,Retina_Endo,SC_Endo)
save(scRNA,file="Endo_merge.RData")

#整合数据
###--数据标准化
scRNA <- SCTransform(scRNA)
###--PCA
scRNA <- RunPCA(scRNA, npcs= 30, verbose=FALSE)
#ElbowPlot(scRNA, ndims = 30)
###--数据整合-group.by.vars参数是设置按哪个分组整合
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=20)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:20) %>% RunUMAP(reduction="harmony", dims = 1:20)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:20) %>% FindClusters(dims = 1:20, resolution = 0.4)
save(scRNA,file="Endo.RData")



