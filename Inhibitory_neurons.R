#加载R包
library(Seurat);library(tidyverse);library(cowplot);library(patchwork);library(ggplot2);library(limma);library(AnnotationDbi);library(org.Hs.eg.db);library(MySeuratWrappers);library(scRNAtoolVis);library(readxl);library(harmony);library(writexl)
#remotes::install_github("lyc-1995/MySeuratWrappers")
#devtools::install_github('junjunlab/scRNAtoolVis')
#‘geomtextpath’, ‘ggh4x’, ‘ggunchull’, ‘jjAnno’

#load data
load("~/Aging/InN/ACC_InN.RData")
load("~/Aging/InN/AMY_InN.RData")
load("~/Aging/InN/EC_InN.RData")
load("~/Aging/InN/HIP_InN.RData")
load("~/Aging/InN/MB_InN.RData")
load("~/Aging/InN/PFC_InN.RData")
load("~/Aging/InN/SC_InN.RData")

#合并数据
scRNA = merge(ACC_InN, y=c(AMY_InN,EC_InN,HIP_InN,MB_InN,PFC_InN,SC_InN))
rm(ACC_InN,AMY_InN,EC_InN,HIP_InN,MB_InN,PFC_InN,SC_InN)
save(scRNA,file="InN_merge.RData")

#整合数据
###--数据标准化
scRNA <- SCTransform(scRNA)
###--PCA
scRNA <- RunPCA(scRNA, npcs= 30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
###--数据整合-group.by.vars参数是设置按哪个分组整合
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=30)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:30) %>% RunUMAP(reduction="harmony", dims = 1:30)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:30) %>% FindClusters(dims = 1:30, resolution = 0.4)
save(scRNA,file="InN.RData")


#第二次整合数据
scRNA = scRNA[,scRNA$seurat_clusters %in% c("0","1","2","3","4","5","6","7","8","9")]
cellinfo <- subset(scRNA@meta.data, select= c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","percent.rb","percent.HB","PMI","SampleID","Diagnosis","copykat.pred","Age","Sex","Tangle.Stage","APOE","Region","DoubletFinder","Group","datasets","Age_Group"))
scRNA <- CreateSeuratObject(scRNA@assays$RNA@counts, meta.data = cellinfo)
###--数据标准化
scRNA <- SCTransform(scRNA)
###--PCA
scRNA <- RunPCA(scRNA, npcs=30, verbose=FALSE)
ElbowPlot(scRNA, ndims = 30)
###--数据整合-group.by.vars参数是设置按哪个分组整合
scRNA <- RunHarmony(scRNA, group.by.vars="SampleID",assay.use="SCT",max.iter.harmony=30)
scRNA <- RunTSNE(scRNA, reduction="harmony",dims=1:30) %>% RunUMAP(reduction="harmony", dims = 1:30)
DefaultAssay(scRNA) <- "SCT"
scRNA <- FindNeighbors(scRNA, reduction = "harmony",dims = 1:30) %>% FindClusters(dims = 1:30, resolution = 0.2)
save(scRNA,file="Astro2.RData")

##
DimPlot(scRNA,group.by = "SampleID",reduction = "umap")
DimPlot(scRNA,group.by = "seurat_clusters",reduction = "umap",label = TRUE)
DimPlot(scRNA,group.by = "SampleID", reduction = "umap", split.by = "Region")
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["Sex"]])
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["Group"]])
table(scRNA@meta.data[["Age"]],scRNA@meta.data[["Group"]])
table(scRNA@meta.data[["seurat_clusters"]],scRNA@meta.data[["Region"]])

##
DefaultAssay(scRNA) <- "RNA"
pdf("celltype_marker1.pdf", width = 6.8,height = 3.3)
markers.to.plot <- c("GAD1","GAD2","NR2F2","LHX6","CNR1","CCK","CALB2","VIP","SV2C","LAMP5","RELN","KIT","SLIT2","TLE4","SST","NPY","PVALB","WBSCR17")
DotPlot(scRNA, features = markers.to.plot,group.by = "celltype") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))
dev.off()

##
###为分群重新指定细胞类型 
cluster_celltype <-  data.frame(readxl::read_xlsx("cluster_celltype.xlsx"))
current.cluster.ids <- cluster_celltype$cluster
new.cluster.ids <- cluster_celltype$celltype
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$seurat_clusters,scRNA$celltype)
scRNA = scRNA[,scRNA$celltype %in% c("InN_CNR1_CCK","InN_HIP","InN_LAMP5_KIT","InN_LAMP5_RELN","InN_LAMP5_SV2C","InN_LHX6_SLIT2","InN_LHX6_SST",
                                     "InN_LHX6_TLE4","InN_NR2F2_CNR1","InN_SST_NPY","InN_SV2C_KIT","InN_SV2C_TLE4","InN_VIP_CALB2")]
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
scRNA$celltype <- factor(scRNA$celltype,levels = c("InN_CNR1_CCK","InN_NR2F2_CNR1","InN_VIP_CALB2","InN_SV2C_KIT","InN_LAMP5_RELN",
                                                   "InN_LAMP5_KIT","InN_LAMP5_SV2C","InN_SV2C_TLE4","InN_LHX6_SLIT2","InN_LHX6_TLE4","InN_LHX6_SST",
                                                   "InN_SST_NPY","InN_HIP"))
table(scRNA$celltype)
Idents(scRNA) <- 'celltype'
save(scRNA,file="InN_anno.RData")


##
table(scRNA@meta.data[["seurat_clusters"]], scRNA@meta.data[["datasets"]])
table(scRNA$celltype, scRNA$Region)
prop.table(table(scRNA$celltype, scRNA$Region), margin = 2)*100
scRNA$Region_Group <- paste0(scRNA$Region,"_",scRNA$Group)
prop.table(table(scRNA$celltype, scRNA$Region_Group), margin = 1)*100
prop.table(table(scRNA$copykat.pred, scRNA$Sex), margin = 2)


#鉴定细胞类型之后重新寻找差异基因
#Idents(scRNA)="celltype.datasets"
DefaultAssay(scRNA) <- "RNA"
markers <- FindAllMarkers(scRNA, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T, test.use = "wilcox")  ##耗时久
write_xlsx(markers, "markers.xlsx")
markers_df = markers %>% group_by(cluster) %>% top_n(n = 500, wt = avg_log2FC)
write_xlsx(markers_df, "markers500.xlsx")
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write_xlsx(top10, "markers10.xlsx")

##
###-------------------------------------------Brain--------------------------------------------###
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372',
                '#66934c','#65a252','#da8883','#ca5652','#a28cb1','#d06b4f','#e7b375',
                '#765c91','#da8240','#c4a09f','#9cc28d','#765c91','#bcaf93','#f6f09b',
                '#c39e64',"11111")
my36colors <- c('#467ea7', '#70b5c8','#b76374','#cc929e','#3d6775','#729aa0','#6b2f41',
                '#9c727f','#4a8379','#6daea2','#de9588','#ecc2b8','#45537d','#77819e','#a9adaf','#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#91bf74','#a89372')
my36colors <- c('#985a38', '#aac6d2','#be3935','#649394','#f6f09b','#da8883','#ca5652','#765c91','#da8240','#c4a09f','#9cc28d','#e7b375','#467ea7')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 21, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Group", pt.size=0.01, label = F,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_Group.pdf", plot = p1, device = 'pdf', width = 20, height = 16, units = 'cm')
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Region", pt.size=0.01, label = F,repel = TRUE,cols = my36colors)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "umap_Region.pdf", plot = p1, device = 'pdf', width = 18, height = 16, units = 'cm')



#标准流程
DefaultAssay(scRNA) <- "integrated"
table(aging.list[["Control"]]@meta.data[["celltype"]])#查看去除批次效应之后的结果
scRNA$datasets=str_replace(scRNA$orig.ident,"_.*$","")
scRNA$Sex_Age <- paste(scRNA$Sex, scRNA$Age_stage, sep = "_")

#可视化
Idents(scRNA)="seurat_clusters"
{#tsne
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sex", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Sex", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_SampleID_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "SampleID", groupFacet = 'Sex', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_SampleID_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Group", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Group_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "Group", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_Group_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "SampleID", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_Group_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_1.pdf", plot = p1, device = 'pdf', width = 30, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_cluster_2.pdf", plot = p1, device = 'pdf', width = 30, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "seurat_clusters", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_cluster_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  #umap
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sex", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Sex", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_SampleID_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "SampleID", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_SampleID_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Group", pt.size=0.01)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Group_1.pdf", plot = p1, device = 'pdf', width = 10.5, height = 9, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "Group", pt.size=0.01)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_Group_2.pdf", plot = p1, device = 'pdf', width = 11, height = 9, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "SampleID", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = F, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_Group_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_1.pdf", plot = p1, device = 'pdf', width = 28, height = 25, units = 'cm')
  ggsave(filename = "umap_cluster_1.png", plot = p1, device = 'png', width = 28, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", pt.size=0.01, label = TRUE,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_cluster_2.pdf", plot = p1, device = 'pdf', width = 28, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "seurat_clusters", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = T, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_cluster_3.pdf", plot = p1, device = 'pdf', width = 30, height = 9, units = 'cm')
  
  rm('p1')
}
#
Idents(scRNA)="celltype"
{#tsne
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype", groupFacet = 'copykat.pred', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype_copykat_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "copykat.pred", groupFacet = 'Group', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_copykat_group_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype", groupFacet = 'Sex', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype_Sex_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = T,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_1.pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype_2.pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype_3.pdf", plot = p1, device = 'pdf', width = 30, height = 18, units = 'cm')
  
  #umap
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype", groupFacet = 'copykat.pred', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype_copykat_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "copykat.pred", groupFacet = 'Group', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_copykat_group_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype", groupFacet = 'Sex', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype_Sex_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_1.pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype_2.pdf", plot = p1, device = 'pdf', width = 25, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype_3.pdf", plot = p1, device = 'pdf', width = 30, height = 18, units = 'cm')
  
  rm('p1')
}
#
Idents(scRNA)="celltype1"
{#tsne
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype1", groupFacet = 'copykat.pred', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype1_copykat_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "copykat.pred", groupFacet = 'Group', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_copykat_group_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype1", groupFacet = 'Sex', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype1_Sex_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype1", pt.size=0.01, label = T,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype1_1.pdf", plot = p1, device = 'pdf', width = 34, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype1", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "tsne_celltype1_2.pdf", plot = p1, device = 'pdf', width = 34, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'tsne',clusterCol = "celltype1", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "tsne_celltype1_3.pdf", plot = p1, device = 'pdf', width = 30, height = 18, units = 'cm')
  
  #umap
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype1", groupFacet = 'copykat.pred', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype1_copykat_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "copykat.pred", groupFacet = 'Group', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_copykat_group_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype1", groupFacet = 'Sex', noSplit = F, nrow = 1,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype1_Sex_1.pdf", plot = p1, device = 'pdf', width = 40, height = 20, units = 'cm')
  
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype1", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
    axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype1_1.pdf", plot = p1, device = 'pdf', width = 33, height = 25, units = 'cm')
  p1 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype1", pt.size=0.01, label = F,repel = TRUE)+
    theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
  ggsave(filename = "umap_celltype1_2.pdf", plot = p1, device = 'pdf', width = 33, height = 25, units = 'cm')
  p1 <- clusterCornerAxes(scRNA,reduction = 'umap',clusterCol = "celltype1", groupFacet = 'Group', noSplit = F,pSize=0.01,cellLabel = T, cellLabelSize = 3,show.legend = F, aspect.ratio = 1,themebg = 'bwCorner')
  ggsave(filename = "umap_celltype1_3.pdf", plot = p1, device = 'pdf', width = 30, height = 18, units = 'cm')
  
  rm('p1')
}

####
p11 <- DimPlot(scRNA, reduction = "umap",group.by = "ident", 
               cols =c('#476D87', '#53A85F', '#F1BB72',"grey","grey","grey","grey","grey","#F3B1A0","grey","grey","grey",
                       "grey","grey","grey","grey","grey","grey","grey","grey","grey","#E95C59","grey","grey",
                       "grey","grey","grey","grey","grey","grey","grey"), pt.size=0.5, split.by = "datasets")+
  theme(plot.title=element_text(size=0),plot.tag =element_text(size=0), strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "Oligos.pdf", plot = p11, device = 'pdf', width = 35, height = 10, units = 'cm')

my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

# 检查 PCA 分群结果， 这里只展示前 12 个 PC,每个 PC 只显示 3 个基因； 
print(scRNA[["pca"]], dims = 1:12, nfeatures = 3) 
#绘制 pca 散点图； 去除图例
pca_plot<-DimPlot(scRNA, reduction = "pca")
ggsave(filename = "pca.png", plot = pca_plot, device = 'png', width = 16, height = 12, units = 'cm')
#画前 6 个主成分的热图； 
p11<-DimHeatmap(scRNA, dims = 1:12, cells = 500, balanced = TRUE) 
ggsave(filename = "PCA_Heatmap.png", plot = DimHeatmap(scRNA, dims = 1:12, cells = 500, balanced = TRUE) , device = 'png', width = 15, height = 20, units = 'cm')

#这一步的目的是鉴定出细胞与细胞之间表达量相差很大的基因，用于后续鉴定细胞类型，
#我们使用默认参数，即“vst”方法选取2000个高变基因。
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# 提取表达量变变化最高的 10 个基因； 
top10 <- head(VariableFeatures(scRNA), 10)
top10
# 绘制带有和不带有标签的变量特征的散点图
p11 <- VariableFeaturePlot(scRNA)+NoLegend()
p12 <- LabelPoints(plot = p11, points = top10, repel = TRUE, xnudge=0, ynudge=0)
p11+p12
ggsave(filename = "variablegene.png", plot = p11+p12, device = 'png', width = 20, height = 12, units = 'cm')
rm('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','pca_plot','top10')
rm('tsne1','tsne2','umap2','umap1','tsne_umap1','tsne_umap2','tsne_umap3','tsne_umap4')
#而对所有基因进行标准化的方法如下： 
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes, vars.to.regress = "percent.mt") ##耗时
#线性降维（PCA）,默认用高变基因集，但也可通过 features 参数自己指定； 
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA)) 

##堆叠小提琴图显示差异基因鉴定细胞marker gene
DefaultAssay(scRNA) <- "RNA"
Astro_marker<-VlnPlot(scRNA, features = c("AQP4","GFAP","GPR98","MASS1","SOX9","SLC1A2",	"SLC1A3",	"GLUL",	"VIM","GPC5","RYR3"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "Astro_marker.pdf", plot = Astro_marker, device = 'pdf', width = 10, height = 12, units = 'cm')
rm('Astro_marker')                    

RGL_marker<-VlnPlot(scRNA, features = c("SLC1A3","SOX2",	"AQP4","GFAP","ASCL1","VIM","PAX6","HOPX","NES","PDGFRB","PROX1"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "RGL_marker.pdf", plot = RGL_marker, device = 'pdf', width = 8, height = 12, units = 'cm')
rm('RGL_marker')  

IPC_marker<-VlnPlot(scRNA, features = c("SLC1A3","SOX2",	"AQP4","GFAP","SOX6","ASCL1","RFC4","MAP2","FOXG1"), stacked=T,pt.size=0,combine = FALSE)
ggsave(filename = "IPC_marker.pdf", plot = IPC_marker, device = 'pdf', width = 8, height = 8, units = 'cm')
rm('IPC_marker')  

NSC_marker<-VlnPlot(scRNA,features = c("SLC1A3","SOX2","GFAP","SOX6","SOX5","PAX6","HOPX","NEUROD1","SOX4","SOX11"), stacked=T, pt.size=0,combine = FALSE)
ggsave(filename = "NSC_marker_3.pdf", plot = NSC_marker, device = 'pdf', width = 10, height = 9, units = 'cm')
rm('NSC_marker') 

all_marker<-VlnPlot(scRNA, features = c("CLDN5","RELN","ADARB2","SRGN","ROBO2","CNTNAP2","GRIK1"), stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=4))
ggsave(filename = "all_marker.pdf", plot = all_marker, device = 'pdf', width = 10, height = 10, units = 'cm')
rm('all_marker')

VlnPlot(scRNA, features = markers.to.plot, stacked=T,pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=4))

VlnPlot(scRNA_M, features = c("CLDN5","RELN","ADARB2","SRGN","ROBO2","CNTNAP2","GRIK1"), split.by = "Age_stage",pt.size=0,combine = FALSE)+
  theme(axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=4))

#单个基因差异基因小提琴图
p1<-VlnPlot(scRNA, features = c("MT3"),idents = 'RGL',pt.size = 1,split.by = "datasets",ncol = 1)+theme(axis.title.x=element_text(size=0))
ggsave(filename = "DEG.pdf", plot = p1, device = 'pdf', width = 10, height = 10)


######
##
modify_vlnplot <- function(scRNA, feature, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),...) {
  p <- VlnPlot(scRNA, features = feature, pt.size = pt.size, ... ) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

## main function
StackedVlnPlot <- function(scRNA, features, pt.size = 0, plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) modify_vlnplot(scRNA = scRNA,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45), axis.ticks.x = element_line())
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

#配色方案
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')

cDAM_1 <- StackedVlnPlot(scRNA, c("C3","CSF1R","TYROBP","CTSB"), pt.size=0, split.by = 'datasets',cols=my36colors)
cDAM_2 <- StackedVlnPlot(scRNA, c("CTSD","APOE","B2M","FTH1"), pt.size=0, split.by = 'datasets',cols=my36colors)
fig_marker <- plot_grid(cDAM_1,cDAM_2,ncol = 2)
ggsave(filename = "DAM1.png", plot = cDAM_1, device = 'png', width = 30, height = 20, units = 'cm')
ggsave(filename = "DAM2.png", plot = cDAM_2, device = 'png', width = 30, height = 20, units = 'cm')

#基于细胞marker的cluster鉴定
celltype_marker_neuron=c("PDE1C","ATP13A4", "ZNF536",  "CALD1") #neuron
neuron<-VlnPlot(scRNA,features = celltype_marker_neuron,pt.size = 0,ncol = 1)
celltype_marker_neuron=c("MAN2A1","FTH1", "NLGN4X","ARPP21","SPOCK3" ) #neuron
neuron<-VlnPlot(scRNA,features = celltype_marker_neuron,pt.size = 0,ncol = 1)

#绘制 Marker 基因的 tsne 图
NPCs <- FeaturePlot(scRNA, features = c("SIX1","CITED1","TMEM100","ITGA8","PAX8","HIST1H4C","LYPD1","DAPL1"),ncol = 3) 
ggsave(filename = "NPCs.pdf", plot = NPCs, device = 'pdf', width = 27, height = 20, units = 'cm')
rm(NPCs)


###为分群重新指定细胞类型 
current.cluster.ids <- c("GSM4120425","GSM4120426","GSM4120427","GSM4120428","GSM4432645","GSM4432646","GSM5687879","GSM5687887","GSM5687889")
new.cluster.ids <- c("Aging","Aging","Aging","Aging","Adult","Adult","Aging","Adult","Adult") 
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
scRNA@meta.data$Group <- plyr::mapvalues(x = as.character(scRNA$SampleID), from = current.cluster.ids, to = new.cluster.ids)

###为分群重新指定细胞类型 
cluster_celltype <-  data.frame(readxl::read_xlsx("cluster_celltype.xlsx"))
#scRNA = scRNA[,scRNA@meta.data[["seurat_clusters"]] %in% c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")]
current.cluster.ids <- cluster_celltype$cluster
new.cluster.ids <- cluster_celltype$celltype
scRNA@meta.data$celltype <- plyr::mapvalues(x = as.integer(as.character(scRNA@meta.data$seurat_clusters)), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA$seurat_clusters,scRNA$celltype)
scRNA@meta.data[["celltype"]]<-factor(scRNA@meta.data[["celltype"]], levels=c("ExN","InN","Astro","Astro_Oligo", "Astro","Astro_Astro","Astro_Oligo",
                                                                              "OPC","OPC_Astro","NFOL","MOL","Endo","Peri"))

##
###为分群重新指定细胞类型 
current.cluster.ids <- c("Br5161_Astro", "Br5212_Astro", "Br5400_Astro", "SRR14808188", "SRR14808189", "SRR14808211", "SRR14808212", "SRR14808213", "SRR14808214", "SRR14808215","SRR9262919","SRR9262941")
new.cluster.ids <- c("M","M","F","M","M","M","M","M","F","M","M","M") 
scRNA$Sex <- plyr::mapvalues(x = as.character(scRNA$SampleID), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$SampleID,scRNA$Sex)
current.cluster.ids <- c("Br5161_Astro", "Br5212_Astro", "Br5400_Astro", "SRR14808188", "SRR14808189", "SRR14808211", "SRR14808212", "SRR14808213", "SRR14808214", "SRR14808215","SRR9262919","SRR9262941")
new.cluster.ids <- c("54","52","53","68","68","58","65","90","92","85","19","22") 
scRNA$Age <- plyr::mapvalues(x = as.character(scRNA$SampleID), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$SampleID,scRNA$Age)
current.cluster.ids <- c("Br5161_Astro", "Br5212_Astro", "Br5400_Astro", "SRR14808188", "SRR14808189", "SRR14808211", "SRR14808212", "SRR14808213", "SRR14808214", "SRR14808215","SRR9262919","SRR9262941")
new.cluster.ids <- c("Adult","Adult","Adult","Aging","Aging","Adult","Aging","Longevity","Longevity","Aging","Adult","Adult") 
scRNA$Group <- plyr::mapvalues(x = as.character(scRNA$SampleID), from = current.cluster.ids, to = new.cluster.ids)
table(scRNA@meta.data$SampleID,scRNA$Group)

#绘制 tsne 图(修改标签后的)； 
#tsne
p1 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = F,repel = F,max.overlaps=10)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p2 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE,max.overlaps=10)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
tsne1 <- plot_grid(p1, p2,align = "v",ncol = 2)
ggsave(filename = "tsne1.pdf", plot = tsne1, device = 'pdf', width = 23.5, height = 8, units = 'cm')

p3 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p4 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
tsne2 <- plot_grid(p3, p4,align = "v",ncol = 2)
ggsave(filename = "tsne2-3.pdf", plot = tsne2, device = 'pdf', width = 23.5, height = 8, units = 'cm')

#umap
p5 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = F,repel = F)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p6 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
umap1 <- plot_grid(p5, p6,align = "v",ncol = 2)
ggsave(filename = "umap1.pdf", plot = umap1, device = 'pdf', width = 23.5, height = 8, units = 'cm')

p7 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.1, label =F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
p8 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", pt.size=0.1, label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
umap2 <- plot_grid(p7, p8,align = "v",ncol = 2)
ggsave(filename = "umap2.pdf", plot = umap2, device = 'pdf', width = 23.5, height = 8, units = 'cm')

#tsne_umap
tsne_umap1 <- plot_grid(p1, p5,align = "v",ncol = 2)
ggsave(filename = "tsne_umap1.pdf", plot = tsne_umap1, device = 'pdf', width = 23.5, height = 8, units = 'cm')
tsne_umap2 <- plot_grid(p2, p6,align = "v",ncol = 2)
ggsave(filename = "tsne_umap2.pdf", plot = tsne_umap2, device = 'pdf', width = 23.5, height = 8, units = 'cm')
tsne_umap3 <- plot_grid(p3, p7,align = "v",ncol = 2)
ggsave(filename = "tsne_umap3.pdf", plot = tsne_umap3, device = 'pdf', width = 23.5, height = 8, units = 'cm')
tsne_umap4 <- plot_grid(p4, p8,align = "v",ncol = 2)
ggsave(filename = "tsne_umap4.pdf", plot = tsne_umap4, device = 'pdf', width = 23.5, height = 8, units = 'cm')

#tsne.by datasets
p9 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, split.by = "orig.ident", label = TRUE,repel = TRUE)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "tsne3.pdf", plot = p9, device = 'pdf', width = 25, height = 9.5, units = 'cm')
p9 <- DimPlot(scRNA, reduction = "tsne", group.by = "celltype", pt.size=0.01, split.by = "orig.ident", label = F,repel = F)+
  theme(plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20))
ggsave(filename = "tsne4.pdf", plot = p9, device = 'pdf', width = 25, height = 9.5, units = 'cm')
#umap.by datasets
p10 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", split.by = "Age_stage", pt.size=0.01, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
ggsave(filename = "umap3.pdf", plot = p10, device = 'pdf', width = 25, height =8.5, units = 'cm')
p10 <- DimPlot(scRNA, reduction = "umap", group.by = "celltype", split.by = "Age_stage",pt.size=0.01, label = F,repel = F)+theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),plot.title=element_text(size=0),strip.text=element_text(size=10),axis.title=element_text(size=10))
ggsave(filename = "umap4.pdf", plot = p10, device = 'pdf', width = 25, height =8.5, units = 'cm')

rm('p1','p2','p3','p4','p5','p6','p7','p8','p9','p10')
rm('tsne1','tsne2','umap2','umap1','tsne_umap1','tsne_umap2','tsne_umap3','tsne_umap4')

scRNA$datasets=str_replace(scRNA$orig.ident,"_.*$","")
scRNA$datasets <- sub('\\d+', '', scRNA$datasets)
table(scRNA@meta.data[["datasets"]])
#####
DefaultAssay(scRNA) <- "RNA"
cell_marker<-VlnPlot(scRNA,direction = c("horizontal"),
                     features = c("GFAP","AQP4","SOX2","SLC1A3","NETO1","PROX1","BCL11B","DCX","SEMA3C","MPPED1","CPNE4","SATB2","PCP4","RGS4",
                                  "STMN2","CAMK2A","SNAP25","RBFOX3","GAD1","MOG","PDGFRA","OLIG2","CX3CR1","PTPRC","FLI1","FN1","CCK","LAMP5","SV2C","CNR1",
                                  "SST","CALB2","PVALB","RELN","DNAH9","CFAP54","DCN","COL1A2"),group.by = "seurat_clusters",stacked=T,pt.size=0,combine = FALSE)+
  theme(strip.text=element_text(size=6),
        axis.title=element_text(size=6,vjust = 1),
        axis.ticks=element_line(size=0),
        axis.text.y=element_text(size=6,colour ='black'),
        axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=4),
        axis.title.y=element_text(size=0),
        panel.spacing = unit(-0.1, "lines"))
ggsave(filename = "cluster_marker_3.pdf", plot = cell_marker, device = 'pdf', width = 20, height = 15, units = 'cm')
rm('cell_marker') 

cell_marker<-VlnPlot(scRNA,
                     features = c("MOBP","PLP1","MBP","CNTNAP2","SYT1","RBFOX1","SNAP25","SLC17A7","SATB2","CAMK2A","STMN2","GAD1","GAD2","FLT1","CLDN5","VCAN","SOX6","PDGFRA","GFAP","AQP4","SOX2","SLC1A3","GPC5",
                                  "SLC1A2","CD74","CX3CR1","PTPRC","LPAR6","DOCK8"),group.by = "seurat_clusters",stacked=T,pt.size=0,combine = FALSE)+
  theme(strip.text=element_text(size=6),
        axis.title=element_text(size=6,vjust = 1),
        axis.ticks=element_line(size=0),
        axis.text.y=element_text(size=0,colour ='black'),
        axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=6),
        axis.title.y=element_text(size=0),
        panel.spacing = unit(-0.1, "lines"))
ggsave(filename = "cluster_marker_4.pdf", plot = cell_marker, device = 'pdf', width = 15, height = 20,units = 'cm')
rm('cell_marker') 

cell_marker<-VlnPlot(scRNA,
                     features = c("SLC17A7","GAD1","MOBP","PLP1","MBP","CNTNAP2","SYT1","RBFOX1","FLT1","CLDN5","PDGFRA","GFAP","AQP4","GPC5","SLC1A2","CD74","LPAR6","DOCK8"),
                     group.by = "seurat_clusters",stacked=T,pt.size=0,combine = FALSE)+
  theme(strip.text=element_text(size=6),
        axis.title=element_text(size=6,vjust = 1),
        axis.ticks=element_line(size=0),
        axis.text.y=element_text(size=0,colour ='black'),
        axis.text.x=element_text(vjust = 0, hjust = 0.5,angle=0,size=6),
        axis.title.y=element_text(size=0),
        panel.spacing = unit(-0.1, "lines"))
ggsave(filename = "cluster_marker_5.pdf", plot = cell_marker, device = 'pdf', width = 6, height = 6,units = 'cm')
rm('cell_marker') 


###Dotplot
DefaultAssay(scRNA) <- "RNA"
table(scRNA@meta.data$celltype)

markers.to.plot <- c("MOBP","PLP1","MBP","OPALIN","CNTNAP2","SYT1",
                     "RBFOX1","SNAP25","SLC17A7","MAP1B","SATB2","CAMK2A","STMN2","GAD1","GAD2","DCX",
                     "FLT1","CLDN5","DCN","COL1A2","VCAN","SOX6","PDGFRA","GFAP","AQP4","SOX2","SLC1A3","GPC5","SLC1A2","ATP1B2",
                     "CSF1R","C3","P2RY12","CD74","CX3CR1","PTPRC","LPAR6","DOCK8","CCK","LAMP5","SV2C","CNR1","SST")
markers.to.plot <- c("RBFOX1","SNAP25","SLC17A7","CAMK2A","GAD1","GAD2","DCX","RELN",
                     "SLC1A3","SLC1A2","GFAP","AQP4","CSF1R","CD74","VCAN","SOX6","OLIG1","OLIG2","BCAS1","PLP1","MBP","FLT1","CLDN5","DCN","COL1A2","PTPRC","CD247","CLIC6","DNAH9","CFAP54")

pdf("cluster_marker_1.pdf", width = 16,height = 8)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("cluster_marker_2.pdf", width = 16,height = 8)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
pdf("cluster_marker_3.pdf", width = 9,height = 8)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "seurat_clusters", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("cluster_marker_4.pdf", width = 9,height = 8)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

markers.to.plot <- c("RBFOX1","SNAP25","SLC17A7","CAMK2A","GAD1","GAD2",
                     "SLC1A3","SLC1A2","GFAP","AQP4","CSF1R","CD74","VCAN","SOX6","BCAS1","PLP1","MBP","FLT1","CLDN5","DCN","COL1A2")
pdf("celltype_marker_1.pdf", width = 8,height = 3.5)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("celltype_marker_2.pdf", width = 8,height = 3.5)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

markers.to.plot <- c("RBFOX1","SNAP25","SLC17A7","CAMK2A","CUX2","RORB","FOXP2","PCP4","HS3ST4","THEMIS","TLE4",
                     "GAD1","GAD2","NR2F1","LHX6","LAMP5","KIT","NOS1","PVALB","MEPE","PLEKHH2","SST","VIP","CNR1",
                     "GFAP","AQP4","CSF1R","CD74","VCAN","SOX6","PLP1","MBP","FLT1","CLDN5","DCN","COL1A2","PTPRC","CD247")
pdf("celltype1_marker_1.pdf", width = 15.5,height = 7)
DotPlot(scRNA, features = markers.to.plot, cols = c("blue", "red","green"), group.by = "celltype1", dot.scale = 8) +  RotatedAxis()
dev.off()
pdf("celltype1_marker_2.pdf", width = 15.5,height = 7)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "celltype1", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()

markers.to.plot <- c("COL1A1","COL1A2","ACTA2","TAGLN","ABCC9","P2RY14","CD247","PTPRC","FYN","BCAS1","C1QC","C1QB","C1QA","MEGF10","MERTK","CLIC6","DNAH9","CFAP54")

#CA
markers.to.plot <- c("SNAP25","SATB2","MPPED1","PCP4","RGS4","CPNE4","NETO1")
pdf("CA_marker_1.pdf", width = 5,height = 14)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
#ImN
markers.to.plot <- c("SNAP25","MEIS2","UNC5D","FOXG1","EIF1B","MEF2C","STMN2","NNAT","BCL11B","PROX1","DCX","SEMA3C","NETO1")
pdf("ImN_marker_1.pdf", width = 6,height = 14)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
#L2_3
markers.to.plot <- c("SNAP25","CUX2","RASGRF2","CAMK2A","DCX")
pdf("L2_3_marker_1.pdf", width = 5,height = 14)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
#L4
markers.to.plot <- c("SNAP25","SYT1","RORB","SATB2","FOXP2","FOXP1","DCX")
pdf("L4_marker_1.pdf", width = 5,height = 14)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
#L5_6_THEMIS
markers.to.plot <- c("SNAP25","PCP4","THEMIS","SATB2","FOXP1","CNR1","NTNG2","SULF1","RGS12","TLL1","ADRA1A","TLE4","SOX5","FOXP2","HS3ST4","TUBA1A","NFIB","BCL11B","DCX")
pdf("L5_6_marker_1.pdf", width = 10,height = 14)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()
#InN
markers.to.plot <- c("LHX6","NR2F2","SOX6","DLX1","DLX2","GAD2","GAD1","DCX","PROX1","PVALB","MEPE","PLEKHH2","PLCL1","NPY","SLIT2","SST","ERBB4","RELN","SV2C","KIT","CNR1","LAMP5","CCK","NOS1","CA1","VIP","NOX4","SCTR","ABI3BP","SCML4","CHRNA2","PENK","CALB2")
pdf("InN_marker_1.pdf", width = 11,height = 14)
DotPlot(scRNA, features = markers.to.plot, cols = c("lightgrey", "red"),group.by = "seurat_clusters", col.min = 0,col.max = 2.5, dot.scale = 8) + RotatedAxis()
dev.off()


#绘制marker基因的 umap 图
DefaultAssay(scRNA) <- "integrated"
DefaultAssay(scRNA) <- "RNA"

celltype_umap<- FeaturePlot(scRNA, features = c("AQP4","GFAP"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "Astro_umap.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("MFGE8","ASCL1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_1.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("VIM","PDGFRB"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_2.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)
celltype_umap<- FeaturePlot(scRNA, features = c("NES","PROX1"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 
ggsave(filename = "NSC_umap_3.pdf", plot = celltype_umap, device = 'pdf', width = 9.5, height = 16, units = 'cm')
rm(celltype_umap)

FeaturePlot(scRNA, features = c("HOPX","SLC1A3"),ncol = 1,pt.size = 0.01,cols =c("lightgrey", "#CD0000")) 

DefaultAssay(scRNA) <- "RNA"
p <- Nebulosa::plot_density(scRNA, features = c("SLC1A3","HOPX"),joint = T,reduction = "umap",size = 0.5)[[3]]
ggsave(filename = "RGL_1.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("GFAP","SOX2"),joint = T,reduction = "umap",size = 0.5)[[3]]
ggsave(filename = "RGL_2.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SLC1A3","ASCL1"),joint = T,reduction = "umap",size = 0.5)[[3]]
ggsave(filename = "IPC.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("GFAP","AQP4"),joint = T,reduction = "umap",size = 0.5)[[3]]
ggsave(filename = "Astro.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("PTPRC","CSF1R"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "Astroglia.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("VCAN","PDGFRA"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "OPC.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("MOG","PLP1","BCAS1"),joint = T,reduction = "umap",size = 0.01)[[4]]
ggsave(filename = "NFOL.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("MOG","PLP1"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "MOL.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25"),joint = T,reduction = "umap",size = 0.01)
ggsave(filename = "Neuron.pdf", plot = p, device = 'pdf', width = 10, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25","DCX"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "ImN.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25","CAMK2A"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "ExN.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25","CAMK2A","PROX1"),joint = T,reduction = "umap",size = 0.01)[[4]]
ggsave(filename = "DG_ExN.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25","SLC6A1"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "InN.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25","GAD1","GAD2","NR2F2"),joint = T,reduction = "umap",size = 0.01)[[4]]
ggsave(filename = "CGE_InN.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("SNAP25","GAD1","GAD2","LHX6"),joint = T,reduction = "umap",size = 0.01)[[5]]
ggsave(filename = "MGE_InN.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')
p <- Nebulosa::plot_density(scRNA, features = c("EBF1","RGS5"),joint = T,reduction = "umap",size = 0.01)[[3]]
ggsave(filename = "End.pdf", plot = p, device = 'pdf', width = 11, height = 8, units = 'cm')


##相关性
table(scRNA$seurat_clusters)  
av<-AverageExpression(scRNA,group.by = "seurat_clusters", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_seurat_clusters_1.pdf", width = 7.5,height = 7)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'))
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_seurat_clusters_1.csv") #保存结果
#
table(scRNA$celltype)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype.pdf", width = 3,height = 2.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv") #保存结果
#
table(scRNA$celltype1)  
av<-AverageExpression(scRNA,group.by = "celltype1", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype1.pdf", width = 6,height = 5.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype1.csv") #保存结果
#
table(scRNA$celltype)  
av<-AverageExpression(scRNA,group.by = "celltype", assays = "RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pdf("cor_celltype.pdf", width = 8,height = 7.5)
pheatmap::pheatmap(cor(av[cg,],method = 'spearman'),treeheight_row = 10,treeheight_col=10)
dev.off()
write.csv(cor(av[cg,],method = "spearman"),"cor_celltype.csv") #保存结果


##细胞比例计算做图
Idents(scRNA)="celltype"
table(scRNA$Group)#查看各组细胞数
prop.table(table(Idents(scRNA)))
table(Idents(scRNA), scRNA$Group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(scRNA), scRNA$Group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
#Cellratio$Var2<-factor(Cellratio$Var2, levels=c("BCH","BCL","BCN","Normal"))
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = my36colors)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x=element_text(vjust = 1, hjust = 1,angle=45,size=10))
ggsave(filename = "Cell proportion_2.pdf", device = 'pdf', width = 9, height = 10, units = 'cm')

#
AB<- subset(scRNA,Mki67>1,slot="counts")
##细胞比例计算做图
Idents(scRNA)="seurat_clusters"
Idents(scRNA)="celltype"
prop.table(table(Idents(scRNA)))
table(Idents(scRNA), scRNA$Group)#各组不同细胞群细胞数
Cellratio <- prop.table(table(Idents(scRNA), scRNA$Group), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
library(ggplot2)
ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  scale_fill_manual(values = my36colors)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
ggsave(filename = "Cellratio_all.pdf", device = 'pdf', width = 8, height = 12, units = 'cm')

##保存数据
scRNA = scRNA[,scRNA$celltype %in% c("ExN","InN","Astro","Astro","OPC","MOL","Endo","Peri")]
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
scRNA$celltype <- factor(scRNA$celltype,levels=c("ExN","InN","Astro","Astro","OPC","MOL","Endo","Peri"))

scRNA$celltype1 <- factor(scRNA$celltype1,levels=c("L2_3_CUX2","L4_6_RORB_FOXP2","L5_PCP4","L5_6_HS3ST4_THEMIS","L5_6_HS3ST4_TLE4","L6_THEMIS",
                                                   "InN_LAMP5_KIT","InN_LAMP5_KIT_NOS1","InN_LAMP5_SST","InN_PVALB_MEPE","InN_PVALB_PLEKHH2","InN_SST","InN_VIP","InN_CNR1",
                                                   "Astro","Astro","OPC","MOL","Endo","Peri","T"))
save(scRNA,file="Astro_seurat.RData")

#对所有基因进行scale
DefaultAssay(scRNA)="RNA"
scRNA <- NormalizeData(scRNA)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)

save(scRNA,file="scRNA_Astro_seurat.RData")




