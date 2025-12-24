#scTenifoldKnk
#scTenifoldKnk is a machine learning workflow performing virtual KO experiments to predict gene function.
#It constructs gene regulatory networks using single-cell RNA sequencing data from wild-type samples and then computationally deletes target genes. 
#The source code of scTenifoldKnk is available at https://github.com/cailab-tamu/scTenifoldKnk
#scTenifoldKnk has been implemented in R, Python, Julia, and Matlab. 
#The R package is available at the CRAN repository at https://cran.r-project.org/web/packages/scTenifoldKnk/. 
#The Matlab application is available in scGEAToolbox
#安装包---install.packages('scTenifoldKnk') install_github('cailab-tamu/scTenifoldKnk')
#https://github.com/cailab-tamu/scTenifoldKnk/tree/master/R
library(scTenifoldKnk);library(Seurat);
#devtools::install_github("MarioniLab/miloR", red = "devel")
library(miloR);library(Seurat);library(ggplot2);library(SingleCellExperiment);library(SeuratWrappers);library(scuttle);
library(scater);library(scales);library(forcats);library(data.table);library(stringr);library(dplyr);library(scran);library(patchwork)

#scRNA=subset(scRNA,downsample=50)
# ExN
DotPlot(scRNA1, features = c("LRRTM4"),group.by = "Group") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))

# Micro
DotPlot(scRNA1, features = c("SPP1"),group.by = "Group") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))

# MOL
DotPlot(scRNA1, features = c("VCAN","CD81","APOD"),group.by = "Group") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))

# OPC
DotPlot(scRNA1, features = c("ANK2"),group.by = "Group") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))

# Astro
DotPlot(scRNA1, features = c("CLU","DPP10","DPYSL2"),group.by = "Group") + 
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(hjust = 1,vjust=1,angle=90))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#3F6699","#FFFFFF","#923331"))

##
###---实际数据
table(scRNA$celltype)
table(scRNA$Group)
scRNA1 <- scRNA[,scRNA$celltype %in% c("OPC")]
scRNA_seq <- scRNA1[,scRNA1$Group %in% c("Aging")]
scRNA_seq$Group <- as.factor(as.character(scRNA_seq$Group))
X <- data.frame(scRNA_seq@assays$RNA@counts)
save(X ,file="HIP_OPC_Aging_KO_counts.RData")

Idents(scRNA_seq) <- "Group"
scRNA_seq <- subset(scRNA_seq,downsample = 10000)






