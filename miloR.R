
#devtools::install_github("MarioniLab/miloR", red = "devel")
library(miloR);library(Seurat);library(ggplot2);library(SingleCellExperiment);library(SeuratWrappers);library(scuttle);
library(scater);library(scales);library(forcats);library(data.table);library(stringr);library(dplyr);library(scran);library(patchwork)

unique(scRNA$Group)
#scRNA$Group <- factor(scRNA$Group)
#scRNA<- scRNA[,scRNA$Age %in% c("20W","21W")]
#scRNA<- scRNA[,scRNA$Group %in% c("DS")]
scRNA1 <- as.SingleCellExperiment(scRNA)
reducedDim(scRNA1, "UMAP") <- scRNA@reductions$umap@cell.embeddings
traj_milo <- miloR::Milo(scRNA1)
reducedDim(traj_milo, "UMAP") <- scRNA@reductions$umap@cell.embeddings
reducedDim(traj_milo, "PCA") <- scRNA@reductions$pca@cell.embeddings
reducedDim(traj_milo, "TSNE") <- scRNA@reductions$tsne@cell.embeddings
traj_milo <- miloR::buildGraph(traj_milo, k = 30, d = 50)
traj_milo <- makeNhoods(traj_milo, prop = 0.2, k = 30, d = 50, refined = TRUE)
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), sample = "SampleID")

# Differential abundance testing
traj_design <- data.frame(colData(traj_milo))[, c("SampleID","Group")]
traj_design$SampleID <- as.factor(traj_design$SampleID)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$SampleID
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
traj_design
traj_milo <- calcNhoodDistance(traj_milo, d=30)
rownames(traj_design) <- traj_design$SampleID
da_results <- testNhoods(traj_milo, design= ~Group, design.df = traj_design)
head(da_results)
da_results %>% arrange(- SpatialFDR) %>% head() 
table(da_results$SpatialFDR < 0.1)

#Visualize neighbourhoods displaying DA
traj_milo <- buildNhoodGraph(traj_milo)
plotUMAP(traj_milo, colour_by = "celltype") + plotNhoodGraphDA(traj_milo, da_results, alpha=0.05) + plot_layout(guides="collect")

da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "celltype")
plotDAbeeswarm(da_results, group.by = "celltype") +
  scale_color_gradient2(low="#070091", mid="lightgrey",high="$910000",limits=c(-5,5),oob=squish) +
  lab(x="", y="Log2 Fold Change") +
  theme_bw(base_size = 10) +
  theme(axis.text = element_text(colour = 'black'))


## Plot single-cell UMAP
umap_pl <- plotReducedDim(traj_milo, dimred = "umap", colour_by="stage", text_by = "celltype", text_size = 3, point_size=0.5) + guides(fill="none")
## Plot neighbourhood graph
nh_graph_pl <- plotNhoodGraphDA(traj_milo, da_results, layout="umap",alpha=0.1) 
umap_pl + nh_graph_pl + plot_layout(guides="collect")

###
ggplot(da_results, aes(celltype_fraction)) + geom_histogram(bins=50)
da_results$celltype <- ifelse(da_results$celltype_fraction < 0.7, "Mixed", da_results$celltype)
plotDAbeeswarm(da_results, group.by = "celltype")

#########################################------------multiple groups----------#########################################
#--Define cell neighbourhoods
scRNA1 <- as.SingleCellExperiment(scRNA)
traj_milo <- miloR::Milo(scRNA1)
reducedDim(scRNA1, "UMAP") <- scRNA@reductions$umap@cell.embeddings
reducedDim(traj_milo, "UMAP") <- scRNA@reductions$umap@cell.embeddings
reducedDim(traj_milo, "PCA") <- scRNA@reductions$pca@cell.embeddings
reducedDim(traj_milo, "TSNE") <- scRNA@reductions$tsne@cell.embeddings

#plotUMAP(traj_milo, colour_by="celltype") + plotUMAP(traj_milo, colour_by="Group")
traj_milo <- buildGraph(traj_milo, k = 10, d = 20)
traj_milo <- makeNhoods(traj_milo, prop = 0.9, k = 10, d=20, refined = TRUE, refinement_scheme="graph")
traj_milo <- countCells(traj_milo, meta.data = data.frame(colData(traj_milo)), samples="SampleID")
#--Differential abundance testing with contrasts
traj_design <- data.frame(colData(traj_milo))[,c("SampleID", "Group")]
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$SampleID
## Reorder rownames to match columns of nhoodCounts(milo)
traj_design <- traj_design[colnames(nhoodCounts(traj_milo)), , drop=FALSE]
table(traj_design$Group)
rownames(traj_design) <- traj_design$SampleID
contrast.1 <- c("GroupAging - GroupAdult")
contrast.1 <- c("GroupLongevity - GroupAdult")
contrast.1 <- c("GroupLongevity - GroupAging")
da_results <- testNhoods(traj_milo, design = ~ 0 + Group, design.df = traj_design, model.contrasts = contrast.1, fdr.weighting="graph-overlap")
da_results %>% arrange(- SpatialFDR) %>% head() 
table(da_results$SpatialFDR < 0.1)

#Visualize neighbourhoods displaying DA
traj_milo <- buildNhoodGraph(traj_milo)
#可视化1----------------------------------------------------------------------------
#UMAP图展示neighborhoods
pdf('milo_PMN_Longevity_Aging-Adult.pdf', width=18, height=8)
plotUMAP(traj_milo, colour_by = "celltype",point_size = 0.01) + theme(
  axis.line = element_blank(),axis.ticks = element_blank(),axis.text = element_blank(),
  plot.title=element_text(size=0),strip.text=element_text(size=20),axis.title=element_text(size=20)) +
  plotNhoodGraphDA(traj_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",#修改颜色
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish) + 
  plot_layout(guides="collect")
dev.off()
#可视化2----------------------------------------------------------------------------
#蜂群图展示celltype logFC变化
da_results <- annotateNhoods(traj_milo, da_results, coldata_col = "celltype")
pdf('milo_beeswarm_Aging-Adult.pdf', width=5, height=4)
plotDAbeeswarm(da_results, group.by = "celltype") +
  scale_color_gradient2(low="#070091",
                        mid="lightgrey",
                        high="#910000",
                        #limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black')) 
dev.off()

###----multiple comparisons
contrast.all <- c("GroupAdult - GroupInfancy", "GroupAging - GroupAdult")
# this is the edgeR code called by `testNhoods`
model <- model.matrix(~ 0 + Group, data=traj_design)
mod.constrast <- makeContrasts(contrasts=contrast.all, levels=model)
mod.constrast
##then simple pair-wise comparisons


### group levels analysis
model <- model.matrix(~ 0 + Group, data=traj_design)
ave.contrast <- c("(GroupAging + GroupLongevity)/2  - GroupAdult")
ave.contrast <- c("GroupLongevity - (GroupAging + GroupAdult)/2")
mod.constrast <- makeContrasts(contrasts=ave.contrast, levels=model)
mod.constrast
da_results <- testNhoods(traj_milo, design = ~ 0 + Group, design.df = traj_design, model.contrasts = ave.contrast, fdr.weighting="graph-overlap")

table(da_results$SpatialFDR < 0.1)
traj_milo <- buildNhoodGraph(traj_milo)
plotUMAP(traj_milo, colour_by="SubType") + plotNhoodGraphDA(traj_milo, da_results, alpha=0.1) + plot_layout(guides="auto" )









