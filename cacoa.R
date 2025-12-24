library(cacoa);library(conos);library(Matrix);library(ggplot2);library(ggrastr);library(cowplot);library(Seurat)
library(Cairo);library(sccore);library(tidyverse);library(tidyr);library(data.table);library(ggsignif)

#
table(scRNA$Group)
DefaultAssay(scRNA) <- "integrated"
scRNA <- FindNeighbors(scRNA, graph.name= c("graph1","graph2"),dims = 1:30)
DefaultAssay(scRNA) <- "RNA"
#
scRNA_all <- scRNA
table(scRNA_all$Group)
scRNA <- scRNA_all[,scRNA_all$Group %in% c("Adult","Aging")]
table(scRNA$SampleID,scRNA$Group)

Adult <- readxl::read_xlsx("Adult.xlsx")
Aging <- readxl::read_xlsx("Aging.xlsx")
scRNA$Group <- as.factor(as.character(scRNA$Group))
# sample.groups: vector with condition labels per sample named with sample ids
samplegroups <- list(
  Adult = c("Br5161_HPC","Br5212_HPC","Br5287_HPC","GSM5618238_16","GSM5618239_14","GSM5618240_17","GSM5967900_47","GSM5967901_48"),
  Aging = c("GSE198323_24","GSE198323_26","GSE198323_28","GSE198323_29","GSE198323_30","GSE198323_31","GSE198323_34","GSE198323_38",
            "GSM5618239_21","GSM5618240_18","GSM5618241_19","GSM5967902_49","GSM5967903_50","GSM5967904_51","GSM5967905_52","GSM6280585","GSM6280586"))
diseasef <- as.factor(setNames(rep(names(samplegroups),unlist(lapply(samplegroups,length))),unlist(samplegroups)))
# cell.groups: cell type annotation vector named by cell ids
scRNA$celltype <- as.character(scRNA$celltype)
scRNA$SampleID <- as.character(scRNA$SampleID)
cell.groups <- as.character(scRNA$celltype)
names(cell.groups) <- colnames(scRNA)
# sample.per.cell: vector with sample labels per cell named with cell ids
sample.per.cell <- as.character(scRNA$SampleID)
names(sample.per.cell) <- colnames(scRNA)
# ref.level: id of the condition, corresponding to the reference (i.e. Adult)
ref.level <- "Adult"
# target.level: id of the condition, corresponding to the target (i.e. case)
target.level <- "Aging"
#
embedding <- data.frame(scRNA@reductions$umap@cell.embeddings)
#
cao <- Cacoa$new(scRNA, sample.groups=diseasef, cell.groups=cell.groups, sample.per.cell=sample.per.cell, 
                 ref.level=ref.level, target.level=target.level, graph.name= "graph1", embedding = embedding, n.cores = 20)
cao$plot.params <- list(size=0.1, alpha=0.1, font.size=c(2, 3))

##
N_CORES <- 20;FORCE <- TRUE;verb1 <- TRUE;verb2 <- TRUE
cao$estimateExpressionShiftMagnitudes(n.permutations=5000, verbose=verb2, min.samp.per.type=2)
cao$estimateExpressionShiftMagnitudes(n.permutations=2500, top.n.genes=500, n.pcs=5,min.samp.per.type=2, name='es.top.de', verbose=verb2)
cao$estimateCellLoadings(n.seed = 1, n.boot = 1000, n.cores = 20, filter.empty.cell.types = F)
for (met in c('kde', 'graph')) {
  cn <- paste0('cell.density.', met)
  cao$estimateCellDensity(method=met, estimate.variation=FALSE, verbose=verb2, name=cn, beta=10)
  #cao$estimateDiffCellDensity(type='wilcox', adjust.pvalues=TRUE, verbose=verb2, n.permutations=500, name=cn)
  #cao$estimateDiffCellDensity(type='subtract', adjust.pvalues=FALSE, verbose=verb2, name=cn)
  cao$estimateDiffCellDensity(type = "t.test", verbose=TRUE, n.cores = 10, name=cn, smooth = TRUE)
}
cao$estimateClusterFreeExpressionShifts(n.top.genes=3000, gene.selection="expression", verbose=verb1)
cao$estimateClusterFreeDE(n.top.genes=1500, min.expr.frac=0.01, adjust.pvalues=TRUE, smooth=TRUE, verbose=verb1)
cao$smoothClusterFreeZScores(progress.chunks=10, z.adj=TRUE, verbose=verb1)
cao$estimateDEPerCellType(independent.filtering=TRUE, test='DESeq2.Wald', verbose=verb1)
cao$estimateOntology(type="GSEA", org.db=org.Hs.eg.db::org.Hs.eg.db, verbose=verb1, n.cores=1)
#write_rds(cao)
#message(d.name, " done!")

####----------Compositional differences
##
pdf("CellDensity_1.pdf", width = 8,height = 4)
cao$plotCellDensity(name='cell.density.kde') %>% plot_grid(plotlist=., ncol=2)
dev.off()
#
g0 <- cao$plotEmbedding(color.by='cell.groups')
plot_grid(g0, cao$plotDiffCellDensity(name='cell.density.kde', legend.position=c(0, 1)), ncol=2)
cao$plotCellDensity(type='t.test',name='cell.density.graph') %>% plot_grid(plotlist=., ncol=2)
#
pdf("CellDensity_2.pdf", width = 8,height = 4)
plot_grid(g0, cao$plotDiffCellDensity(name='cell.density.kde',size = 1, legend.position=c(0, 1),adjust.pvalues = F,show.legend = T), ncol=2)
dev.off()
#
pdf("CellDensity_3.pdf", width = 6,height = 5)
cao$plotDiffCellDensity(type='t.test', name='cell.density.kde',size = 1, adjust.pvalues = F,show.legend = T)
dev.off()
#
pdf("CellDensity_4.pdf", width = 4,height = 4)
cao$plotCellLoadings()
dev.off()
  
####----------Expression differences
###
pdf("Exp_distances_1.pdf", width = 6,height = 6)
lapply(c("coda", "expression.shifts"), function(sp) {
  cao$plotSampleDistances(space=sp, legend.position=c(1, 1))
}) %>% 
  plot_grid(plotlist=., ncol=2, labels=c("CoDA", "Expression"), hjust=0, label_x=0.02, label_y=0.99)
dev.off()

#
my36colors <- c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                '#968175')
sample_meta <- cao$data.object@meta.data
lapply(c('coda', 'expression'), function(sp) {
  cao$plotSampleDistances(space=sp, legend.position=c(1, 1), palette=my36colors,sample.colors=sample_meta$SampleID)
}) %>% plot_grid(plotlist=., ncol=2)

#
pdf("Exp_distances_2.pdf", width = 6,height = 6)
cao$plotExpressionDistance(joint = TRUE, notch = F, show.significance = TRUE)+ 
  scale_fill_manual(labels = c("Adult", "Aging"), values = c("#980e5c","#435790")) + 
  theme(legend.text = element_text(size = 15), axis.title.y = element_text(size = 15))
dev.off()

#
pdf("Exp_distances_3.pdf", width = 3,height = 3)
cao$plotExpressionDistance(notch = F, show.significance = TRUE, alpha = 0.1)+ 
  scale_fill_manual(labels = c("Adult", "Aging"), values = c("#980e5c","#435790"))
dev.off()

#
#cao$plotExpressionShiftMagnitudes()
pdf("Exp_distances_4.pdf", width = 2,height = 4)
cao$plotExpressionShiftMagnitudes(notch = T, name='es.top.de')
dev.off()
#
#cao$estimateCommonExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, min.samp.per.type=3, min.gene.frac=0.05)
#caom$estimateCommonExpressionShiftMagnitudes(n.permutations=2500, min.cells.per.sample=10, min.samp.per.type=3, min.gene.frac=0.05)

#
cao$plotVolcano(xlim=c(-3.5, 3.5), ylim=c(0, 3.5), lf.cutoff=1)

## 报错
sample_meta <- cao$data.object@meta.data
lapply(c("coda", "expression.shifts"), function(sp) {
  smd <- as.data.frame(sample_meta) %>% dplyr::select(Group,Age,PMI,Sex)
  sep.res <- cao$estimateMetadataSeparation(smd, space=sp, dist="l1", name=paste0("md.", sp),
                                            show.warning=FALSE)
  (-log10(sep.res$padjust)) %>% {tibble(Type=names(.), value=.)} %>%
    cacoa:::plotMeanMedValuesPerCellType(type="bar", yline=-log10(0.05), ylab="-log10(separation p-value)") +
    scale_y_continuous(expand=c(0, 0, 0.05, 0)) +
    scale_fill_manual(values=rep("#2b8cbe", length(sample_meta))) +
    theme(axis.title.y=element_text(size=13))
}) %>% 
  plot_grid(plotlist=., ncol=2, labels=c("CoDA", "Expression"), hjust=0, label_x=0.22, label_y=0.98)


##
cao$estimateDEPerCellType(
  independent.filtering=TRUE, test='DESeq2.Wald', verbose=FALSE, resampling.method='fix.samples', 
  fix.n.samples=6, n.cells.subsample=30, name='de.fix.samples', n.resamplings=50
)
cao$plotNumberOfDEGenes(
  name="de.fix.samples", type="box", show.resampling.results=TRUE, jitter.alpha=0.5,
  show.jitter=TRUE, y.offset=1
) + scale_y_log10(labels=c(0, 10, 100), breaks=c(1, 11, 101), expand=c(0, 0), limits=c(1, 300))

save(cao,file="Aging_PFC_ExN_cao.RData")



