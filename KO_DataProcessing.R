
library(fgsea);library(ggplot2);library(enrichR);library(igraph)
source('/Users/niuruize/Downloads/scRNA/R_code/plotDR.R')
source('/Users/niuruize/Downloads/scRNA/R_code/plotKO.R')
source('/Users/niuruize/Downloads/scRNA/R_code/hsa2mmu_SYMBOL.R')

load('HIP_ExN_Aging_LRRTM4_counts.RData')
ExN_LRRTM4 = oX
write.csv(ExN_LRRTM4$diffRegulation, 'ExN_LRRTM4_counts.csv')
dGenes <- ExN_LRRTM4$diffRegulation$gene[ExN_LRRTM4$diffRegulation$p.adj < 0.05]

#png('dr2_ExN_LRRTM4.png', width = 2000, height = 2000, res = 300)
pdf("dr2_ExN_LRRTM4.pdf", width = 5,height = 5)
plotDR(ExN_LRRTM4)
dev.off()

#png('ego2_ExN_LRRTM4.png', width = 3000, height = 3000, res = 300, bg = NA)
pdf("ego2_ExN_LRRTM4.pdf", width = 10,height = 10)
X <- ExN_LRRTM4
gKO <- 'LRRTM4'
q <- 0.995
gList <- unique(c(gKO, X$diffRegulation$gene[X$diffRegulation$p.adj < 0.05]))
sCluster <- as.matrix(X$tensorNetworks$WT[gList,gList])
koInfo <- sCluster[gKO,]
sCluster[abs(sCluster) <= quantile(abs(sCluster), q)] <- 0
sCluster[gKO,] <- koInfo
diag(sCluster) <- 0
sCluster <-  reshape2::melt(as.matrix(sCluster))
colnames(sCluster) <- c('from', 'to', 'W')
sCluster <- sCluster[sCluster$W != 0,]
netPlot <- graph_from_data_frame(sCluster, directed = TRUE)
dPlot <- centr_degree(netPlot)$res
W <- rep(1,nrow(sCluster))
sG   <- (names(V(netPlot))[dPlot > 1])[-1]
W[sCluster$from %in% sG] <- 0.2
W[sCluster$to %in% sG] <- 0.2
W[sCluster$from %in% gKO] <- 1
W[sCluster$from %in% gKO & sCluster$to %in% sG] <- 0.8
set.seed(1)
layPlot <- layout_with_fr(netPlot, weights = W)
dPlot <- (dPlot/max(dPlot))*20
dbs <- listEnrichrDbs()
E <- enrichr(gList, c("GO_Biological_Process_2023","KEGG_2021_Human"))
E <- do.call(rbind.data.frame, E)
#E <- E[E$Adjusted.P.value < 0.05,]
E <- E[order(E$Adjusted.P.value),]
E$Term <- unlist(lapply(strsplit(E$Term,''), function(X){
  X[1] <- toupper(X[1])
  X <- paste0(X,collapse = '')
  X <- gsub('\\([[:print:]]+\\)|Homo[[:print:]]+|WP[[:digit:]]+','',X)
  X <- gsub("'s",'',X)
  X <- unlist(strsplit(X,','))[1]
  X <- gsub('[[:blank:]]$','',X)
  return(X)
}))
write.csv(E,"LRRTM4_EGO.csv") #保存结果
E <- E[E$Term %in% c('Nervous System Development','Positive Regulation Of Neuron Differentiation','Regulation Of Synapse Assembly',
                     'Regulation Of Neural Precursor Cell Proliferation'),]
#E <- E[c(1,2,3,4,5,6),]
tPlot <- strsplit(E$Genes, ';')
pPlot <- matrix(0,nrow = length(V(netPlot)), ncol = nrow(E))
rownames(pPlot) <- toupper(names(V(netPlot)))
for(i in seq_along(tPlot)){
  pPlot[unlist(tPlot[i]),i] <- 1
}
pPlot <- lapply(seq_len(nrow(pPlot)), function(X){as.vector(pPlot[X,])})
names(pPlot) <- names(V(netPlot))
tPlot <- unique(unlist(tPlot))
eGenes <- toupper(names(V(netPlot))) %in% tPlot
vColor <- rgb(195/255, 199/255, 198/255 ,0.3)
pieColors <- list(hcl.colors(nrow(E), palette = 'Zissou 1', alpha = 0.7))
par(mar=c(4,0,0,0), xpd = TRUE)
suppressWarnings(plot(netPlot,
                      layout = layPlot,
                      edge.arrow.size=.2,
                      vertex.label.color="black",
                      vertex.shape = ifelse(eGenes,'pie','circle'),
                      vertex.pie = pPlot,
                      vertex.size = 10+dPlot,
                      vertex.pie.color=pieColors,
                      vertex.label.family="Times",
                      vertex.label.font=ifelse(eGenes,2,1),
                      edge.color = ifelse(E(netPlot)$W > 0, 'red', 'blue'),
                      edge.curved = ifelse(W == 0.2, 0, 0.1),
                      vertex.color = vColor,
                      vertex.frame.color = NA))
sigLevel <- formatC(E$P.value, digits = 2, format = 'g', width = 0, drop0trailing = TRUE)
gSetNames <- lengths(strsplit(E$Genes, ';'))
gSetNames <- paste0('(', gSetNames,') ', E$Term, ' FDR = ', sigLevel)
legend(x = -1.05, y = -1.05, legend = gSetNames, bty = 'n', ncol = 2, cex = 1, col = unlist(pieColors), pch = 16)
dev.off()


MGI <- gmtPathways('https://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=KEGG_2019_Human')

zExN_LRRTM4 <- ExN_LRRTM4$diffRegulation$Z
names(zExN_LRRTM4) <- toupper(ExN_LRRTM4$diffRegulation$gene)
set.seed(1)
E <- fgseaMultilevel(MGI, zExN_LRRTM4)
write.csv(data.frame(lapply(E, as.character)),"LRRTM4_GSEA.csv") #保存结果

#png('gsea1_ExN_LRRTM4.png', width = 1000, height = 1000, res = 300)
pdf("gsea1_ExN_LRRTM4.pdf", width = 5,height = 5)
gSet <- 'Axon guidance'
plotEnrichment(MGI[[gSet]], zExN_LRRTM4) +
  labs(
    title = 'Axon guidance',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_ExN_LRRTM4.png', width = 1000, height = 1000, res = 300)
pdf("gsea2_ExN_LRRTM4.pdf", width = 5,height = 5)
gSet <- 'Neuroactive ligand-receptor interaction'
plotEnrichment(MGI[[gSet]], zExN_LRRTM4) +
  labs(
    title = 'Neuroactive ligand-receptor interaction',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_ExN_LRRTM4.png', width = 1000, height = 1000, res = 300)
pdf("gsea3_ExN_LRRTM4.pdf", width = 5,height = 5)
gSet <- 'Long-term potentiation'
plotEnrichment(MGI[[gSet]], zExN_LRRTM4) +
  labs(
    title = 'Long-term potentiation',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_ExN_LRRTM4.png', width = 1000, height = 1000, res = 300)
pdf("gsea4_ExN_LRRTM4.pdf", width = 5,height = 5)
gSet <- 'Synaptic vesicle cycle'
plotEnrichment(MGI[[gSet]], zExN_LRRTM4) +
  labs(
    title = 'Synaptic vesicle cycle',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_ExN_LRRTM4.png', width = 1000, height = 1000, res = 300)
pdf("gsea5_ExN_LRRTM4.pdf", width = 5,height = 5)
gSet <- 'Mismatch repair'
plotEnrichment(MGI[[gSet]], zExN_LRRTM4) +
  labs(
    title = 'Mismatch repair',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

#png('../Results/gsea2_ExN_LRRTM4.png', width = 1000, height = 1000, res = 300)
pdf("gsea6_ExN_LRRTM4.pdf", width = 5,height = 5)
gSet <- 'Glutamatergic synapse'
plotEnrichment(MGI[[gSet]], zExN_LRRTM4) +
  labs(
    title = 'Glutamatergic synapse',
    subtitle = paste0('FDR = ', formatC(E$padj[E$pathway %in% gSet], digits = 2, format = 'e'))) +
  xlab('Gene rank') +
  ylab('Enrichment Score') + theme(plot.title = element_text(face = 2, size = 25))
dev.off()

