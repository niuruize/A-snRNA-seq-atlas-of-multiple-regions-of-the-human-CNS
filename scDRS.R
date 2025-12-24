
# sample cells
scRNA$Group_celltype <- paste(scRNA$Group, scRNA$celltype, sep = "_")
table(scRNA$Group_celltype)
Idents(scRNA) <- "Group_celltype"

library(SeuratDisk)
scRNA$celltype <- as.character(scRNA$celltype)
count <- scRNA@assays$RNA@counts
scRNA$Group_celltype <- paste(scRNA$Group, scRNA$celltype, sep = "_")
scRNA2<-CreateSeuratObject(counts=count, meta.data = scRNA@meta.data)
SaveH5Seurat(scRNA2,filename="HIP.h5Seurat",overwrite = T)
Convert("HIP.h5seurat",dest="h5ad",assay="RNA",overwrite = T)

HIP_cov <- subset(scRNA@meta.data, select= c("Sex","percent.mt","percent.rb","nCount_RNA"))
HIP_cov$Sex = factor(HIP_cov$Sex)
HIP_cov$percent.mt = as.numeric(HIP_cov$percent.mt)
HIP_cov$percent.rb = as.numeric(HIP_cov$percent.rb)
HIP_cov$nCount_RNA = as.numeric(HIP_cov$nCount_RNA)
write.table(HIP_cov, file='HIP_cov.tsv',sep = '\t')

# save metadata table:
scRNA$barcode <- colnames(scRNA)
scRNA$UMAP_1 <- scRNA@reductions$umap@cell.embeddings[,1]
scRNA$UMAP_2 <- scRNA@reductions$umap@cell.embeddings[,2]
write.csv(scRNA@meta.data, file='metadata.csv', quote=F, row.names=F)

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(scRNA@reductions$pca@cell.embeddings, file='pca.csv', quote=F, row.names=F)

# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(scRNA, assay='RNA', slot='counts')
writeMM(counts_matrix, file='counts.mtx')

# write gene names
write.table(data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',quote=F,row.names=F,col.names=F)


###---analysis in python
import os
os.environ["KMP_DUPLICATE_LIB_OK"]="TRUE"
os.chdir("/Users/niuruize/Downloads/scRNAseq/TS_HIP/15_scDRS")
os.getcwd()

import scanpy as sc
import anndata 
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

type="TS"


# load sparse matrix:
X = io.mmread("counts.mtx")

# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr())

# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")

# load gene names:
with open("gene_names.csv", 'r') as f:
  gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['SampleID'], frameon=False, save=True)

# save dataset as anndata format
adata.write('my_data.h5ad')

# reload dataset
adata = sc.read_h5ad('my_data.h5ad')


