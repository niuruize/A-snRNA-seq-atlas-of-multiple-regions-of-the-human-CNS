
##tmp= AB[,AB@meta.data$celltype %in% c("Uro-Epi", "T-cell","tubular", "Mono","LOH","Endo")]
#expr <- tmp@assays[["RNA"]]@counts
  
# Load R libraries ####
library(Seurat)
library(Matrix)
library(pheatmap)
library(sfsmisc)
library(MASS)
library(hopach)
  
# Load Aging Seurat object ####
scRNA = scRNA[,scRNA$celltype %in% c("ExN","InN","DaN","Astro","Micro","OPC","MOL","Endo","Peri")] #Epend数量较少，不纳入分析
scRNA$celltype <- as.factor(as.character(scRNA$celltype))
scRNA$celltype <- factor(scRNA$celltype,levels=c("ExN","InN","DaN","Astro","Micro","OPC","MOL","Endo","Peri"))

AB = subset(scRNA, downsample = 20000)
table(AB@meta.data[["celltype"]])
AB$grouping = AB$Group
celltypes <- unique(AB@meta.data$celltype)
celltypes <- celltypes[which(!is.na(celltypes))]
celltypes <- setdiff(celltypes, c("ExN","InN","Astro","Micro","OPC","Endo","Peri","T"))
getEuclideanDistance <- function(celltype, lowcv = T){
  print(paste("Working on", celltype))
  library(hopach)
  tmp= AB[,AB@meta.data$celltype %in% c("MOL")]
  expr <- tmp@assays[["RNA"]]@counts
  
  zeros <- which(Matrix::rowSums(expr) == 0)
  expr <- data.matrix(expr[-zeros,])
  
  Down_Sample_Matrix <-function (expr_mat) {
    min_lib_size <- min(colSums(expr_mat))
    down_sample <- function(x) {
      prob <- min_lib_size/sum(x)
      return(unlist(lapply(x, function(y) {
        rbinom(1, y, prob)
      })))
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
  }
  ds_expr <- Down_Sample_Matrix(expr)
  
  nsample <- min(table(tmp@meta.data$grouping)[c("Adult", "Aging","Longevity")])
  if(nsample < 10){
    print("Not enough cells")
    return(NULL)
  } 
  old_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "Longevity")], nsample)
  adult_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "Aging")], nsample)
  young_r <- sample(rownames(tmp@meta.data)[which(tmp@meta.data$grouping == "Adult")], nsample)
  ds_expr_r <- ds_expr[, c(young_r, adult_r,old_r)]
  
  if(lowcv){
    getLowCVgenes <- function(matr){
      means <- Matrix::rowMeans(matr)
      bins <- quantile(means, c(seq(from = 0, to = 1, length = 11)))
      mean_bin <- unlist(lapply(means, function(x) min(which(bins >= x))))
      asplit <- split(names(means), mean_bin)
      genes <- unique(unlist(lapply(asplit[setdiff(names(asplit), c("1", "11"))], function(x){
        coef_var <- apply(matr, 1, function(x) sd(x)/mean(x))
        bottom10percent <- names(head(sort(coef_var), round(10*length(coef_var))))
      })))
      genes
    }
    genes <- getLowCVgenes(ds_expr_r)
  }
  else{
    genes <- rownames(ds_expr_r)
  }
  
  calcEuclDist <- function(matr, young, adult, old){
    tmp <- data.matrix(sqrt(matr[genes, young]))
    mean <- rowMeans(sqrt(matr[genes, young]))
    d_young <- distancevector(t(tmp), mean , d="euclid")
    names(d_young) <- young
    tmp <- data.matrix(sqrt(matr[genes, adult]))
    mean <- rowMeans(sqrt(matr[genes, adult]))
    d_adult <- distancevector(t(tmp), mean , d="euclid")
    names(d_adult) <- adult
    tmp <- data.matrix(sqrt(matr[genes, old]))
    mean <- rowMeans(sqrt(matr[genes, old]))
    d_old <- distancevector(t(tmp), mean , d="euclid")
    names(d_old) <- old
    
    list(young = d_young, adult=d_adult, old = d_old)
  }
  ds <- calcEuclDist(matr = ds_expr_r, old = old_r, adult = adult_r, young = young_r)
  ds
}

# Run for all celltypes ####
res <- lapply(celltypes, function(x) getEuclideanDistance(x, lowcv = F))
names(res) <- celltypes
MOL <- res
res<-c(ExN,InN,Astro,Micro,OPC,MOL,Endo,Peri,T)

res_original <- res
# Calculate mean differences and p-values ####
diffs <- unlist(lapply(res_original, function(x) log2(mean(x[[2]]) / mean(x[[1]]))))
pvals <- unlist(lapply(res_original, function(x) wilcox.test(x[[1]], x[[2]])$p.value))
adj_pvals <- p.adjust(pvals, method = "BH")
sizes <- (-log10(adj_pvals))
sizes[which(sizes < 1)] <- 1
sizes[which(sizes > 4)] <- 4
sizes <- sizes * 0.75
farben <- rep("grey", length(adj_pvals))
farben[which(adj_pvals < 0.05)] <- "purple"
diffs_2 <- unlist(lapply(res_original, function(x) log2(mean(x[[3]]) / mean(x[[2]]))))
pvals_2 <- unlist(lapply(res_original, function(x) wilcox.test(x[[2]], x[[3]])$p.value))
adj_pvals_2 <- p.adjust(pvals_2, method = "BH")

save(res,file="res.RData")

pdf("Transcriptional_noise_1.pdf",width = 4,height = 5)
ord <- rev(order(diffs))
par(mar = c(15,5,2,5))
boxplot(do.call(c, res[ord]), las = 2, outline = F, col = c("blue", "red","green"), ylab = "Transcriptional noise", xaxt = 'n')
axis(1, at = seq(from = 1.5, to = 27.5, by = 3), names(res)[ord], las = 2)
dev.off()

