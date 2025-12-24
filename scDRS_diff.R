library(data.table);library(hash);library(ggraph);library(ggplot2)

##
{
AD <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/AD.score.gz",header=T,sep="\t")
AD_score <- data.frame(AD=AD$norm_score)
rownames(AD_score) <- AD$V1
scRNA <- AddMetaData(scRNA, metadata = AD_score)
ADHD <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/ADHD.score.gz",header=T,sep="\t")
ADHD_score <- data.frame(ADHD=ADHD$norm_score)
rownames(ADHD_score) <- ADHD$V1
scRNA <- AddMetaData(scRNA, metadata = ADHD_score)
BIP <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/BIP.score.gz",header=T,sep="\t")
BIP_score <- data.frame(BIP=BIP$norm_score)
rownames(BIP_score) <- BIP$V1
scRNA <- AddMetaData(scRNA, metadata = BIP_score)
INT <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/INT.score.gz",header=T,sep="\t")
INT_score <- data.frame(INT=INT$norm_score)
rownames(INT_score) <- INT$V1
scRNA <- AddMetaData(scRNA, metadata = INT_score)
MDD <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/MDD.score.gz",header=T,sep="\t")
MDD_score <- data.frame(MDD=MDD$norm_score)
rownames(MDD_score) <- MDD$V1
scRNA <- AddMetaData(scRNA, metadata = MDD_score)
MS <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/MS.score.gz",header=T,sep="\t")
MS_score <- data.frame(MS=MS$norm_score)
rownames(MS_score) <- MS$V1
scRNA <- AddMetaData(scRNA, metadata = MS_score)
NRT <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/NRT.score.gz",header=T,sep="\t")
NRT_score <- data.frame(NRT=NRT$norm_score)
rownames(NRT_score) <- NRT$V1
scRNA <- AddMetaData(scRNA, metadata = NRT_score)
SCZ <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/SCZ.score.gz",header=T,sep="\t")
SCZ_score <- data.frame(SCZ=SCZ$norm_score)
rownames(SCZ_score) <- SCZ$V1
scRNA <- AddMetaData(scRNA, metadata = SCZ_score)
VNR <- fread("/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/data/VNR.score.gz",header=T,sep="\t")
VNR_score <- data.frame(VNR=VNR$norm_score)
rownames(VNR_score) <- VNR$V1
scRNA <- AddMetaData(scRNA, metadata = VNR_score)
}
##p.value
scRNA1 <- scRNA1[,scRNA1$Group %in% c("Adult","Aging")]
{
  {
    #
    data <- FetchData(scRNA1, vars = c("Group","AD","celltype"))
    p_ME <- compare_means(AD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    AD <- data.frame(celltype=p_ME$celltype,AD=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","ADHD","celltype"))
    p_ME <- compare_means(ADHD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    ADHD <- data.frame(celltype=p_ME$celltype,ADHD=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","BIP","celltype"))
    p_ME <- compare_means(BIP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    BIP <- data.frame(celltype=p_ME$celltype,BIP=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","INT","celltype"))
    p_ME <- compare_means(INT~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    INT <- data.frame(celltype=p_ME$celltype,INT=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","MDD","celltype"))
    p_ME <- compare_means(MDD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    MDD <- data.frame(celltype=p_ME$celltype,MDD=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","MS","celltype"))
    p_ME <- compare_means(MS~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    MS <- data.frame(celltype=p_ME$celltype,MS=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","NRT","celltype"))
    p_ME <- compare_means(NRT~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    NRT <- data.frame(celltype=p_ME$celltype,NRT=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","SCZ","celltype"))
    p_ME <- compare_means(SCZ~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    SCZ <- data.frame(celltype=p_ME$celltype,SCZ=p_ME$p.adj)
    #
    data <- FetchData(scRNA1, vars = c("Group","VNR","celltype"))
    p_ME <- compare_means(VNR~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    VNR <- data.frame(celltype=p_ME$celltype,VNR=p_ME$p.adj)
  }
  {
    p_ME_all <- merge(AD,  ADHD, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, BIP, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, INT, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, MDD, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, MS, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, NRT, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, SCZ, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, VNR, by="celltype",all=TRUE)
  }
  head(p_ME_all)
  write.csv(p_ME_all,"p_score_all.csv", na="0",row.names = F) #保存结果
}
##-logp.value
{
  {
    #
    data <- FetchData(scRNA1, vars = c("Group","AD","celltype"))
    p_ME <- compare_means(AD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    AD <- data.frame(celltype=p_ME$celltype,AD=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","ADHD","celltype"))
    p_ME <- compare_means(ADHD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    ADHD <- data.frame(celltype=p_ME$celltype,ADHD=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","BIP","celltype"))
    p_ME <- compare_means(BIP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    BIP <- data.frame(celltype=p_ME$celltype,BIP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","INT","celltype"))
    p_ME <- compare_means(INT~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    INT <- data.frame(celltype=p_ME$celltype,INT=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","MDD","celltype"))
    p_ME <- compare_means(MDD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    MDD <- data.frame(celltype=p_ME$celltype,MDD=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","MS","celltype"))
    p_ME <- compare_means(MS~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    MS <- data.frame(celltype=p_ME$celltype,MS=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","NRT","celltype"))
    p_ME <- compare_means(NRT~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    NRT <- data.frame(celltype=p_ME$celltype,NRT=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","SCZ","celltype"))
    p_ME <- compare_means(SCZ~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    SCZ <- data.frame(celltype=p_ME$celltype,SCZ=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","VNR","celltype"))
    p_ME <- compare_means(VNR~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    VNR <- data.frame(celltype=p_ME$celltype,VNR=-log10(p_ME$p.adj))
  }
  {
    p_ME_all <- merge(AD,  ADHD, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, BIP, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, INT, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, MDD, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, MS, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, NRT, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, SCZ, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, VNR, by="celltype",all=TRUE)
  }
  head(p_ME_all)
  write.csv(p_ME_all,"logp_score_all.csv", na="0",row.names = F) #保存结果
}

##
{
  scoreID <- c("AD","ADHD","SCZ","VNR","BIP","INT","MDD","MS","NRT")
  scRNA1$celltype_Group <- paste(scRNA1$celltype, scRNA1$Group, sep = "_")
  score<-FetchData(scRNA1, vars = c("celltype_Group",scoreID))
  score<-na.omit(score)
  score_mean <- aggregate(score[,2:10], by=list(type=score$celltype_Group),mean)
  rownames(score_mean) <- score_mean$type
  score_mean <- score_mean[,-1]
  score_mean <- na.omit(score_mean)
  mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
  p1 = pheatmap(t(score_mean), show_colnames = T,scale = "row", cluster_cols = F,treeheight_col = 10, treeheight_row=15,color= mycol)
}

##up_down
{
  Astro <- score_mean[2,]-score_mean[1,]
  rownames(Astro) <- 'Astro'
  CR <- score_mean[4,]-score_mean[3,]
  rownames(CR) <- 'CR'
  Endo <- score_mean[6,]-score_mean[5,]
  rownames(Endo) <- 'Endo'
  Epend <- score_mean[8,]-score_mean[7,]
  rownames(Epend) <- 'Epend'
  ExN <- score_mean[10,]-score_mean[9,]
  rownames(ExN) <- 'ExN'
  InN <- score_mean[12,]-score_mean[11,]
  rownames(InN) <- 'InN'
  Micro <- score_mean[14,]-score_mean[13,]
  rownames(Micro) <- 'Micro'
  MOL <- score_mean[16,]-score_mean[15,]
  rownames(MOL) <- 'MOL'
  NFOL <- score_mean[18,]-score_mean[17,]
  rownames(NFOL) <- 'NFOL'
  OPC <- score_mean[20,]-score_mean[19,]
  rownames(OPC) <- 'OPC'
  Peri <- score_mean[22,]-score_mean[21,]
  rownames(Peri) <- 'Peri'
  Tcell <- score_mean[24,]-score_mean[23,]
  rownames(Tcell) <- 'Tcell'
  
  score_dif <- rbind(Astro,CR,Endo,Epend,ExN,InN,Micro,MOL,NFOL,OPC,Peri,Tcell)
  head(score_dif)
  write.csv(score_dif,"score_all_diff.csv", na="0") #保存结果
}

head(score_dif)
write.csv(score_dif,"score_all_diff.csv", na="0") #保存结果


library(pheatmap)
r2 <- read.csv('/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/logp_score_all.csv',header=TRUE,row.names=1)
p2 <- read.csv('/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/p_score_all.csv',header=TRUE,row.names=1)
p3 <- read.csv('/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/15_scDRS/score_all_diff.csv',header=TRUE,row.names=1)

#使用显著性星号标记进行替换；
p2[p2 < 0.001] <- "***"
p2[p2 >= 0.001 & p2 < 0.01] <- "**"
p2[p2 >= 0.01 & p2 < 0.05] <- "*"
p2[p2 > 0.05] <- ""

#
p3[p3 > 0] <- 1
p3[p3 < 0] <- -1

#
r2 <- as.matrix(r2)
r2[is.infinite(r2)] <- 30
r2[r2 >= 30] <- 30
r2 <- r2*p3
#mycol<-colorRampPalette(c( "#0f86a9", "white", "#ed8b10"))(200)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)

#绘制热图；
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=9,cellwidth=16,
              cellheight=16,cluster_rows=T,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "Hip_score1.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')
p1 = pheatmap(t(r2),scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=F,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= t(p2),show_rownames=T) 
ggsave(filename = "Hip_score2.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')



##umap可视化
scRNA$SCZ<-as.numeric(scRNA$SCZ)
ggplot(data.frame(scRNA@meta.data, scRNA@reductions$umap@cell.embeddings), aes(UMAP_1, UMAP_2, color=SCZ)) +
  geom_point(size=1.5) + scale_color_viridis(option="H")+
  theme_light(base_size = 15)+labs(title = "Breast cancer score")+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+theme(plot.title = element_text(hjust = 0.5))
ggsave("scDRS_result/umap_of_breast_cancer_score.pdf",height=6,width=8)

##
#library(RColorBrewer)
#display.brewer.all()
