
library(Seurat);library(tidyverse);library(cowplot);library(Matrix);library(readxl);library(ggpubr);library(pheatmap);library(ggpubr)

table(scRNA$Group)
scRNA1 = scRNA[,scRNA$Group %in% c("Adult","Aging")]
scRNA1$Group <- as.factor(as.character(scRNA1$Group))
table(scRNA1$Group,scRNA1$celltype)
##score
#load("~/Downloads/scRNA/EC/EC/6_score/PCD_score_genesets.RData")
{
scRNA1 <- AddModuleScore(scRNA1, features = markers[1], ctrl = 100, name = "Apoptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[2], ctrl = 100, name = "Pyroptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[3], ctrl = 100, name = "Ferroptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[4], ctrl = 100, name = "Autophagy")
scRNA1 <- AddModuleScore(scRNA1, features = markers[5], ctrl = 100, name = "Necroptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[6], ctrl = 100, name = "Parthanatos")
scRNA1 <- AddModuleScore(scRNA1, features = markers[7], ctrl = 100, name = "Entotic")
scRNA1 <- AddModuleScore(scRNA1, features = markers[8], ctrl = 100, name = "Netotic")
scRNA1 <- AddModuleScore(scRNA1, features = markers[9], ctrl = 100, name = "Cuproptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[10], ctrl = 100, name = "Lysome_dependent")
scRNA1 <- AddModuleScore(scRNA1, features = markers[11], ctrl = 100, name = "Alkaliptosis")
scRNA1 <- AddModuleScore(scRNA1, features = markers[12], ctrl = 100, name = "Oxeiprosis")
}

##celltype
{
  score_ID <- data.frame(scoreID=paste(names(markers), "1", sep = ""))
  score <- FetchData(scRNA1, vars = c("celltype",score_ID$scoreID))
  score_mean <- aggregate(score[,2:13], by=list(type=score$celltype),mean)
  rownames(score_mean) <- score_mean$type
  score_mean <- score_mean[,-1]
  score_mean <- na.omit(score_mean)
  mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
  p1 = pheatmap(score_mean,scale = "none", border_color= "grey", number_color= "black",
                fontsize_number=14,fontsize_row=8,fontsize_col=9,cellwidth=14,
                cellheight=14,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
                color= mycol,show_rownames=T)
  ggsave(filename = "PFC_score_celltype.pdf", p1,device = 'pdf', width = 10, height = 20, units = 'cm')
}

##P.value
{
{
  #
  data <- FetchData(scRNA1, vars = c("Group","Apoptosis1","celltype"))
  p_ME <- compare_means(Apoptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Apoptosis <- data.frame(celltype=p_ME$celltype,Apoptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Pyroptosis1","celltype"))
  p_ME <- compare_means(Pyroptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Pyroptosis <- data.frame(celltype=p_ME$celltype,Pyroptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Ferroptosis1","celltype"))
  p_ME <- compare_means(Ferroptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Ferroptosis <- data.frame(celltype=p_ME$celltype,Ferroptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Autophagy1","celltype"))
  p_ME <- compare_means(Autophagy1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Autophagy <- data.frame(celltype=p_ME$celltype,Autophagy=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Necroptosis1","celltype"))
  p_ME <- compare_means(Necroptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Necroptosis <- data.frame(celltype=p_ME$celltype,Necroptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Cuproptosis1","celltype"))
  p_ME <- compare_means(Cuproptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Cuproptosis <- data.frame(celltype=p_ME$celltype,Cuproptosis=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Parthanatos1","celltype"))
  p_ME <- compare_means(Parthanatos1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Parthanatos <- data.frame(celltype=p_ME$celltype,Parthanatos=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Entotic1","celltype"))
  p_ME <- compare_means(Entotic1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Entotic <- data.frame(celltype=p_ME$celltype,Entotic=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Netotic1","celltype"))
  p_ME <- compare_means(Netotic1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Netotic <- data.frame(celltype=p_ME$celltype,Netotic=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Lysome_dependent1","celltype"))
  p_ME <- compare_means(Lysome_dependent1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Lysome_dependent <- data.frame(celltype=p_ME$celltype,Lysome_dependent=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Alkaliptosis1","celltype"))
  p_ME <- compare_means(Alkaliptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Alkaliptosis <- data.frame(celltype=p_ME$celltype,Alkaliptosis1=p_ME$p.adj)
  #
  data <- FetchData(scRNA1, vars = c("Group","Oxeiprosis1","celltype"))
  p_ME <- compare_means(Oxeiprosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  Oxeiprosis <- data.frame(celltype=p_ME$celltype,Oxeiprosis=p_ME$p.adj)
}
{
  p_ME_all <- merge(Apoptosis,  Pyroptosis, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Ferroptosis, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Autophagy, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Necroptosis, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Cuproptosis, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Parthanatos, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Entotic, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Netotic, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Lysome_dependent, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Alkaliptosis, by="celltype",all=TRUE)
  p_ME_all <- merge(p_ME_all, Oxeiprosis, by="celltype",all=TRUE)
}
head(p_ME_all)
write.csv(p_ME_all,"p_score_all.csv", na="0",row.names = F) #保存结果
}
##-logP.value
{
  {
    #
    data <- FetchData(scRNA1, vars = c("Group","Apoptosis1","celltype"))
    p_ME <- compare_means(Apoptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Apoptosis <- data.frame(celltype=p_ME$celltype,Apoptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Pyroptosis1","celltype"))
    p_ME <- compare_means(Pyroptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Pyroptosis <- data.frame(celltype=p_ME$celltype,Pyroptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Ferroptosis1","celltype"))
    p_ME <- compare_means(Ferroptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Ferroptosis <- data.frame(celltype=p_ME$celltype,Ferroptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Autophagy1","celltype"))
    p_ME <- compare_means(Autophagy1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Autophagy <- data.frame(celltype=p_ME$celltype,Autophagy=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Necroptosis1","celltype"))
    p_ME <- compare_means(Necroptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Necroptosis <- data.frame(celltype=p_ME$celltype,Necroptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Cuproptosis1","celltype"))
    p_ME <- compare_means(Cuproptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Cuproptosis <- data.frame(celltype=p_ME$celltype,Cuproptosis=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Parthanatos1","celltype"))
    p_ME <- compare_means(Parthanatos1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Parthanatos <- data.frame(celltype=p_ME$celltype,Parthanatos=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Entotic1","celltype"))
    p_ME <- compare_means(Entotic1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Entotic <- data.frame(celltype=p_ME$celltype,Entotic=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Netotic1","celltype"))
    p_ME <- compare_means(Netotic1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Netotic <- data.frame(celltype=p_ME$celltype,Netotic=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Lysome_dependent1","celltype"))
    p_ME <- compare_means(Lysome_dependent1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Lysome_dependent <- data.frame(celltype=p_ME$celltype,Lysome_dependent=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Alkaliptosis1","celltype"))
    p_ME <- compare_means(Alkaliptosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Alkaliptosis <- data.frame(celltype=p_ME$celltype,Alkaliptosis1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA1, vars = c("Group","Oxeiprosis1","celltype"))
    p_ME <- compare_means(Oxeiprosis1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    Oxeiprosis <- data.frame(celltype=p_ME$celltype,Oxeiprosis=-log10(p_ME$p.adj))
  }
  {
    p_ME_all <- merge(Apoptosis,  Pyroptosis, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Ferroptosis, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Autophagy, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Necroptosis, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Cuproptosis, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Parthanatos, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Entotic, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Netotic, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Lysome_dependent, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Alkaliptosis, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, Oxeiprosis, by="celltype",all=TRUE)
  }
  head(p_ME_all)
  write.csv(p_ME_all,"logp_score_all.csv", na="0",row.names = F) #保存结果
}

##up_down
{
scRNA1$celltype_Group <- paste(scRNA1$celltype, scRNA1$Group, sep = "_")
score_ID <- data.frame(scoreID=paste(names(markers), "1", sep = ""))
score<-FetchData(scRNA1, vars = c("celltype_Group",score_ID$scoreID))
score_mean <- aggregate(score[,2:13], by=list(type=score$celltype_Group),mean)
rownames(score_mean) <- score_mean$type
score_mean <- score_mean[,-1]
score_mean <- na.omit(score_mean)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap(t(score_mean), show_colnames = T,scale = "row", cluster_cols = F,treeheight_col = 10, treeheight_row=15,color= mycol)
}

gc()
##
rownames(score_mean)
{
  Astro <- score_mean[2,]-score_mean[1,]
  rownames(Astro) <- 'Astro'
  CR <- score_mean[4,]-score_mean[3,]
  rownames(CR) <- 'CR'
  Endo <- score_mean[6,]-score_mean[5,]
  rownames(Endo) <- 'Endo'
  ExN <- score_mean[8,]-score_mean[7,]
  rownames(ExN) <- 'ExN'
  InN <- score_mean[10,]-score_mean[9,]
  rownames(InN) <- 'InN'
  Micro <- score_mean[12,]-score_mean[11,]
  rownames(Micro) <- 'Micro'
  MOL <- score_mean[14,]-score_mean[13,]
  rownames(MOL) <- 'MOL'
  OPC <- score_mean[16,]-score_mean[15,]
  rownames(OPC) <- 'OPC'
  Peri <- score_mean[18,]-score_mean[17,]
  rownames(Peri) <- 'Peri'
  scFEA_dif <- rbind(Astro,CR,Endo,ExN,InN,Micro,MOL,OPC,Peri)
}

head(scFEA_dif)
write.csv(scFEA_dif,"scFEA_all_diff.csv", na="0",row.names = T) #保存结果


###
library(pheatmap)
r2 <- read.csv('logp_score_all.csv',header=TRUE,row.names=1)
p2 <- read.csv('p_score_all.csv',header=TRUE,row.names=1)
p3 <- read.csv('scFEA_all_diff.csv',header=TRUE,row.names=1)

#使用显著性星号标记进行替换；
p2[p2 < 0.001] <- "***"
p2[p2 >= 0.001 & p2 < 0.01] <- "**"
p2[p2 >= 0.01 & p2 < 0.05] <- "*"
p2[p2 >= 0.05] <- ""

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
              fontsize_number=12,fontsize_row=8,fontsize_col=9,cellwidth=15,
              cellheight=15,cluster_rows=T,cluster_cols=T,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "PCD1.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')
p1 = pheatmap(r2,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=12,fontsize_row=8,fontsize_col=9,cellwidth=16,
              cellheight=16,cluster_rows=T,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "PCD2.pdf", p1,device = 'pdf', width = 15, height = 15, units = 'cm')


