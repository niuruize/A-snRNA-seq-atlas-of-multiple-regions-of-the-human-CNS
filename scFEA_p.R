library(ggpubr)
predME <- colnames(predFlux)
table(scRNA$Group)
scRNA <- scRNA[,scRNA$Group %in% c("Adult","Aging")]

{
{
#
data <- FetchData(scRNA, vars = c("Group",predME[1],"celltype"))
p_ME <- compare_means(AMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME1 <- data.frame(celltype=p_ME$celltype,AMP=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[2],"celltype"))
p_ME <- compare_means(Pyruvate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME2 <- data.frame(celltype=p_ME$celltype,Pyruvate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[3],"celltype"))
p_ME <- compare_means(Acetyl.CoA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME3 <- data.frame(celltype=p_ME$celltype,Acetyl.CoA=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[4],"celltype"))
p_ME <- compare_means(Glutamate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME4 <- data.frame(celltype=p_ME$celltype,Glutamate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[5],"celltype"))
p_ME <- compare_means(X2OG~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME5 <- data.frame(celltype=p_ME$celltype,X2OG=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[6],"celltype"))
p_ME <- compare_means(Oxaloacetate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME6 <- data.frame(celltype=p_ME$celltype,Oxaloacetate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[7],"celltype"))
p_ME <- compare_means(Glycine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME7 <- data.frame(celltype=p_ME$celltype,Glycine=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[8],"celltype"))
p_ME <- compare_means(Succinate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME8 <- data.frame(celltype=p_ME$celltype,Succinate=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[9],"celltype"))
p_ME <- compare_means(UDP.N.acetylglucosamine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME9 <- data.frame(celltype=p_ME$celltype,UDP.N.acetylglucosamine=p_ME$p.adj)
#
data <- FetchData(scRNA, vars = c("Group",predME[10],"celltype"))
p_ME <- compare_means(lysine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
p_ME10 <- data.frame(celltype=p_ME$celltype,lysine=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[11],"celltype"))
  p_ME <- compare_means(Aspartate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME11 <- data.frame(celltype=p_ME$celltype,Aspartate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[12],"celltype"))
  p_ME <- compare_means(Glutathione~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME12 <- data.frame(celltype=p_ME$celltype,Glutathione=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[13],"celltype"))
  p_ME <- compare_means(Arginine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME13 <- data.frame(celltype=p_ME$celltype,Arginine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[14],"celltype"))
  p_ME <- compare_means(Glutamine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME14 <- data.frame(celltype=p_ME$celltype,Glutamine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[15],"celltype"))
  p_ME <- compare_means(Serine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME15 <- data.frame(celltype=p_ME$celltype,Serine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[16],"celltype"))
  p_ME <- compare_means(Methionine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME16 <- data.frame(celltype=p_ME$celltype,Methionine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[17],"celltype"))
  p_ME <- compare_means(Ornithine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME17 <- data.frame(celltype=p_ME$celltype,Ornithine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[18],"celltype"))
  p_ME <- compare_means(Phenylalanine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME18 <- data.frame(celltype=p_ME$celltype,Phenylalanine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[19],"celltype"))
  p_ME <- compare_means(Tyrosine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME19 <- data.frame(celltype=p_ME$celltype,Tyrosine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[20],"celltype"))
  p_ME <- compare_means(Succinyl.CoA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME20 <- data.frame(celltype=p_ME$celltype,Succinyl.CoA=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[21],"celltype"))
  p_ME <- compare_means(Cysteine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME21 <- data.frame(celltype=p_ME$celltype,Cysteine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[22],"celltype"))
  p_ME <- compare_means(B.Alanine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME22 <- data.frame(celltype=p_ME$celltype,B.Alanine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[23],"celltype"))
  p_ME <- compare_means(Propanoyl.CoA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME23 <- data.frame(celltype=p_ME$celltype,Propanoyl.CoA=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[24],"celltype"))
  p_ME <- compare_means(Glucose.1.phosphate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME24 <- data.frame(celltype=p_ME$celltype,Glucose.1.phosphate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[25],"celltype"))
  p_ME <- compare_means(UMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME25 <- data.frame(celltype=p_ME$celltype,UMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[26],"celltype"))
  p_ME <- compare_means(Dolichyl.phosphate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME26 <- data.frame(celltype=p_ME$celltype,Dolichyl.phosphate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[27],"celltype"))
  p_ME <- compare_means(CDP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME27 <- data.frame(celltype=p_ME$celltype,CDP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[28],"celltype"))
  p_ME <- compare_means(Choline~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME28 <- data.frame(celltype=p_ME$celltype,Choline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[29],"celltype"))
  p_ME <- compare_means(G3P~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME29 <- data.frame(celltype=p_ME$celltype,G3P=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[30],"celltype"))
  p_ME <- compare_means(PRPP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME30 <- data.frame(celltype=p_ME$celltype,PRPP=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[31],"celltype"))
  p_ME <- compare_means(Fumarate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME31 <- data.frame(celltype=p_ME$celltype,Fumarate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[32],"celltype"))
  p_ME <- compare_means(Leucine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME32 <- data.frame(celltype=p_ME$celltype,Leucine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[33],"celltype"))
  p_ME <- compare_means(IMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME33 <- data.frame(celltype=p_ME$celltype,IMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[34],"celltype"))
  p_ME <- compare_means(Putrescine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME34 <- data.frame(celltype=p_ME$celltype,Putrescine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[35],"celltype"))
  p_ME <- compare_means(Histidine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME35 <- data.frame(celltype=p_ME$celltype,Histidine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[36],"celltype"))
  p_ME <- compare_means(GMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME36 <- data.frame(celltype=p_ME$celltype,GMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[37],"celltype"))
  p_ME <- compare_means(Proline~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME37 <- data.frame(celltype=p_ME$celltype,Proline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[38],"celltype"))
  p_ME <- compare_means(Malate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME38 <- data.frame(celltype=p_ME$celltype,Malate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[39],"celltype"))
  p_ME <- compare_means(Citrate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME39 <- data.frame(celltype=p_ME$celltype,Citrate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[40],"celltype"))
  p_ME <- compare_means(UDP.glucuronic.acid~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME40 <- data.frame(celltype=p_ME$celltype,UDP.glucuronic.acid=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[41],"celltype"))
  p_ME <- compare_means(Valine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME41 <- data.frame(celltype=p_ME$celltype,Valine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[42],"celltype"))
  p_ME <- compare_means(Cholesterol~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME42 <- data.frame(celltype=p_ME$celltype,Cholesterol=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[43],"celltype"))
  p_ME <- compare_means(Threonine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME43 <- data.frame(celltype=p_ME$celltype,Threonine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[44],"celltype"))
  p_ME <- compare_means(X3PD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME44 <- data.frame(celltype=p_ME$celltype,X3PD=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[45],"celltype"))
  p_ME <- compare_means(dCMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME45 <- data.frame(celltype=p_ME$celltype,dCMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[46],"celltype"))
  p_ME <- compare_means(Lactate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME46 <- data.frame(celltype=p_ME$celltype,Lactate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[47],"celltype"))
  p_ME <- compare_means(Hypoxanthine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME47 <- data.frame(celltype=p_ME$celltype,Hypoxanthine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[48],"celltype"))
  p_ME <- compare_means(Glucose~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME48 <- data.frame(celltype=p_ME$celltype,Glucose=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[49],"celltype"))
  p_ME <- compare_means(Citrulline~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME49 <- data.frame(celltype=p_ME$celltype,Citrulline=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[50],"celltype"))
  p_ME <- compare_means(GABA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME50 <- data.frame(celltype=p_ME$celltype,GABA=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[51],"celltype"))
  p_ME <- compare_means(dTMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME51 <- data.frame(celltype=p_ME$celltype,dTMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[52],"celltype"))
  p_ME <- compare_means(dUMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME52 <- data.frame(celltype=p_ME$celltype,dUMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[53],"celltype"))
  p_ME <- compare_means(Xanthine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME53 <- data.frame(celltype=p_ME$celltype,Xanthine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[54],"celltype"))
  p_ME <- compare_means(Chondroitin~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME54 <- data.frame(celltype=p_ME$celltype,Chondroitin=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[55],"celltype"))
  p_ME <- compare_means(Isoleucine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME55 <- data.frame(celltype=p_ME$celltype,Isoleucine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[56],"celltype"))
  p_ME <- compare_means(X.E.E..Farnesyl.PP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME56 <- data.frame(celltype=p_ME$celltype,X.E.E..Farnesyl.PP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[57],"celltype"))
  p_ME <- compare_means(Deoxyadenosine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME57 <- data.frame(celltype=p_ME$celltype,Deoxyadenosine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[58],"celltype"))
  p_ME <- compare_means(XMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME58 <- data.frame(celltype=p_ME$celltype,XMP=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[59],"celltype"))
  p_ME <- compare_means(G6P~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME59 <- data.frame(celltype=p_ME$celltype,G6P=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[60],"celltype"))
  p_ME <- compare_means(dCDP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME60 <- data.frame(celltype=p_ME$celltype,dCDP=p_ME$p.adj)
}
{
  #
  data <- FetchData(scRNA, vars = c("Group",predME[61],"celltype"))
  p_ME <- compare_means(Spermine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME61 <- data.frame(celltype=p_ME$celltype,Spermine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[62],"celltype"))
  p_ME <- compare_means(Argininosuccinate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME62 <- data.frame(celltype=p_ME$celltype,Argininosuccinate=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[63],"celltype"))
  p_ME <- compare_means(Dolichyl.phosphate.D.mannose~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME63 <- data.frame(celltype=p_ME$celltype,Dolichyl.phosphate.D.mannose=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[64],"celltype"))
  p_ME <- compare_means(AICAR~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME64 <- data.frame(celltype=p_ME$celltype,AICAR=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[65],"celltype"))
  p_ME <- compare_means(Fatty.Acid~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME65 <- data.frame(celltype=p_ME$celltype,Fatty.Acid=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[66],"celltype"))
  p_ME <- compare_means(Pyrimidine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME66 <- data.frame(celltype=p_ME$celltype,Pyrimidine=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[67],"celltype"))
  p_ME <- compare_means(X.Glc.3..GlcNAc.2..Man.9..Asn.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME67 <- data.frame(celltype=p_ME$celltype,X.Glc.3..GlcNAc.2..Man.9..Asn.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[68],"celltype"))
  p_ME <- compare_means(X.GlcNAc.4..Man.3..Asn.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME68 <- data.frame(celltype=p_ME$celltype,X.GlcNAc.4..Man.3..Asn.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[69],"celltype"))
  p_ME <- compare_means(X.Gal.2..GlcA.1..Xyl.1..Ser.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME69 <- data.frame(celltype=p_ME$celltype,X.Gal.2..GlcA.1..Xyl.1..Ser.1=p_ME$p.adj)
  #
  data <- FetchData(scRNA, vars = c("Group",predME[70],"celltype"))
  p_ME <- compare_means(X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
  p_ME70 <- data.frame(celltype=p_ME$celltype,X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1=p_ME$p.adj)
}

{
p_ME_all <- merge(p_ME1,  p_ME2, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME3, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME4, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME5, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME6, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME7, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME8, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME9, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME10, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME11, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME12, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME13, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME14, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME15, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME16, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME17, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME18, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME19, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME20, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME21, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME22, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME23, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME24, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME25, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME26, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME27, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME28, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME29, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME30, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME31, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME32, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME33, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME34, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME35, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME36, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME37, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME38, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME39, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME40, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME41, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME42, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME43, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME44, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME45, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME46, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME47, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME48, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME49, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME50, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME51, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME52, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME53, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME54, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME55, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME56, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME57, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME58, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME59, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME60, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME61, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME62, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME63, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME64, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME65, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME66, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME67, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME68, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME69, by="celltype",all=TRUE)
p_ME_all <- merge(p_ME_all, p_ME70, by="celltype",all=TRUE)
}

head(p_ME_all)

write.csv(p_ME_all,"p_ME_all.csv", na="0",row.names = F) #保存结果
}

{
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[1],"celltype"))
    p_ME <- compare_means(AMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME1 <- data.frame(celltype=p_ME$celltype,AMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[2],"celltype"))
    p_ME <- compare_means(Pyruvate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME2 <- data.frame(celltype=p_ME$celltype,Pyruvate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[3],"celltype"))
    p_ME <- compare_means(Acetyl.CoA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME3 <- data.frame(celltype=p_ME$celltype,Acetyl.CoA=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[4],"celltype"))
    p_ME <- compare_means(Glutamate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME4 <- data.frame(celltype=p_ME$celltype,Glutamate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[5],"celltype"))
    p_ME <- compare_means(X2OG~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME5 <- data.frame(celltype=p_ME$celltype,X2OG=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[6],"celltype"))
    p_ME <- compare_means(Oxaloacetate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME6 <- data.frame(celltype=p_ME$celltype,Oxaloacetate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[7],"celltype"))
    p_ME <- compare_means(Glycine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME7 <- data.frame(celltype=p_ME$celltype,Glycine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[8],"celltype"))
    p_ME <- compare_means(Succinate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME8 <- data.frame(celltype=p_ME$celltype,Succinate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[9],"celltype"))
    p_ME <- compare_means(UDP.N.acetylglucosamine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME9 <- data.frame(celltype=p_ME$celltype,UDP.N.acetylglucosamine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[10],"celltype"))
    p_ME <- compare_means(lysine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME10 <- data.frame(celltype=p_ME$celltype,lysine=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[11],"celltype"))
    p_ME <- compare_means(Aspartate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME11 <- data.frame(celltype=p_ME$celltype,Aspartate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[12],"celltype"))
    p_ME <- compare_means(Glutathione~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME12 <- data.frame(celltype=p_ME$celltype,Glutathione=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[13],"celltype"))
    p_ME <- compare_means(Arginine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME13 <- data.frame(celltype=p_ME$celltype,Arginine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[14],"celltype"))
    p_ME <- compare_means(Glutamine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME14 <- data.frame(celltype=p_ME$celltype,Glutamine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[15],"celltype"))
    p_ME <- compare_means(Serine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME15 <- data.frame(celltype=p_ME$celltype,Serine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[16],"celltype"))
    p_ME <- compare_means(Methionine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME16 <- data.frame(celltype=p_ME$celltype,Methionine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[17],"celltype"))
    p_ME <- compare_means(Ornithine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME17 <- data.frame(celltype=p_ME$celltype,Ornithine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[18],"celltype"))
    p_ME <- compare_means(Phenylalanine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME18 <- data.frame(celltype=p_ME$celltype,Phenylalanine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[19],"celltype"))
    p_ME <- compare_means(Tyrosine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME19 <- data.frame(celltype=p_ME$celltype,Tyrosine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[20],"celltype"))
    p_ME <- compare_means(Succinyl.CoA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME20 <- data.frame(celltype=p_ME$celltype,Succinyl.CoA=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[21],"celltype"))
    p_ME <- compare_means(Cysteine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME21 <- data.frame(celltype=p_ME$celltype,Cysteine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[22],"celltype"))
    p_ME <- compare_means(B.Alanine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME22 <- data.frame(celltype=p_ME$celltype,B.Alanine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[23],"celltype"))
    p_ME <- compare_means(Propanoyl.CoA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME23 <- data.frame(celltype=p_ME$celltype,Propanoyl.CoA=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[24],"celltype"))
    p_ME <- compare_means(Glucose.1.phosphate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME24 <- data.frame(celltype=p_ME$celltype,Glucose.1.phosphate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[25],"celltype"))
    p_ME <- compare_means(UMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME25 <- data.frame(celltype=p_ME$celltype,UMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[26],"celltype"))
    p_ME <- compare_means(Dolichyl.phosphate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME26 <- data.frame(celltype=p_ME$celltype,Dolichyl.phosphate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[27],"celltype"))
    p_ME <- compare_means(CDP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME27 <- data.frame(celltype=p_ME$celltype,CDP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[28],"celltype"))
    p_ME <- compare_means(Choline~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME28 <- data.frame(celltype=p_ME$celltype,Choline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[29],"celltype"))
    p_ME <- compare_means(G3P~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME29 <- data.frame(celltype=p_ME$celltype,G3P=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[30],"celltype"))
    p_ME <- compare_means(PRPP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME30 <- data.frame(celltype=p_ME$celltype,PRPP=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[31],"celltype"))
    p_ME <- compare_means(Fumarate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME31 <- data.frame(celltype=p_ME$celltype,Fumarate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[32],"celltype"))
    p_ME <- compare_means(Leucine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME32 <- data.frame(celltype=p_ME$celltype,Leucine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[33],"celltype"))
    p_ME <- compare_means(IMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME33 <- data.frame(celltype=p_ME$celltype,IMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[34],"celltype"))
    p_ME <- compare_means(Putrescine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME34 <- data.frame(celltype=p_ME$celltype,Putrescine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[35],"celltype"))
    p_ME <- compare_means(Histidine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME35 <- data.frame(celltype=p_ME$celltype,Histidine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[36],"celltype"))
    p_ME <- compare_means(GMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME36 <- data.frame(celltype=p_ME$celltype,GMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[37],"celltype"))
    p_ME <- compare_means(Proline~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME37 <- data.frame(celltype=p_ME$celltype,Proline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[38],"celltype"))
    p_ME <- compare_means(Malate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME38 <- data.frame(celltype=p_ME$celltype,Malate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[39],"celltype"))
    p_ME <- compare_means(Citrate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME39 <- data.frame(celltype=p_ME$celltype,Citrate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[40],"celltype"))
    p_ME <- compare_means(UDP.glucuronic.acid~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME40 <- data.frame(celltype=p_ME$celltype,UDP.glucuronic.acid=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[41],"celltype"))
    p_ME <- compare_means(Valine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME41 <- data.frame(celltype=p_ME$celltype,Valine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[42],"celltype"))
    p_ME <- compare_means(Cholesterol~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME42 <- data.frame(celltype=p_ME$celltype,Cholesterol=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[43],"celltype"))
    p_ME <- compare_means(Threonine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME43 <- data.frame(celltype=p_ME$celltype,Threonine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[44],"celltype"))
    p_ME <- compare_means(X3PD~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME44 <- data.frame(celltype=p_ME$celltype,X3PD=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[45],"celltype"))
    p_ME <- compare_means(dCMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME45 <- data.frame(celltype=p_ME$celltype,dCMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[46],"celltype"))
    p_ME <- compare_means(Lactate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME46 <- data.frame(celltype=p_ME$celltype,Lactate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[47],"celltype"))
    p_ME <- compare_means(Hypoxanthine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME47 <- data.frame(celltype=p_ME$celltype,Hypoxanthine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[48],"celltype"))
    p_ME <- compare_means(Glucose~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME48 <- data.frame(celltype=p_ME$celltype,Glucose=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[49],"celltype"))
    p_ME <- compare_means(Citrulline~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME49 <- data.frame(celltype=p_ME$celltype,Citrulline=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[50],"celltype"))
    p_ME <- compare_means(GABA~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME50 <- data.frame(celltype=p_ME$celltype,GABA=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[51],"celltype"))
    p_ME <- compare_means(dTMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME51 <- data.frame(celltype=p_ME$celltype,dTMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[52],"celltype"))
    p_ME <- compare_means(dUMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME52 <- data.frame(celltype=p_ME$celltype,dUMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[53],"celltype"))
    p_ME <- compare_means(Xanthine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME53 <- data.frame(celltype=p_ME$celltype,Xanthine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[54],"celltype"))
    p_ME <- compare_means(Chondroitin~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME54 <- data.frame(celltype=p_ME$celltype,Chondroitin=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[55],"celltype"))
    p_ME <- compare_means(Isoleucine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME55 <- data.frame(celltype=p_ME$celltype,Isoleucine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[56],"celltype"))
    p_ME <- compare_means(X.E.E..Farnesyl.PP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME56 <- data.frame(celltype=p_ME$celltype,X.E.E..Farnesyl.PP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[57],"celltype"))
    p_ME <- compare_means(Deoxyadenosine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME57 <- data.frame(celltype=p_ME$celltype,Deoxyadenosine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[58],"celltype"))
    p_ME <- compare_means(XMP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME58 <- data.frame(celltype=p_ME$celltype,XMP=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[59],"celltype"))
    p_ME <- compare_means(G6P~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME59 <- data.frame(celltype=p_ME$celltype,G6P=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[60],"celltype"))
    p_ME <- compare_means(dCDP~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME60 <- data.frame(celltype=p_ME$celltype,dCDP=-log10(p_ME$p.adj))
  }
  {
    #
    data <- FetchData(scRNA, vars = c("Group",predME[61],"celltype"))
    p_ME <- compare_means(Spermine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME61 <- data.frame(celltype=p_ME$celltype,Spermine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[62],"celltype"))
    p_ME <- compare_means(Argininosuccinate~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME62 <- data.frame(celltype=p_ME$celltype,Argininosuccinate=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[63],"celltype"))
    p_ME <- compare_means(Dolichyl.phosphate.D.mannose~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME63 <- data.frame(celltype=p_ME$celltype,Dolichyl.phosphate.D.mannose=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[64],"celltype"))
    p_ME <- compare_means(AICAR~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME64 <- data.frame(celltype=p_ME$celltype,AICAR=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[65],"celltype"))
    p_ME <- compare_means(Fatty.Acid~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME65 <- data.frame(celltype=p_ME$celltype,Fatty.Acid=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[66],"celltype"))
    p_ME <- compare_means(Pyrimidine~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME66 <- data.frame(celltype=p_ME$celltype,Pyrimidine=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[67],"celltype"))
    p_ME <- compare_means(X.Glc.3..GlcNAc.2..Man.9..Asn.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME67 <- data.frame(celltype=p_ME$celltype,X.Glc.3..GlcNAc.2..Man.9..Asn.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[68],"celltype"))
    p_ME <- compare_means(X.GlcNAc.4..Man.3..Asn.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME68 <- data.frame(celltype=p_ME$celltype,X.GlcNAc.4..Man.3..Asn.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[69],"celltype"))
    p_ME <- compare_means(X.Gal.2..GlcA.1..Xyl.1..Ser.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME69 <- data.frame(celltype=p_ME$celltype,X.Gal.2..GlcA.1..Xyl.1..Ser.1=-log10(p_ME$p.adj))
    #
    data <- FetchData(scRNA, vars = c("Group",predME[70],"celltype"))
    p_ME <- compare_means(X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1~Group, data, group.by = "celltype",method = "wilcox.test", paired = FALSE)
    p_ME70 <- data.frame(celltype=p_ME$celltype,X.Gal.1..GlcNAc.1..Man.1..Ser.Thr.1=-log10(p_ME$p.adj))
  }
  
  {
    p_ME_all <- merge(p_ME1,  p_ME2, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME3, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME4, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME5, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME6, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME7, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME8, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME9, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME10, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME11, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME12, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME13, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME14, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME15, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME16, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME17, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME18, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME19, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME20, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME21, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME22, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME23, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME24, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME25, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME26, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME27, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME28, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME29, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME30, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME31, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME32, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME33, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME34, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME35, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME36, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME37, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME38, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME39, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME40, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME41, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME42, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME43, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME44, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME45, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME46, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME47, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME48, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME49, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME50, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME51, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME52, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME53, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME54, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME55, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME56, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME57, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME58, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME59, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME60, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME61, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME62, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME63, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME64, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME65, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME66, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME67, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME68, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME69, by="celltype",all=TRUE)
    p_ME_all <- merge(p_ME_all, p_ME70, by="celltype",all=TRUE)
  }
  
  head(p_ME_all)
  
  write.csv(p_ME_all,"logp_ME_all.csv", na="0",row.names = F) #保存结果
}

r2 <- read.csv('/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/14_scFEA/logp_ME_all.csv', header = T, row.names = 1)
r2[r2 >= 15] <- 15
r3 <- read.csv('/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/14_scFEA/scFEA_all_diff.csv', header = T, row.names = 1)
r3[r3 > 0] <- 1
r3[r3 < 0] <- -1
r4 <- r3*r2
p2 <- read.csv('/Users/niuruize/Downloads/scRNA/AD_aging/aging/Hip/Seurat/14_scFEA/p_ME_all.csv', header = T, row.names = 1)
p2[p2 < 0.01] <- "**"
p2[p2 < 0.001] <- "***"
p2[p2 < 0.05 & p2 >= 0.01] <- "*"
p2[p2 >= 0.05] <- ""
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap(r4,scale = "none", border_color= "grey", number_color= "black",
              fontsize_number=10,fontsize_row=8,fontsize_col=9,cellwidth=14,
              cellheight=12,cluster_rows=T,cluster_cols=F,treeheight_row = 10,treeheight_col = 10,
              color= mycol,display_numbers= p2,show_rownames=T) 
ggsave(filename = "p_ME_all1.pdf", p1,device = 'pdf', width = 50, height = 20, units = 'cm')

