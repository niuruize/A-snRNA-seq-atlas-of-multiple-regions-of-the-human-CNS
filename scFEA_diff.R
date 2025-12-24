{
  Astro <- celltype_scFEA_mean[2,]-celltype_scFEA_mean[1,]
  rownames(Astro) <- 'Astro'
  CR <- celltype_scFEA_mean[5,]-celltype_scFEA_mean[4,]
  rownames(CR) <- 'CR'
  Endo <- celltype_scFEA_mean[8,]-celltype_scFEA_mean[7,]
  rownames(Endo) <- 'Endo'
  Epend <- celltype_scFEA_mean[11,]-celltype_scFEA_mean[10,]
  rownames(Epend) <- 'Epend'
  ExN <- celltype_scFEA_mean[14,]-celltype_scFEA_mean[13,]
  rownames(ExN) <- 'ExN'
  InN <- celltype_scFEA_mean[17,]-celltype_scFEA_mean[16,]
  rownames(InN) <- 'InN'
  Micro <- celltype_scFEA_mean[20,]-celltype_scFEA_mean[19,]
  rownames(Micro) <- 'Micro'
  MOL <- celltype_scFEA_mean[23,]-celltype_scFEA_mean[22,]
  rownames(MOL) <- 'MOL'
  NFOL <- celltype_scFEA_mean[26,]-celltype_scFEA_mean[25,]
  rownames(NFOL) <- 'NFOL'
  OPC <- celltype_scFEA_mean[29,]-celltype_scFEA_mean[28,]
  rownames(OPC) <- 'OPC'
  Peri <- celltype_scFEA_mean[32,]-celltype_scFEA_mean[31,]
  rownames(Peri) <- 'Peri'
  Tcell <- celltype_scFEA_mean[35,]-celltype_scFEA_mean[34,]
  rownames(Tcell) <- 'Tcell'
  
  scFEA_dif <- rbind(Astro,CR,Endo,Epend,ExN,InN,Micro,MOL,NFOL,OPC,Peri,Tcell)
}

head(scFEA_dif)
write.csv(scFEA_dif,"scFEA_all_diff.csv", na="0") #保存结果

{
  Astro <- celltype_scFEA_mean[3,]-celltype_scFEA_mean[2,]
  rownames(Astro) <- 'Astro'
  CR <- celltype_scFEA_mean[6,]-celltype_scFEA_mean[5,]
  rownames(CR) <- 'CR'
  Endo <- celltype_scFEA_mean[9,]-celltype_scFEA_mean[8,]
  rownames(Endo) <- 'Endo'
  Epend <- celltype_scFEA_mean[12,]-celltype_scFEA_mean[11,]
  rownames(Epend) <- 'Epend'
  ExN <- celltype_scFEA_mean[15,]-celltype_scFEA_mean[14,]
  rownames(ExN) <- 'ExN'
  InN <- celltype_scFEA_mean[18,]-celltype_scFEA_mean[17,]
  rownames(InN) <- 'InN'
  Micro <- celltype_scFEA_mean[21,]-celltype_scFEA_mean[20,]
  rownames(Micro) <- 'Micro'
  MOL <- celltype_scFEA_mean[24,]-celltype_scFEA_mean[23,]
  rownames(MOL) <- 'MOL'
  NFOL <- celltype_scFEA_mean[27,]-celltype_scFEA_mean[26,]
  rownames(NFOL) <- 'NFOL'
  OPC <- celltype_scFEA_mean[30,]-celltype_scFEA_mean[29,]
  rownames(OPC) <- 'OPC'
  Peri <- celltype_scFEA_mean[33,]-celltype_scFEA_mean[32,]
  rownames(Peri) <- 'Peri'
  Tcell <- celltype_scFEA_mean[36,]-celltype_scFEA_mean[35,]
  rownames(Tcell) <- 'Tcell'
  
  scFEA_dif <- rbind(Astro,CR,Endo,Epend,ExN,InN,Micro,MOL,NFOL,OPC,Peri,Tcell)
}

head(scFEA_dif)
write.csv(scFEA_dif,"scFEA_Longevity_Aging_diff.csv", na="0") #保存结果

