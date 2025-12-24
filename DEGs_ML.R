###------------------------------------------------------AutoML-H2O------------------------------------------------------------###
library(h2o)
h2o.init()
library(Seurat);library(dplyr);library(randomForest);library(caret);library(pROC);library(caret);library(caret)
library(modeldata);library(tidymodels);library(recipes) 
# load data
#load("~/aging/EC/EC_Harmony.RData")

# Extracting feature information
table(scRNA$Group)
scRNA = scRNA[,scRNA$Group %in% c("Adult","Aging")]
DefaultAssay(scRNA) <- "RNA"
sce.markers <- read.csv("markers.txt",sep = "")
Idents(scRNA) <- "Group"
#sce <- subset(scRNA, subset = nFeature_RNA > 500)
#sce <- subset(sce,downsample = 200)
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers$gene)) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]

# train data and test data
inTrain <- createDataPartition(y= Idents(sce),p=0.25,list=F)
test_expr <- t_expr[inTrain,]
train_expr <- t_expr[-inTrain,]
test_y <- Idents(sce)[inTrain]
train_y <- Idents(sce)[-inTrain]
save(test_y,train_y, test_expr,train_expr,file = 'input.Rdata') 
train_expr[1:4,1:4]
table(train_y)
table(test_y)
table(scRNA$Group)

# Building data
test_expr = cbind(data.frame(test_expr),data.frame(Group=test_y))
train_expr = cbind(data.frame(train_expr),data.frame(Group=train_y))
#train_expr = cbind(data.frame(t_expr),data.frame(Group=Idents(sce)))
rec <- recipe(Group~.,train_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() # 转换为h2o 数据框
rec$Group <- h2o.asfactor(rec$Group)
y <- "Group"  # 因变量名        
x <- setdiff(names(rec),y)   # 预测变量名

# Building model
am <- h2o.automl(x,y, training_frame = rec,max_models = 20, seed=1234)
save(am,file="h2o_automl.RData")
#am <- h2o.automl(y="Group", training_frame = train_expr, max_models = 10)

#
b <- h2o.get_leaderboard(am) # 默认为排名前6的模型
b
#best <- h2o.get_best_model(am)
best <- h2o.get_best_model(am,"gbm")
perf <- h2o.performance(best)
perf
Confusion_Matrix <- as.matrix(perf@metrics$cm$table[1:2,1:2])
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix1.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')

#测试
test <- recipe(Group~.,test_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() # 转换为h2o 数据框
test$Group <- h2o.asfactor(test$Group)
y <- "Group"  # 因变量名        
x <- setdiff(names(test),y)   # 预测变量名
pred <- h2o.predict(best, test)
## 制作测试结果的混淆矩阵
Confusion_Matrix <- as.matrix(table(cbind(as.data.frame(pred$predict), as.data.frame(test$Group))))
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix_test.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')


#绘制ROC曲线
#AUC
h2o.auc(perf)
## [1] 0.9979307
# 绘制ROC曲线
plot(perf,type="roc")

#变量重要性
h2o.varimp_plot(best)
#基于排列的变量重要性
pdf("permutation_importance.pdf", width = 5,height = 5)
h2o.permutation_importance_plot(best,rec,num_of_features=10)
dev.off()

#模型解释：SHAP summary plot
pdf("SHAP_summary_plot.pdf", width = 6,height = 3)
h2o.shap_summary_plot(best, rec,top_n_features = 10)
dev.off()

#单个样本的SHAP
h2o.shap_explain_row_plot(best, rec, row_index = 2)
#h2o.explain(best, rec)
#部分依赖图（PDP）
h2o.pd_plot(best,rec,"Fth1")

# 多模型的PDP图
h2o.pd_multi_plot(am@leaderboard, rec, "Fth1")

#个体条件期望图（ICE）
h2o.ice_plot(best,rec,show_pdp = TRUE,"Fth1")

#学习曲线
h2o.learning_curve_plot(best)


###验证
load('/Users/niuruize/Downloads/scRNA/BI/GSE197731/GSE197731_seurat.RData')
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) <- "Group"
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers)) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]
##
validation_expr = cbind(data.frame(t_expr),data.frame(Group=Idents(sce)))
rec <- recipe(Group~.,validation_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() # 转换为h2o 数据框
rec$Group <- h2o.asfactor(rec$Group)
y <- "Group"  # 因变量名        
x <- setdiff(names(rec),y)   # 预测变量名
pred <- h2o.predict(best, rec)
## 制作测试结果的混淆矩阵
Confusion_Matrix <- as.matrix(table(cbind(as.data.frame(pred$predict), as.data.frame(validation_expr$Group))))
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix_GSE197731.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')


##
###验证
#load('/Users/niuruize/Downloads/scRNA/BI/GSE174574/GSE174574_seurat.RData')
scRNA = scRNA[,scRNA$celltype %in% c("ExN")]
DefaultAssay(scRNA) <- "RNA"
Idents(scRNA) <- "Group"
sce <- scRNA
table(Idents(sce))
sce <- ScaleData(sce,features = unique(sce.markers$gene)) 
t_expr <- t(as.matrix(sce@assays$RNA@scale.data))
dim(t_expr) 
t_expr[1:4,1:4]
##
validation_expr = cbind(data.frame(t_expr),data.frame(Group=Idents(sce)))
rec <- recipe(Group~.,validation_expr) %>%           
  step_dummy(all_nominal_predictors()) %>%           
  prep() %>%           
  bake(new_data=NULL) %>%          
  as.h2o() # 转换为h2o 数据框
rec$Group <- h2o.asfactor(rec$Group)
y <- "Group"  # 因变量名        
x <- setdiff(names(rec),y)   # 预测变量名
pred <- h2o.predict(best, rec)
## 制作测试结果的混淆矩阵
Confusion_Matrix <- as.matrix(table(cbind(as.data.frame(pred$predict), as.data.frame(validation_expr$Group))))
Confusion_Matrix <- prop.table(Confusion_Matrix, margin = 1)
mycol<-colorRampPalette(c( "#104E8B", "white", "#8B0000"))(200)
p1 = pheatmap::pheatmap(Confusion_Matrix,scale = "none", border_color= "grey", number_color= "black",
                        fontsize_number=7,fontsize_row=10,fontsize_col=10,cellwidth=20,
                        cellheight=20,cluster_rows=F,cluster_cols=F,display_numbers = T,treeheight_row = 10,treeheight_col = 10,
                        color= mycol,show_rownames=T) 
ggsave(filename = "Normalized confusion matrix_Neuron.pdf", p1,device = 'pdf', width = 10, height = 10, units = 'cm')


#单个基因差异基因小提琴图
Idents(scRNA) <- "celltype"
table(scRNA$celltype)
scRNA$Group <- factor(scRNA$Group,levels=c("Adult","Aging"))
p1 <- VlnPlot(scRNA, features = c("MEG3"),pt.size = 0.01, split.by = "Group",ncol = 2)+theme(axis.title.x=element_text(size=0))
ggsave(filename = "MEG3_all.pdf", plot = p1, device = 'pdf', width = 20, height = 6)
Idents(scRNA) <- "Group"
p1 <- VlnPlot(scRNA, features = c("MEG3"),pt.size = 0.01,split.by = "Group",ncol = 1)+theme(axis.title.x=element_text(size=0))
ggsave(filename = "MEG3_all1.pdf", plot = p1, device = 'pdf', width = 5, height = 5)





