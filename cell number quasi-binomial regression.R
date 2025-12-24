# 加载dplyr包
library(dplyr)

# 提取metadata
seurat_objet <- scRNA[,scRNA$Group %in% c("Adult","Aging")]
metadata <- seurat_objet@meta.data

# 查看metadata的列名，确保列名正确
head(metadata)

# 提取所需列
data <- metadata[, c("SampleID", "Group", "celltype", "Region")]


# 计算每个样本中每个细胞亚群的比例
celltype_proportions <- data %>%
  group_by(SampleID, Group, Region, celltype) %>%  # 按样本、疾病状态、脑区和细胞亚群分组
  summarise(count = n()) %>%                       # 计算每个亚群的细胞数
  mutate(proportion = count / sum(count))          # 计算比例

# 查看结果
head(celltype_proportions)

# 计算每个样本的总细胞数
total_cells <- data %>%
  group_by(SampleID, Region) %>%
  summarise(total_cells = n())

# 合并比例和总细胞数
celltype_proportions <- celltype_proportions %>%
  left_join(total_cells, by = c("SampleID", "Region"))

# 查看最终数据
head(celltype_proportions)

library(writexl)

# 导出为 Excel 文件
write_xlsx(celltype_proportions, "celltype_proportions.xlsx")


# 读取Excel文件
data <- read_excel("celltype_proportions.xlsx", sheet = "Sheet1")

# 确保Group和celltype为因子
df=data
df$Group <- factor(df$Group)
df$celltype <- factor(df$celltype)
df$Region <- factor(df$Region)
table(scRNA$celltype,scRNA$Region)

# 创建一个空的数据框来存储结果
results <- data.frame()
celltype <- c("Astro1","Astro2","Astro4","Astro5","Astro6","Astro7");Region <- c("ACC","AMY","EC","HIP","MB","PFC","SC")
# 对每个细胞亚群和每个脑区分别拟合模型
for (cell in celltype) {
  for (region in Region) {
    
    # 针对每个细胞亚群和脑区，子集化数据
    region_cell_data <- subset(df, celltype == cell & Region == region)
    
    # 拟合quasibinomial回归模型，计算成年和老年之间的差异
    model <- glm(proportion ~ Group, 
                 family = quasibinomial(), 
                 data = region_cell_data)
    
    # 获取模型摘要
    summary_model <- summary(model)
    
    # 提取回归系数和P值
    coefficients <- summary_model$coefficients
    # 计算log2(系数)
    log2_coefficients <- log2(abs(coefficients[, "Estimate"]))
    p_values <- coefficients[, "Pr(>|t|)"]  # 提取P值
    
    # 创建一个包含细胞亚群、脑区、log2(系数)和P值的数据框
    region_cell_results <- data.frame(
      Region = region,
      celltype = cell,
      log2_Coefficient = log2_coefficients,
      P_Value = p_values,
      Coefficient = coefficients[, "Estimate"]
    )
    
    # 计算成年组和老年组中每个细胞亚群的比例
    df_grouped <- region_cell_data %>%
      group_by(Group, celltype) %>%
      summarise(mean_proportion = mean(count)) %>%
      spread(Group, mean_proportion)
    
    # 查看成年组与老年组的比例变化
    head(df_grouped)
    
    # 计算老年组与成年组的比例差异
    df_grouped <- df_grouped %>%
      mutate(proportion_diff = `Aging` - `Adult`)  # 老年组与成年组的比例差异
    
    # 查看比例差异
    head(df_grouped)
    
    # 根据变化方向调整log2(系数) 即衰老组与成年组的差值
    region_cell_results$Adjusted_log2_Coefficient <- df_grouped$proportion_diff
    #region_cell_results$log2_Coefficient <- region_cell_results$log2_Coefficient * region_cell_results$Adjusted_log2_Coefficient
    
    # 将每个细胞亚群和脑区的结果追加到总结果数据框中
    results <- rbind(results, region_cell_results)
  }
}

# 查看结果
print(results)


# FDR修正
results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "fdr")

# 查看FDR修正后的结果
print(results)

# 创建log2(系数)矩阵
heatmap_matrix <- reshape(results, 
                          idvar = c("celltype", "Region"), 
                          timevar = "celltype", 
                          direction = "wide")

# 将P值也重塑为矩阵
p_value_matrix <- reshape(results, 
                          idvar = c("celltype", "Region"), 
                          timevar = "celltype", 
                          direction = "wide", 
                          v.names = "Adjusted_P_Value")

# 使用ggplot2绘制热图
library(ggplot2)
# 标准化log2_Coefficient到-1到1范围
min_val <- min(results$log2_Coefficient, na.rm = TRUE)  # 获取log2(系数)的最小值
max_val <- max(results$log2_Coefficient, na.rm = TRUE)  # 获取log2(系数)的最大值
results$Normalized_log2_Coefficient <- 
  (results$log2_Coefficient - min_val) / (max_val - min_val) * 2 - 1
# 标准化Adjusted_log2_Coefficient到-100到100范围
# 限制Adjusted_log2_Coefficient的范围
results$Adjusted_log2_Coefficient <- ifelse(results$Adjusted_log2_Coefficient > 100, 100, 
                                            ifelse(results$Adjusted_log2_Coefficient < -100, -100, 
                                                   results$Adjusted_log2_Coefficient))

# 查看处理后的结果
head(results)

# 查看处理后的结果
head(results)

# 使用ggplot2绘制标准化后的log2(Coefficient)热图
library(ggplot2)

ggplot(results, aes(x = celltype, y = Region, fill = Adjusted_log2_Coefficient)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  labs(title = "Relative enrichment of major cell types across regions",
       x = "Cell Type",
       y = "Region") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(aes(label = ifelse(Adjusted_P_Value < 0.05, "*", "")), color = "black", size = 5) # 显示星号

ggsave(filename = "Astro_subtype_number.pdf",device = 'pdf', width = 11, height = 7, units = 'cm')

