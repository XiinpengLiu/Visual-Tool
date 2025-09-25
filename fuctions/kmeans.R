se_rna <- NormalizeData(se_rna)
se_rna <- FindVariableFeatures(se_rna)
se_rna <- ScaleData(se_rna)
se_rna <- RunPCA(se_rna)

# 假设 se_rna 已经做了 PCA
pc <- Embeddings(se_rna, "pca")[, 1:50]  # 用前 10 个主成分
km <- kmeans(pc, centers = 5, nstart = 20)
se_rna$kmeans5 <- as.factor(km$cluster)

# 用 Seurat 的 DimPlot
DimPlot(se_rna, reduction = "pca", group.by = "kmeans5")

# 或者自己用 ggplot2
library(ggplot2)
pc_df <- as.data.frame(pc)
pc_df$cluster <- as.factor(km$cluster)
ggplot(pc_df, aes(PC_1, PC_2, color = cluster)) +
  geom_point() +
  theme_classic()
