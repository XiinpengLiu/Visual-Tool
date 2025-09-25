# 安装一次即可
# install.packages("tidyseurat")   # CRAN
library(tidyseurat)
library(dplyr)
library(tidyr)
library(ggplot2)
se_rna <- RunUMAP(se_rna, dims = 1:30)
# 先看你对象里有哪些降维结果
names(Reductions(se_rna))

# —— 更稳的版本：构建 {UMAP,PCA,TSNE} 矩阵列表 —— 
get_embed_mat <- function(se, red = c("umap","pca","tsne"), cluster_col = NULL) {
  # 取一个 reduction 的 x,y + cluster_id 矩阵；抓不到就返回 NULL
  grab_one <- function(reduction) {
    # 1) 优先用 Embeddings(); 失败再用 se[[reduction]]@cell.embeddings
    emb <- tryCatch(Embeddings(se, reduction = reduction), error = function(e) NULL)
    if (is.null(emb)) {
      emb <- tryCatch(se[[reduction]]@cell.embeddings, error = function(e) NULL)
    }
    if (is.null(emb)) return(NULL)
    emb <- emb[, 1:2, drop = FALSE]
    colnames(emb) <- c("x","y")
    
    # 2) 聚类号：优先指定 meta 列，否则用 Idents；按 rownames(emb) 对齐
    cl <- if (!is.null(cluster_col) && cluster_col %in% colnames(se@meta.data)) {
      se@meta.data[rownames(emb), cluster_col, drop = TRUE]
    } else {
      Idents(se)[rownames(emb)]
    }
    cl_fac <- factor(as.character(cl))
    mat <- cbind(emb, cluster_id = as.integer(cl_fac))
    rownames(mat) <- rownames(emb)
    attr(mat, "cluster_levels") <- levels(cl_fac)  # 保存映射
    mat
  }
  
  # 逐个抓，统一大写命名，自动剔除 NULL
  mats <- setNames(lapply(tolower(red), grab_one), toupper(tolower(red)))
  Filter(Negate(is.null), mats)
}

# —— 用法 —— 
em_lists <- get_embed_mat(se_rna, red = c("umap","pca","tsne"), cluster_col = NULL)

df_from_mat <- function(mat) {
  labs <- attr(mat, "cluster_levels")
  df <- as.data.frame(mat)
  df$cluster <- if (!is.null(labs)) factor(labs[df$cluster_id], levels = labs)
  else factor(df$cluster_id)
  df
}

# 1) 单个降维：按cluster上色
plot_clusters <- function(em_lists, red = "UMAP", point_size = 0.25) {
  stopifnot(red %in% names(em_lists))
  df <- df_from_mat(em_lists[[red]])
  ggplot(df, aes(x, y, color = cluster)) +
    geom_point(size = point_size, alpha = 0.8) +
    coord_equal() + theme_classic() +
    labs(title = paste0(red, " by cluster"),
         x = paste0(red, "_1"), y = paste0(red, "_2"), color = "cluster")
}

p1 <- plot_clusters(em_lists, red = "UMAP")
p2 <- plot_clusters(em_lists, red = "PCA")
p1
p2
