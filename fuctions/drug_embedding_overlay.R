#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
})

setwd("E:/Project")

source("https://raw.githubusercontent.com/YevhenAkimov/DatasetLT/main/DatasetLT.R")

config <- list(
  rna_rds = "E:/Project/seurat_obj_rna_QC.rds",
  atac_rds = "E:/Project/seurat_obj_atac_QC.rds",
  drug_rds = "E:/Project/results_v2/H180_H160_scRNA_to_phenotype_mapping/mapped_phenotypes/H160_r1_mapped_cGR_smoothed.rds",
  output_dir = file.path(getwd(), "seurat_plots"),
  overlay_dir = file.path(getwd(), "drug_overlays"),
  selections = list(
    seurat_assay = "RNA",
    embedding = "UMAP",
    drug_dataset = "cGR_mapped",
    drug_name = "adavosertib"
  ),
  dims = 30,
  clustering_resolution = 0.5,
  seed = 1234
)

stopifnot(file.exists(config$rna_rds), file.exists(config$atac_rds), file.exists(config$drug_rds))

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$overlay_dir, showWarnings = FALSE, recursive = TRUE)

preprocess_seurat <- function(obj, dims, resolution, seed, label) {
  DefaultAssay(obj) <- DefaultAssay(obj)
  
  keep_cells <- grep("-1$", colnames(obj), value = TRUE)
  if (!length(keep_cells)) stop("未找到以 -1 结尾的细胞条形码。")
  if (length(keep_cells) < ncol(obj)) {
    obj <- subset(obj, cells = keep_cells)
  }

  set.seed(seed)
  obj <- RunPCA(obj, npcs = dims)
  obj <- FindNeighbors(obj, dims = 1:dims)
  obj <- FindClusters(obj, resolution = resolution)
  obj <- RunUMAP(obj, dims = 1:dims)
  obj <- RunTSNE(obj, dims = 1:dims)

  message(sprintf("%s 对象预处理完成。", label))
  obj
}

cluster_matrix <- function(seurat_obj) {
  clusters <- Idents(seurat_obj)
  levels_vec <- levels(clusters)
  mat <- matrix(as.numeric(clusters), ncol = 1,
                dimnames = list(names(clusters), "cluster_id"))
  attr(mat, "cluster_levels") <- levels_vec
  mat
}

save_reduction_plot <- function(seurat_obj, label, reduction, output_dir) {
  if (!reduction %in% names(seurat_obj@reductions)) {
    warning(sprintf("%s 缺少 %s 降维，跳过绘图。", label, toupper(reduction)))
    return(invisible(NULL))
  }
  plt <- DimPlot(seurat_obj, reduction = reduction, label = TRUE) +
    ggtitle(sprintf("%s - %s", label, toupper(reduction))) +
    theme_classic()
  file_path <- file.path(output_dir, sprintf("%s_%s.png", label, toupper(reduction)))
  ggsave(file_path, plt, width = 6, height = 5, dpi = 300)
  message(sprintf("保存 %s 图：%s", toupper(reduction), file_path))
}

add_embeddings <- function(dslt, assay_name, seurat_obj, alias_map) {
  for (alias in names(alias_map)) {
    reduction <- alias_map[[alias]]
    if (!reduction %in% names(seurat_obj@reductions)) {
      next
    }
    emb <- Embeddings(seurat_obj, reduction = reduction)
    if (ncol(emb) < 2) {
      next
    }
    emb <- emb[, 1:2, drop = FALSE]
    colnames(emb) <- paste0(alias, "_", seq_len(ncol(emb)))
    dslt$addEmbedding(assay = assay_name, names = alias, input = emb)
  }
}

plot_drug_overlay <- function(
    dslt, assay_name, embedding_name, cluster_levels,
    drug_dataset, drug_name, output_path,
    z_cap = 3,                # 截断 z-score 的范围
    bg_pt_size = 0.6,         # 背景点大小（灰）
    cl_pt_size = 0.8,         # cluster 彩色点大小
    drug_pt_size = 1.6,       # 药物强度点大小（中空圆）
    drug_pt_stroke = 0.25     # 药物点描边
){
  library(ggplot2)
  library(scales)
  
  dslt$clearFilters()
  dslt$setActiveAssays(c(assay_name, drug_dataset))
  
  embedding <- dslt$getEmbedding(assay = assay_name, name = embedding_name)
  clusters  <- dslt$getAssay(assay_name)
  drug_mat  <- dslt$getAssay(drug_dataset)
  
  common <- Reduce(intersect, list(rownames(embedding), rownames(clusters), rownames(drug_mat)))
  embedding <- embedding[common, , drop = FALSE]
  clusters  <- clusters[common, , drop = FALSE]
  drug_vals <- as.numeric(drug_mat[common, drug_name])
  
  # z-score + 截断，避免极端值压扁色带
  mu <- mean(drug_vals, na.rm = TRUE)
  sdv <- sd(drug_vals, na.rm = TRUE)
  drug_z <- if (is.na(sdv) || sdv == 0) rep(0, length(drug_vals)) else (drug_vals - mu) / sdv
  drug_z <- pmax(pmin(drug_z, z_cap), -z_cap)
  
  # cluster 因子
  cl_id <- as.integer(round(clusters[, "cluster_id"]))
  cl_lbls <- cluster_levels[[assay_name]]
  lab <- rep(NA_character_, length(cl_id))
  ok  <- !is.na(cl_id) & cl_id >= 1 & cl_id <= length(cl_lbls)
  lab[ok] <- cl_lbls[cl_id[ok]]
  cl_fac <- factor(lab, levels = cl_lbls)
  
  df <- data.frame(embedding, cluster = cl_fac, drug_value = drug_z)
  x_col <- colnames(embedding)[1]; y_col <- colnames(embedding)[2]
  
  # Okabe–Ito 调色板（7 色，与你图中 0–6 个簇匹配；簇更多时会循环）
  okabe_ito <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
  
  plt <- ggplot(df, aes_string(x = x_col, y = y_col)) +
    
    # 1) 背景灰点（让整体轮廓更清楚）
    geom_point(color = "#D9D9D9", size = bg_pt_size, alpha = 0.8, stroke = 0) +
    
    geom_point(
      aes(fill = drug_value, color = cluster),
      shape = 21,
      size  = drug_pt_size,
      alpha = 0.95,
      stroke = drug_pt_stroke   # 轮廓线粗细，>0 才能看到簇颜色
    ) +
    
    # cluster 离散色
    scale_color_manual(values = rep(okabe_ito, length.out = nlevels(df$cluster)),
                       na.translate = FALSE, guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    
    # 药物强度连续色（对称发散、等距刻度；极值用 z_cap）
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      limits = c(-z_cap, z_cap), oob = squish,
      breaks = pretty(c(-z_cap, z_cap), n = 5),
      labels = number_format(accuracy = 0.1),
      midpoint = 0, name = sprintf("%s z-score", drug_name)
    ) +
    
    coord_equal() +
    labs(
      title = sprintf("%s — %s overlay %s", assay_name, embedding_name, drug_name),
      subtitle = "Grey: all cells • Color: clusters • Fill: drug z-score (clamped to ±3)",
      x = x_col, y = y_col
    ) +
    
    theme_classic(base_size = 12) +
    theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 10, margin = margin(b = 6)),
      axis.title    = element_text(size = 11),
      axis.text     = element_text(size = 9),
      legend.title  = element_text(face = "bold"),
      legend.text   = element_text(size = 9),
      legend.box    = "vertical",
      legend.spacing.y = unit(2, "pt"),
      legend.key.height = unit(10, "pt"),
      legend.key.width  = unit(12, "pt"),
      plot.margin  = margin(8, 8, 8, 8)
    )
  
  ggsave(output_path, plt, width = 7.2, height = 6.0, dpi = 320)
  message(sprintf("叠加图已保存：%s", output_path))
  plt
}


rna_obj <- readRDS(config$rna_rds)
# atac_obj <- readRDS(config$atac_rds)
# stopifnot(inherits(rna_obj, "Seurat"), inherits(atac_obj, "Seurat"))

drug_list <- readRDS(config$drug_rds)
stopifnot(is.list(drug_list), length(drug_list) > 0)
drug_list <- lapply(drug_list, function(mat) {
  mat <- as.matrix(mat)
  stopifnot(is.numeric(mat))
  mat
})

rna_obj <- preprocess_seurat(rna_obj, config$dims, config$clustering_resolution, config$seed, "RNA")
# atac_obj <- preprocess_seurat(atac_obj, config$dims, config$clustering_resolution, config$seed, "ATAC")

save_reduction_plot(rna_obj, "RNA", "umap", config$output_dir)
save_reduction_plot(rna_obj, "RNA", "tsne", config$output_dir)
# save_reduction_plot(atac_obj, "ATAC", "umap", config$output_dir)
# save_reduction_plot(atac_obj, "ATAC", "tsne", config$output_dir)

dslt <- DatasetLT$new()
cluster_levels <- list()

rna_clusters <- cluster_matrix(rna_obj)
cluster_levels[["RNA"]] <- attr(rna_clusters, "cluster_levels")
dslt$addAssay(names = "RNA", input = rna_clusters)
add_embeddings(dslt, "RNA", rna_obj, c(UMAP = "umap", TSNE = "tsne"))

# atac_clusters <- cluster_matrix(atac_obj)
# cluster_levels[["ATAC"]] <- attr(atac_clusters, "cluster_levels")
# dslt$addAssay(names = "ATAC", input = atac_clusters)
# add_embeddings(dslt, "ATAC", atac_obj, c(UMAP = "umap", TSNE = "tsne"))

for (assay_name in names(drug_list)) {
  dslt$addAssay(names = assay_name, input = drug_list[[assay_name]])
}

sel <- config$selections
stopifnot(sel$seurat_assay %in% c("RNA", "ATAC"))
stopifnot(sel$drug_dataset %in% names(drug_list))
stopifnot(sel$drug_name %in% colnames(drug_list[[sel$drug_dataset]]))
stopifnot(sel$embedding %in% c("UMAP", "TSNE"))

output_path <- file.path(
  config$overlay_dir,
  sprintf("%s_%s_%s_%s.png",
          sel$seurat_assay,
          sel$embedding,
          sel$drug_dataset,
          gsub("[^[:alnum:]]+", "_", sel$drug_name))
)

plot_drug_overlay(dslt,
                  assay_name = sel$seurat_assay,
                  embedding_name = sel$embedding,
                  cluster_levels = cluster_levels,
                  drug_dataset = sel$drug_dataset,
                  drug_name = sel$drug_name,
                  output_path = output_path)


saveRDS(dslt,"dslt.rds")
message("脚本执行完成。")