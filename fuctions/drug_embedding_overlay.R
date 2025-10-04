#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
})

resolve_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_index <- grep(file_arg, cmd_args)
  if (length(file_index) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", cmd_args[file_index]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

script_dir <- resolve_script_dir()
source(file.path(script_dir, "DatasetLT.R"))

config <- list(
  rna_rds = "path/to/rna_seurat.rds",
  atac_rds = "path/to/atac_seurat.rds",
  drug_rds = "path/to/drug_response_list.rds",
  output_dir = file.path(getwd(), "seurat_plots"),
  overlay_dir = file.path(getwd(), "drug_overlays"),
  selections = list(
    seurat_assay = "RNA",
    embedding = "UMAP",
    drug_dataset = "cGR",
    drug_name = "DrugA"
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
  set.seed(seed)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
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

plot_drug_overlay <- function(dslt, assay_name, embedding_name, cluster_levels, drug_dataset, drug_name, output_path) {
  dslt$clearFilters()
  dslt$setActiveAssays(c(assay_name, drug_dataset))
  embedding <- dslt$getEmbedding(assay = assay_name, name = embedding_name)
  clusters <- dslt$getAssay(assay_name)
  drug_matrix <- dslt$getAssay(drug_dataset)
  common <- Reduce(intersect, list(rownames(embedding), rownames(clusters), rownames(drug_matrix)))
  embedding <- embedding[common, , drop = FALSE]
  clusters <- clusters[common, , drop = FALSE]
  drug_values <- drug_matrix[common, drug_name]
  cluster_ids <- as.integer(round(clusters[, "cluster_id"]))
  cluster_labels <- cluster_levels[[assay_name]]
  label_vec <- rep(NA_character_, length(cluster_ids))
  valid_idx <- !is.na(cluster_ids) & cluster_ids >= 1 & cluster_ids <= length(cluster_labels)
  label_vec[valid_idx] <- cluster_labels[cluster_ids[valid_idx]]
  cluster_factor <- factor(label_vec, levels = cluster_labels)
  df <- data.frame(embedding, cluster = cluster_factor, drug_value = as.numeric(drug_values))
  x_col <- colnames(embedding)[1]
  y_col <- colnames(embedding)[2]
  plt <- ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = cluster), size = 1, alpha = 0.5) +
    geom_point(aes(fill = drug_value), shape = 21, color = "black", size = 1.5, alpha = 0.85) +
    scale_fill_viridis_c(option = "magma", direction = -1, name = sprintf("%s 值", drug_name)) +
    guides(color = guide_legend(title = "Cluster")) +
    coord_equal() +
    theme_classic() +
    labs(title = sprintf("%s - %s overlay %s", assay_name, embedding_name, drug_name),
         x = x_col, y = y_col)
  ggsave(output_path, plt, width = 6, height = 5, dpi = 300)
  message(sprintf("叠加图已保存：%s", output_path))
  plt
}

rna_obj <- readRDS(config$rna_rds)
atac_obj <- readRDS(config$atac_rds)
stopifnot(inherits(rna_obj, "Seurat"), inherits(atac_obj, "Seurat"))

drug_list <- readRDS(config$drug_rds)
stopifnot(is.list(drug_list), length(drug_list) > 0)
drug_list <- lapply(drug_list, function(mat) {
  mat <- as.matrix(mat)
  stopifnot(is.numeric(mat))
  mat
})

rna_obj <- preprocess_seurat(rna_obj, config$dims, config$clustering_resolution, config$seed, "RNA")
atac_obj <- preprocess_seurat(atac_obj, config$dims, config$clustering_resolution, config$seed, "ATAC")

save_reduction_plot(rna_obj, "RNA", "umap", config$output_dir)
save_reduction_plot(rna_obj, "RNA", "tsne", config$output_dir)
save_reduction_plot(atac_obj, "ATAC", "umap", config$output_dir)
save_reduction_plot(atac_obj, "ATAC", "tsne", config$output_dir)

dslt <- DatasetLT$new()
cluster_levels <- list()

rna_clusters <- cluster_matrix(rna_obj)
cluster_levels[["RNA"]] <- attr(rna_clusters, "cluster_levels")
dslt$addAssay(names = "RNA", input = rna_clusters)
add_embeddings(dslt, "RNA", rna_obj, c(UMAP = "umap", TSNE = "tsne"))

atac_clusters <- cluster_matrix(atac_obj)
cluster_levels[["ATAC"]] <- attr(atac_clusters, "cluster_levels")
dslt$addAssay(names = "ATAC", input = atac_clusters)
add_embeddings(dslt, "ATAC", atac_obj, c(UMAP = "umap", TSNE = "tsne"))

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

message("脚本执行完成。")
