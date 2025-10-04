#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
})

# Locate DatasetLT implementation shipped with the repo and source it -----------------------------
resolve_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  file_index <- grep(file_arg, cmd_args)
  if (length(file_index) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", cmd_args[file_index]))))
  }
  # Fallback for interactive sourcing
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  getwd()
}

script_dir <- resolve_script_dir()
source(file.path(script_dir, "DatasetLT.R"))

# -----------------------------------------------------------------------------------------------
# Helper utilities ------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
ask_path <- function(prompt) {
  repeat {
    path <- readline(prompt = prompt)
    if (nzchar(path) && file.exists(path)) {
      return(normalizePath(path))
    }
    message("路径无效，请重新输入。")
  }
}

choose_option <- function(options, title) {
  if (length(options) == 0) {
    stop(paste("没有可供选择的", title))
  }
  choice <- utils::menu(options, title = title)
  if (choice == 0) {
    stop("用户取消选择，脚本终止。")
  }
  options[[choice]]
}

load_seurat_object <- function(path) {
  obj <- tryCatch(readRDS(path), error = function(e) {
    stop(sprintf("无法读取 Seurat 对象：%s", e$message))
  })
  if (!inherits(obj, "Seurat")) {
    stop(sprintf("文件 %s 不是 Seurat 对象。", path))
  }
  obj
}

find_reduction_name <- function(seurat_obj, target) {
  red_names <- names(seurat_obj@reductions)
  if (length(red_names) == 0) {
    return(NULL)
  }
  target_lower <- tolower(target)
  exact <- red_names[tolower(red_names) == target_lower]
  if (length(exact) > 0) {
    return(exact[[1]])
  }
  partial <- red_names[grepl(target_lower, red_names, ignore.case = TRUE)]
  if (length(partial) > 0) {
    return(partial[[1]])
  }
  NULL
}

extract_embedding_matrix <- function(seurat_obj, reduction_name, label) {
  emb <- tryCatch(Embeddings(seurat_obj, reduction = reduction_name),
                  error = function(e) NULL)
  if (is.null(emb)) {
    return(NULL)
  }
  emb <- emb[, 1:2, drop = FALSE]
  if (ncol(emb) < 2) {
    return(NULL)
  }
  colnames(emb) <- paste0(label, "_", seq_len(ncol(emb)))
  emb
}

cluster_matrix_from_seurat <- function(seurat_obj) {
  clusters <- Idents(seurat_obj)
  cluster_levels <- levels(clusters)
  matrix_data <- matrix(as.numeric(clusters), ncol = 1,
                        dimnames = list(names(clusters), "cluster_id"))
  attr(matrix_data, "cluster_levels") <- cluster_levels
  matrix_data
}

save_dimensional_plots <- function(seurat_obj, object_label, output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  reductions_to_plot <- c("umap", "tsne")
  generated <- list()
  for (red in reductions_to_plot) {
    red_name <- find_reduction_name(seurat_obj, red)
    if (is.null(red_name)) {
      message(sprintf("%s 对象未找到 %s 降维结果，跳过绘图。", object_label, toupper(red)))
      next
    }
    plt <- tryCatch({
      DimPlot(seurat_obj, reduction = red_name, label = TRUE) +
        ggtitle(sprintf("%s - %s", object_label, toupper(red))) +
        theme_classic()
    }, error = function(e) {
      message(sprintf("绘制 %s %s 失败：%s", object_label, toupper(red), e$message))
      NULL
    })
    if (!is.null(plt)) {
      file_path <- file.path(output_dir, sprintf("%s_%s.png", object_label, toupper(red)))
      ggsave(file_path, plt, width = 6, height = 5, dpi = 300)
      generated[[toupper(red)]] <- file_path
      message(sprintf("已保存 %s 图：%s", toupper(red), file_path))
    }
  }
  generated
}

add_embeddings_to_dataset <- function(dslt, assay_name, seurat_obj, cluster_info, embedding_alias) {
  dslt$addAssay(names = assay_name, input = cluster_info)
  available <- character()
  for (alias in names(embedding_alias)) {
    emb <- extract_embedding_matrix(seurat_obj, embedding_alias[[alias]], alias)
    if (is.null(emb)) {
      next
    }
    dslt$addEmbedding(assay = assay_name, names = alias, input = emb)
    available <- c(available, alias)
  }
  available
}

plot_drug_overlay <- function(dslt, assay_name, embedding_name, cluster_levels, drug_assay, drug_name, output_path = NULL) {
  dslt$clearFilters()
  dslt$setActiveAssays(c(assay_name, drug_assay))

  embedding <- dslt$getEmbedding(assay = assay_name, name = embedding_name)
  assay_data <- dslt$getAssay(assay_name)
  drug_matrix <- dslt$getAssay(drug_assay)

  common_cells <- Reduce(intersect, list(rownames(embedding), rownames(assay_data), rownames(drug_matrix)))
  embedding <- embedding[common_cells, , drop = FALSE]
  assay_data <- assay_data[common_cells, , drop = FALSE]
  drug_values <- drug_matrix[common_cells, drug_name]

  cluster_ids <- as.integer(round(assay_data[, "cluster_id"]))
  cluster_levels_vec <- cluster_levels[[assay_name]]
  cluster_factor <- if (!is.null(cluster_levels_vec)) {
    labels <- rep(NA_character_, length(cluster_ids))
    valid_idx <- !is.na(cluster_ids) & cluster_ids >= 1 & cluster_ids <= length(cluster_levels_vec)
    labels[valid_idx] <- cluster_levels_vec[cluster_ids[valid_idx]]
    factor(labels, levels = cluster_levels_vec)
  } else {
    factor(cluster_ids)
  }

  df <- data.frame(embedding, cluster = cluster_factor, drug_value = as.numeric(drug_values))
  x_col <- colnames(embedding)[1]
  y_col <- colnames(embedding)[2]

  plt <- ggplot(df, aes_string(x = x_col, y = y_col)) +
    geom_point(aes(color = cluster), size = 1, alpha = 0.4) +
    geom_point(aes(fill = drug_value), shape = 21, color = "black", size = 1.6, alpha = 0.85) +
    scale_fill_viridis_c(option = "magma", direction = -1, name = sprintf("%s 值", drug_name)) +
    guides(color = guide_legend(title = "Cluster")) +
    coord_equal() +
    theme_classic() +
    labs(title = sprintf("%s - %s overlay %s", assay_name, embedding_name, drug_name),
         x = x_col, y = y_col)

  if (!is.null(output_path)) {
    ggsave(output_path, plt, width = 6, height = 5, dpi = 300)
    message(sprintf("叠加图已保存：%s", output_path))
  }
  plt
}

# -----------------------------------------------------------------------------------------------
# Main workflow ---------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------
message("请提供 QC 过的 RNA Seurat 对象 (.rds) 路径：")
rna_path <- ask_path("RNA Seurat 对象路径: ")
rna_obj <- load_seurat_object(rna_path)

message("请提供 QC 过的 ATAC Seurat 对象 (.rds) 路径：")
atac_path <- ask_path("ATAC Seurat 对象路径: ")
atac_obj <- load_seurat_object(atac_path)

output_dir <- file.path(getwd(), "seurat_plots")
message(sprintf("降维图将输出到目录：%s", output_dir))

rna_plots <- save_dimensional_plots(rna_obj, "RNA", output_dir)
atac_plots <- save_dimensional_plots(atac_obj, "ATAC", output_dir)

message("请提供药物反应 RDS 文件路径 (列表对象)：")
drug_rds_path <- ask_path("药物反应 RDS 路径: ")
drug_list <- tryCatch(readRDS(drug_rds_path), error = function(e) {
  stop(sprintf("无法读取药物反应数据：%s", e$message))
})
if (!is.list(drug_list) || length(drug_list) == 0) {
  stop("药物反应 RDS 必须是一个非空列表，每个元素为矩阵。")
}

if (is.null(names(drug_list)) || any(names(drug_list) == "")) {
  names(drug_list) <- paste0("assay_", seq_along(drug_list))
}

# Validate matrices
for (nm in names(drug_list)) {
  mat <- as.matrix(drug_list[[nm]])
  if (!is.numeric(mat)) {
    stop(sprintf("药物反应数据 %s 不是数值矩阵。", nm))
  }
  drug_list[[nm]] <- mat
}

dslt <- DatasetLT$new()
cluster_levels <- list()
embedding_records <- list()

rna_clusters <- cluster_matrix_from_seurat(rna_obj)
cluster_levels[["RNA"]] <- attr(rna_clusters, "cluster_levels")
emb_alias_rna <- list(UMAP = find_reduction_name(rna_obj, "umap"),
                      TSNE = find_reduction_name(rna_obj, "tsne"))
emb_alias_rna <- Filter(Negate(is.null), emb_alias_rna)
if (length(emb_alias_rna) > 0) {
  embedding_records[["RNA"]] <- add_embeddings_to_dataset(dslt, "RNA", rna_obj, rna_clusters, emb_alias_rna)
} else {
  dslt$addAssay(names = "RNA", input = rna_clusters)
  embedding_records[["RNA"]] <- character(0)
}

atac_clusters <- cluster_matrix_from_seurat(atac_obj)
cluster_levels[["ATAC"]] <- attr(atac_clusters, "cluster_levels")
emb_alias_atac <- list(UMAP = find_reduction_name(atac_obj, "umap"),
                       TSNE = find_reduction_name(atac_obj, "tsne"))
emb_alias_atac <- Filter(Negate(is.null), emb_alias_atac)
if (length(emb_alias_atac) > 0) {
  embedding_records[["ATAC"]] <- add_embeddings_to_dataset(dslt, "ATAC", atac_obj, atac_clusters, emb_alias_atac)
} else {
  dslt$addAssay(names = "ATAC", input = atac_clusters)
  embedding_records[["ATAC"]] <- character(0)
}

for (assay_name in names(drug_list)) {
  dslt$addAssay(names = assay_name, input = drug_list[[assay_name]])
}

message("成功将 Seurat 嵌入与药物数据载入 DatasetLT。")

# 用户选择药物与数据类型 ---------------------------------------------------------------
drug_dataset <- choose_option(names(drug_list), "请选择药物数据类型：")
drug_names <- colnames(drug_list[[drug_dataset]])
if (is.null(drug_names) || length(drug_names) == 0) {
  stop(sprintf("药物数据 %s 缺少列名，无法继续。", drug_dataset))
}
drug_name <- choose_option(drug_names, sprintf("请选择 %s 中的药物：", drug_dataset))

seurat_assay_choice <- choose_option(names(embedding_records), "请选择要使用的 Seurat 对象 (RNA/ATAC)：")
available_emb <- embedding_records[[seurat_assay_choice]]
if (length(available_emb) == 0) {
  stop(sprintf("选择的 %s 对象没有可用的嵌入。", seurat_assay_choice))
}
embedding_choice <- choose_option(available_emb, sprintf("请选择 %s 对象的嵌入：", seurat_assay_choice))

overlay_dir <- file.path(getwd(), "drug_overlays")
dir.create(overlay_dir, showWarnings = FALSE, recursive = TRUE)
safe_drug_name <- gsub("[^[:alnum:]]+", "_", drug_name)
overlay_path <- file.path(overlay_dir, sprintf("%s_%s_%s_%s.png",
                                              seurat_assay_choice,
                                              embedding_choice,
                                              drug_dataset,
                                              safe_drug_name))

overlay_plot <- plot_drug_overlay(dslt, seurat_assay_choice, embedding_choice,
                                  cluster_levels, drug_dataset, drug_name,
                                  output_path = overlay_path)

print(overlay_plot)

message("脚本执行完成。")
