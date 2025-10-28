source("https://raw.githubusercontent.com/YevhenAkimov/graphics-R/main/graphics_functions.R")
dslt_local_path <- file.path("fuctions", "DatasetLT_modi.R")
if (file.exists(dslt_local_path)) {
  source(dslt_local_path)
} else {
  source("https://raw.githubusercontent.com/YevhenAkimov/Visual-Tool/main/fuctions/final/DatasetLT_modi.R")
}
source("https://raw.githubusercontent.com/YevhenAkimov/phenomics_scripts/main/phenomics_helpers.R")
source("https://raw.githubusercontent.com/YevhenAkimov/general_purpose_R/main/general_helpers.R")
source("https://raw.githubusercontent.com/YevhenAkimov/graphics-R/main/colors.R")

options(shiny.maxRequestSize = 10240 * 1024^2) # 1000 MB
options(future.globals.maxSize = 2 * 1024^3)  # 2 GB


# -------------------------------------------------------------------------
# Helper utilities --------------------------------------------------------
# -------------------------------------------------------------------------

read_table_like <- function(path) {
  stopifnot(length(path) == 1)
  ext <- tools::file_ext(path)
  switch(tolower(ext),
    rds = readRDS(path),
    rdata = {
      e <- new.env(parent = emptyenv())
      load(path, envir = e)
      e[[ls(e)[1]]]
    },
    csv = readr::read_csv(path, show_col_types = FALSE) %>% tibble::column_to_rownames(1),
    tsv = readr::read_tsv(path, show_col_types = FALSE) %>% tibble::column_to_rownames(1),
    txt = readr::read_tsv(path, show_col_types = FALSE) %>% tibble::column_to_rownames(1),
    stop(glue::glue("Unsupported file extension: {ext}"))
  )
}

read_mapping_file <- function(path) {
  obj <- read_table_like(path)
  if (is.data.frame(obj)) {
    obj <- as.data.frame(obj)
  }
  stopifnot(is.data.frame(obj))
  if (!"cell_barcode" %in% names(obj) && !is.null(rownames(obj))) {
    obj$cell_barcode <- rownames(obj)
  }
  original_names <- names(obj)
  lower_names <- tolower(original_names)
  locate_column <- function(target) {
    idx <- which(lower_names == tolower(target))
    if (length(idx)) {
      original_names[idx[1]]
    } else {
      NA_character_
    }
  }

  rename_map <- c(
    cell_barcode = locate_column("cell_barcode"),
    UMI = locate_column("umi"),
    Header = locate_column("header"),
    Barcode = locate_column("barcode")
  )

  mapping <- obj
  for (nm in names(rename_map)) {
    col_name <- rename_map[[nm]]
    if (!is.na(col_name)) {
      mapping[[nm]] <- mapping[[col_name]]
    }
  }

  required <- c("cell_barcode", "Barcode", "UMI", "Header")
  missing_required <- setdiff(required, names(mapping))
  if (length(missing_required)) {
    stop(paste("Mapping file missing required columns:", paste(missing_required, collapse = ", ")))
  }

  mapping$cell_barcode <- as.character(mapping$cell_barcode)
  mapping$Barcode <- as.character(mapping$Barcode)
  mapping$UMI <- as.character(mapping$UMI)
  mapping$Header <- as.character(mapping$Header)

  mapping
}

validate_10x_components <- function(file_names) {
  lower <- tolower(file_names)
  required_patterns <- list(
    matrix = "matrix",
    mtx = "mtx",
    barcodes = "barcode",
    features = "feature|gene|peak"
  )
  checks <- vapply(required_patterns, function(pattern) {
    any(grepl(pattern, lower))
  }, logical(1))
  if (!all(checks[c("matrix", "mtx", "barcodes", "features")])) {
    stop("Please upload matrix.mtx(.gz), features/peaks.tsv(.gz), and barcodes.tsv(.gz) files, or provide a single 10x HDF5 (.h5/.hdf5) file.")
  }
  invisible(TRUE)
}

copy_upload_to_temp_dir <- function(file_info) {
  if (is.null(file_info) || nrow(file_info) == 0) {
    return(NULL)
  }
  validate_10x_components(file_info$name)
  temp_dir <- tempfile("tenx_")
  dir.create(temp_dir)
  for (i in seq_len(nrow(file_info))) {
    dest <- file.path(temp_dir, basename(file_info$name[i]))
    file.copy(file_info$datapath[i], dest, overwrite = TRUE)
  }
  temp_dir
}

read_10x_from_upload <- function(file_info) {
  if (is.null(file_info) || nrow(file_info) == 0) {
    return(NULL)
  }

  lower_names <- tolower(file_info$name)
  if (nrow(file_info) == 1 && grepl("\\.(h5|hdf5)$", lower_names)) {
    temp_dir <- tempfile("tenx_h5_")
    dir.create(temp_dir)
    dest <- file.path(temp_dir, basename(file_info$name[1]))
    file.copy(file_info$datapath[1], dest, overwrite = TRUE)
    on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
    return(Read10X_h5(dest))
  }

  temp_dir <- copy_upload_to_temp_dir(file_info)
  if (is.null(temp_dir)) {
    return(NULL)
  }
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)
  Read10X(data.dir = temp_dir)
}

preserve_uploaded_file <- function(file_info) {
  if (is.null(file_info)) {
    return(NULL)
  }
  stopifnot(nrow(file_info) == 1)
  temp_dir <- tempfile("upload_")
  dir.create(temp_dir)
  dest <- file.path(temp_dir, basename(file_info$name))
  file.copy(file_info$datapath, dest, overwrite = TRUE)
  dest
}

preserve_fragments_with_index <- function(file_info) {
  if (is.null(file_info) || nrow(file_info) == 0) {
    return(list(fragments = NULL, index = NULL))
  }

  temp_dir <- tempfile("fragments_")
  dir.create(temp_dir)

  lower_names <- tolower(file_info$name)
  dests <- character(nrow(file_info))
  for (i in seq_len(nrow(file_info))) {
    dests[i] <- file.path(temp_dir, basename(file_info$name[i]))
    file.copy(file_info$datapath[i], dests[i], overwrite = TRUE)
  }

  fragment_idx <- which(grepl("\\.(tsv|tsv\\.gz)$", lower_names))
  if (!length(fragment_idx)) {
    return(list(fragments = NULL, index = NULL))
  }

  fragments_path <- dests[fragment_idx[1]]
  index_candidates <- dests[grepl("\\.tbi(\\.gz)?$", lower_names)]

  if (length(index_candidates) && grepl("\\.tbi\\.gz$", index_candidates[1])) {
    if (requireNamespace("R.utils", quietly = TRUE)) {
      dest_index <- sub("\\.gz$", "", index_candidates[1])
      tryCatch({
        R.utils::gunzip(index_candidates[1], destname = dest_index, overwrite = TRUE, remove = TRUE)
        index_candidates[1] <- dest_index
      }, error = function(e) {
        NULL
      })
    }
  }

  list(
    fragments = fragments_path,
    index = if (length(index_candidates)) index_candidates[1] else NULL
  )
}

ensure_tabix_index <- function(fragments_path) {
  if (is.null(fragments_path) || !file.exists(fragments_path)) {
    return(FALSE)
  }

  index_path <- paste0(fragments_path, ".tbi")
  if (file.exists(index_path)) {
    return(TRUE)
  }

  if (!requireNamespace("Rsamtools", quietly = TRUE)) {
    return(FALSE)
  }

  if (!grepl("\\.gz$", fragments_path)) {
    return(FALSE)
  }

  success <- tryCatch({
    Rsamtools::indexTabix(fragments_path, format = "bed")
    file.exists(index_path)
  }, error = function(e) {
    FALSE
  })

  success
}

normalize_matrix <- function(mat) {
  if (is.data.frame(mat)) mat <- as.matrix(mat)
  if (!is.matrix(mat)) {
    stop("Expected a matrix-like object")
  }
  if (is.null(rownames(mat))) {
    rownames(mat) <- make.unique(rep("sample", nrow(mat)))
  }
  if (is.null(colnames(mat))) {
    colnames(mat) <- paste0("feature_", seq_len(ncol(mat)))
  }
  mode(mat) <- "numeric"
  mat
}

ensure_numeric_vector <- function(x, default = NA_real_) {
  if (length(x) == 0) return(default)
  suppressWarnings(v <- as.numeric(x))
  if (anyNA(v)) default else v
}

make_heatmap_df <- function(mat, rows = NULL, cols = NULL) {
  if (is.null(mat)) return(NULL)
  mat <- normalize_matrix(mat)
  if (!is.null(rows)) mat <- mat[rows, , drop = FALSE]
  if (!is.null(cols)) mat <- mat[, cols, drop = FALSE]
  df <- as.data.frame(mat)
  df$barcode <- rownames(df)
  tidyr::pivot_longer(df, -barcode, names_to = "feature", values_to = "value")
}

combine_feature_plots <- function(result) {
  if (is.null(result) || !is.list(result)) return(NULL)
  violin <- result[["violin"]]
  feature <- result[["feature_plot"]]
  if (inherits(violin, "ggplot") && inherits(feature, "ggplot")) {
    violin + feature
  } else if (inherits(violin, "ggplot")) {
    violin
  } else if (inherits(feature, "ggplot")) {
    feature
  } else {
    NULL
  }
}

update_export_history <- function(state, action) {
  state$export_history <- rbind(
    state$export_history,
    data.frame(timestamp = Sys.time(), action = action, stringsAsFactors = FALSE)
  )
}

extract_drug_choices <- function(assay_list) {
  if (is.null(assay_list)) return(character())
  unique(unlist(lapply(assay_list, colnames)))
}

compose_upload_summary <- function(state) {
  lineage_status <- if (is.null(state$dslt)) "Not loaded" else "Loaded"
  sc_rna_status <- if (!is.null(state$seurat$pb_rna)) {
    "Pseudo-bulk ready"
  } else if (!is.null(state$seurat$sc_rna_raw)) {
    "Uploaded"
  } else {
    "Not uploaded"
  }
  sc_atac_status <- if (!is.null(state$seurat$pb_atac)) {
    "Pseudo-bulk ready"
  } else if (!is.null(state$seurat$sc_atac_raw)) {
    "Uploaded"
  } else {
    "Not uploaded"
  }
  mapping_status <- if (!is.null(state$mapping$rna) || !is.null(state$mapping$atac)) "Provided" else "Not provided"
  drug_matrix_status <- if (is.null(state$drug_matrix)) "Not provided" else "Provided"
  denoised_status <- if (length(state$denoised_assays)) paste(state$denoised_assays, collapse = ", ") else "None"

  ready <- c()
  if (!is.null(state$seurat$pb_rna)) ready <- c(ready, "RNA")
  if (!is.null(state$seurat$pb_atac)) ready <- c(ready, "ATAC")

  lines <- c(
    paste("Lineage data:", lineage_status),
    paste("Single-cell RNA:", sc_rna_status),
    paste("Single-cell ATAC:", sc_atac_status),
    paste("Mapping files:", mapping_status),
    paste("External drug matrix:", drug_matrix_status),
    paste("Denoised assays:", denoised_status)
  )

  if (length(ready)) {
    lines <- c(lines, paste("Pseudo-bulk generated for:", paste(ready, collapse = ", ")))
  }

  paste(lines, collapse = "\n")
}

get_drug_scores <- function(assays, assay_name, drugs) {
  if (length(drugs) == 0 || is.null(assays)) return(NULL)
  mat <- assays[[assay_name]]
  if (is.null(mat)) return(NULL)
  mat <- normalize_matrix(mat)
  drugs <- intersect(colnames(mat), drugs)
  if (!length(drugs)) return(NULL)
  vals <- rowMeans(mat[, drugs, drop = FALSE], na.rm = TRUE)
  data.frame(barcode = rownames(mat), score = vals, stringsAsFactors = FALSE)
}

get_dslt_embedding_safe <- function(dslt, level, assay_name, slot) {
  if (is.null(dslt) || !nzchar(level) || !nzchar(assay_name) || !nzchar(slot)) {
    return(NULL)
  }

  tryCatch({
    embeddings <- dslt[["embeddings"]]
    if (is.null(embeddings)) return(NULL)
    level_embeddings <- embeddings[[level]]
    if (is.null(level_embeddings)) return(NULL)
    assay_embeddings <- level_embeddings[[assay_name]]
    if (is.null(assay_embeddings)) return(NULL)
    assay_embeddings[[slot]]
  }, error = function(e) NULL)
}

get_dslt_assay_safe <- function(dslt, level, assay_name) {
  if (is.null(dslt) || !nzchar(level) || !nzchar(assay_name)) {
    return(NULL)
  }

  tryCatch({
    assays <- dslt[["assays"]]
    if (is.null(assays)) return(NULL)
    level_assays <- assays[[level]]
    if (is.null(level_assays)) return(NULL)
    level_assays[[assay_name]]
  }, error = function(e) NULL)
}

get_dslt_column_metadata_safe <- function(dslt, assay_name, field, level = NULL) {
  if (is.null(dslt) || !nzchar(assay_name) || !nzchar(field)) {
    return(NULL)
  }

  column_meta <- tryCatch(dslt[["columnMetadata"]], error = function(e) NULL)
  if (is.null(column_meta)) {
    return(NULL)
  }

  levels_to_check <- if (is.null(level)) list(NULL) else as.list(level)

  for (lvl in levels_to_check) {
    container <- if (is.null(lvl)) {
      column_meta[[assay_name]]
    } else {
      level_key <- as.character(lvl)
      if (is.null(column_meta[[level_key]])) {
        next
      }
      column_meta[[level_key]][[assay_name]]
    }

    if (is.null(container)) {
      next
    }

    # Support both list-style access and data.frames where the field is a column
    value <- container[[field]]
    if (is.null(value) && is.data.frame(container) && field %in% colnames(container)) {
      value <- container[[field]]
    }

    if (!is.null(value)) {
      return(value)
    }
  }

  NULL
}

ensure_archetype_metadata <- function(state, assay_name, level = "lineage") {
  stopifnot(identical(level, "lineage"))
  if (is.null(state$dslt) || !nzchar(assay_name)) {
    return(NULL)
  }

  arch <- get_dslt_column_metadata_safe(state$dslt, assay_name, "archetypes", level = level)
  if (!is.null(arch)) {
    return(arch)
  }

  dslt_result <- tryCatch(
    performArchetypeAndAnnotate(state$dslt, assay_name, levels = level),
    error = function(e) {
      showNotification(paste("Archetype analysis failed:", e$message), type = "error")
      NULL
    }
  )

  if (!is.null(dslt_result)) {
    if (inherits(dslt_result, "DatasetLT") ||
        (is.list(dslt_result) && all(c("assays", "embeddings") %in% names(dslt_result)))) {
      state$dslt <- dslt_result
    } else if (is.matrix(dslt_result) || is.data.frame(dslt_result)) {
      arch <- as.matrix(dslt_result)
  level_key <- "lineage"
      if (is.null(state$dslt[["columnMetadata"]])) state$dslt[["columnMetadata"]] <- list()
      if (is.null(state$dslt[["columnMetadata"]][[level_key]])) state$dslt[["columnMetadata"]][[level_key]] <- list()
      if (is.null(state$dslt[["columnMetadata"]][[level_key]][[assay_name]])) {
        state$dslt[["columnMetadata"]][[level_key]][[assay_name]] <- list()
      }
      state$dslt[["columnMetadata"]][[level_key]][[assay_name]][["archetypes"]] <- arch
    }
  }

  get_dslt_column_metadata_safe(state$dslt, assay_name, "archetypes", level = level)
}

ensure_knn_embeddings <- function(state, assay_name, level = "lineage") {
  if (is.null(state$dslt) || !nzchar(assay_name)) return(list(coords = NULL, clusters = NULL))
  stopifnot(identical(level, "lineage"))
  assay_present <- !is.null(get_dslt_assay_safe(state$dslt, level, assay_name))
  if (!assay_present) {
    return(list(coords = NULL, clusters = NULL))
  }
  coords <- get_dslt_embedding_safe(state$dslt, level, assay_name, "umap")
  clusters <- get_dslt_embedding_safe(state$dslt, level, assay_name, "louvain_clusters")
  if (is.null(coords) || is.null(clusters)) {
    dslt_updated <- tryCatch(
      runKnnAnalysisAndStore(state$dslt, assay_name, levels = level),
      error = function(e) {
        showNotification(paste("KNN analysis failed:", e$message), type = "error")
        NULL
      }
    )
    if (is.null(dslt_updated)) {
      return(list(coords = NULL, clusters = NULL))
    }
    state$dslt <- dslt_updated
    coords <- get_dslt_embedding_safe(state$dslt, level, assay_name, "umap")
    clusters <- get_dslt_embedding_safe(state$dslt, level, assay_name, "louvain_clusters")
  }
  if (is.data.frame(clusters) && "louvain_clusters" %in% colnames(clusters)) {
    clusters <- clusters$louvain_clusters
  }
  list(coords = coords, clusters = clusters)
}

cluster_palette <- c(
  "#ebce2b", "#702c8c", "#db6917", "#96cde6", "#ba1c30", "#c0bd7f",
  "#7f7e80", "#5fa641", "#d485b2", "#4277b6", "#df8461", "#463397",
  "#e1a11a", "#91218c", "#e8e948", "#7e1510", "#92ae31", "#6f340d",
  "#d32b1e", "#2b3514"
)

drug_palette <- rev(c('#67001f','#b2182b',"#f4a582",'#fddbc7',"#ffffff",'#d1e5f0','#92c5de','#2166ac','#053061'))

get_dimensional_reduction <- function(assays) {
  if (tolower(assays) == "atac") {
    "lsi"
  } else {
    "pca"
  }
}

ensure_pca <- function(seu, dslt = NULL, npcs = 30, assays = "RNA", level = "lineage") {
  reduction <- get_dimensional_reduction(assays)
  needs_reduction <- !reduction %in% Reductions(seu)
  if (!needs_reduction) {
    existing <- Embeddings(seu, reduction)
    needs_reduction <- ncol(existing) < npcs
  }
  if (needs_reduction) {
    if (tolower(assays) == "atac") {
      res <- run_svd(seu, dslt = dslt, npcs = npcs, assays = assays, level = level)
    } else {
      res <- run_pca(seu, dslt = dslt, npcs = npcs, assays = assays, level = level)
    }
    seu <- res$seu
    dslt <- res$dslt
  } else if (!is.null(dslt)) {
    dslt <- update_dslt_embedding(dslt, assays, level, reduction, Embeddings(seu, reduction))
  }
  list(seu = seu, dslt = dslt)
}

ensure_umap <- function(seu, dims, dslt = NULL, assays = "RNA", level = "lineage", reduction_name = "umap") {
  dims_seq <- seq(dims[1], dims[2])
  pca_res <- ensure_pca(seu, dslt = dslt, npcs = max(dims_seq), assays = assays, level = level)
  seu <- pca_res$seu
  dslt <- pca_res$dslt
  reduction <- get_dimensional_reduction(assays)
  if (!reduction_name %in% names(Reductions(seu))) {
    res <- run_umap(seu,
      dslt = dslt,
      dimsl = dims_seq[1],
      dimsh = dims_seq[length(dims_seq)],
      assays = assays,
      level = level,
      reduction = reduction
    )
    seu <- res$seu
    dslt <- res$dslt
  } else if (!is.null(dslt)) {
    dslt <- update_dslt_embedding(dslt, assays, level, reduction_name, Embeddings(seu, reduction_name))
  }
  list(seu = seu, reduction = reduction_name, dslt = dslt)
}

ensure_tsne <- function(seu, dims, dslt = NULL, assays = "RNA", level = "lineage", reduction_name = "tsne") {
  dims_seq <- seq(dims[1], dims[2])
  pca_res <- ensure_pca(seu, dslt = dslt, npcs = max(dims_seq), assays = assays, level = level)
  seu <- pca_res$seu
  dslt <- pca_res$dslt
  reduction <- get_dimensional_reduction(assays)
  if (!reduction_name %in% names(Reductions(seu))) {
    res <- run_tsne(seu,
      dslt = dslt,
      dimsl = dims_seq[1],
      dimsh = dims_seq[length(dims_seq)],
      assays = assays,
      level = level,
      reduction = reduction
    )
    seu <- res$seu
    dslt <- res$dslt
  } else if (!is.null(dslt)) {
    dslt <- update_dslt_embedding(dslt, assays, level, reduction_name, Embeddings(seu, reduction_name))
  }
  list(seu = seu, reduction = reduction_name, dslt = dslt)
}

ensure_clusters <- function(seu, dslt = NULL, res = 0.8, dims = 1:30, assays = "RNA", level = "lineage") {
  dims_seq <- seq(dims[1], dims[length(dims)])
  pca_res <- ensure_pca(seu, dslt = dslt, npcs = max(dims_seq), assays = assays, level = level)
  seu <- pca_res$seu
  dslt <- pca_res$dslt
  reduction <- get_dimensional_reduction(assays)
  graph_name <- paste0(DefaultAssay(seu), "_snn")
  if (!graph_name %in% names(seu@graphs) || !"seurat_clusters" %in% colnames(seu@meta.data)) {
    assay_name <- if (!is.null(assays)) toupper(assays) else toupper(assays)
    cluster_args <- list(
      seu = seu,
      dimsl = dims_seq[1],
      dimsh = dims_seq[length(dims_seq)],
      reduction = reduction,
      assays = assays
    )
    if (assay_name != "ATAC") {
      cluster_args$res <- res
    }
    seu <- do.call(run_neighbors_and_clusters, cluster_args)
  }
  if (!is.null(dslt) && "seurat_clusters" %in% colnames(seu@meta.data)) {
    dslt <- add_clusters_to_dslt(dslt, seu, assays = assays, level = level)
  }
  list(seu = seu, dslt = dslt)
}

ensure_kmeans_clusters <- function(seu, k, res = 0.5, dims = 1:30, dslt = NULL, assays = "RNA", level = "lineage") {
  k <- max(2, as.integer(k))
  dims_seq <- seq(dims[1], dims[length(dims)])
  cluster_res <- ensure_clusters(seu,
    dslt = dslt,
    res = res,
    dims = dims_seq,
    assays = assays,
    level = level
  )
  seu <- cluster_res$seu
  dslt <- cluster_res$dslt
  key <- paste0("kmeans", k)
  if (!key %in% colnames(seu@meta.data)) {
    reduction <- get_dimensional_reduction(assays)
    seu <- run_kmeans_clustering(seu,
      knum = k,
      dimsl = dims_seq[1],
      dimsh = dims_seq[length(dims_seq)],
      reduction = reduction
    )
  }
  if (!is.null(dslt)) {
    dslt <- add_clusters_to_dslt(dslt, seu, knum = k, assays = assays, level = level)
  }
  list(seu = seu, column = key, dslt = dslt)
}


#------------------------------------------------seurat object
# ATAC-seq
#' 创建ATAC-seq Seurat对象
#' 
#' @param h5_file 10X h5文件路径
#' @param fragments_file fragments文件路径
#' @param genome 基因组版本，默认"hg38"
#' @param min_cells 最少细胞数，默认10
#' @param min_features 最少特征数，默认200
#' @return Seurat对象，包含ATAC-seq数据和注释信息
create_atac_seurat <- function(h5_file = NULL,
                                counts_data = NULL,
                                fragments_file,
                                genome = "hg38",
                                min_cells = 10,
                                min_features = 200,
                                metadata = NULL) {

  if (!is.null(counts_data)) {
    counts_mat <- counts_data
  } else {
    stopifnot(!is.null(h5_file))
    counts_mat <- Read10X_h5(filename = h5_file)
  }

  chrom_assay <- CreateChromatinAssay(
    counts = counts_mat,
    sep = c(":", "-"),
    genome = genome,
    fragments = fragments_file,
    min.cells = min_cells,
    min.features = min_features
  )

  cell_names <- colnames(chrom_assay)

  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata)

    if (is.null(rownames(metadata))) {
      stop("Metadata must have row names corresponding to cell barcodes.")
    }

    missing_cells <- setdiff(cell_names, rownames(metadata))
    if (length(missing_cells) > 0) {
      warning(
        "Metadata is missing entries for some cells. Corresponding values will be filled with NA.",
        call. = FALSE
      )
    }

    metadata <- metadata[cell_names, , drop = FALSE]
    rownames(metadata) <- cell_names
  }

  # 创建Seurat对象
  seurat_obj <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks",
    meta.data = metadata
  )
  
  # 添加基因注释信息
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- genome
  Annotation(seurat_obj) <- annotations
  
  return(seurat_obj)
}

# 使用示例:
# seurat_obj <- create_atac_seurat(
#   h5_file = "H160_pPer_filtered_peak_bc_matrix.h5",
#   metadata_file = "H160_pPer_singlecell.csv",
#   fragments_file = "H160_pPer_fragments.tsv.gz"
# )

#------------------------------------------------QC
# RNA-seq
#' 创建RNA-seq Seurat对象
#' 
#' @param h5_file 10X RNA h5文件路径
#' @param min_features 最少特征数，默认200
#' @param min_cells 最少细胞数，默认3
#' @return Seurat对象，包含RNA-seq数据
create_rna_seurat <- function(h5_file = NULL,
                               matrix_data = NULL,
                               min_features = 200,
                               min_cells = 3) {

  if (!is.null(matrix_data)) {
    H_RNA <- matrix_data
  } else {
    stopifnot(!is.null(h5_file))
    H_RNA <- Read10X_h5(h5_file)
  }

  # 提取基因表达矩阵
  mat <- if (is.list(H_RNA)) {
    # 常见键名："Gene Expression" 或 "Gene Expression-AB" 等
    nm <- grep("Gene Expression|RNA", names(H_RNA), value = TRUE)[1]
    H_RNA[[nm]]
  } else {
    H_RNA
  }
  
  # 创建Seurat对象
  seurat_obj_rna <- CreateSeuratObject(
    counts = mat,
    assay = "RNA",
    min.features = min_features,
    min.cells = min_cells
  )
  
  return(seurat_obj_rna)
}

# 使用示例:
# seurat_obj_rna <- create_rna_seurat(
#   h5_file = "H160_pPer_GEX_filtered_feature_bc_matrix.h5"
# )

#' 生成RNA质控小提琴图
#'
#' @param se_rna Seurat对象,包含RNA数据
#' @param min_features 最低基因数阈值,默认200
#' @param max_features 最高基因数阈值,默认6000
#' @param min_counts 最低UMI/reads阈值,默认NULL表示不限制
#' @param max_counts 最高UMI/reads阈值,默认NULL表示不限制
#' @param max_percent_mt 最高线粒体基因百分比阈值,默认10
#'
#' @return patchwork对象,包含三个小提琴图的组合
#' @export
#'
#' @examples
#' violin_plots <- plot_rna_qc_violin(se_rna, min_features = 200, max_features = 6000, max_percent_mt = 10)
#' print(violin_plots)
plot_rna_qc_violin <- function(se_rna, 
                               min_features = 200, 
                               max_features = 6000, 
                               max_percent_mt = 10) {
  
  # 计算线粒体基因百分比(如果未计算)
  if (!"percent.mt" %in% colnames(se_rna@meta.data)) {
    se_rna[["percent.mt"]] <- PercentageFeatureSet(se_rna, pattern = "^MT-")
  }
  
  # 1. nFeature_RNA 的小提琴图
  p_features <- VlnPlot(se_rna, features = "nFeature_RNA", pt.size = 0, group.by = "orig.ident") +
    geom_hline(yintercept = min_features, linetype = "dashed", color = "red") +
    geom_hline(yintercept = max_features, linetype = "dashed", color = "red") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend()
  
  # 2. nCount_RNA 的小提琴图
  p_counts <- VlnPlot(se_rna, features = "nCount_RNA", pt.size = 0, group.by = "orig.ident") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
    NoLegend()
  
  # 3. percent.mt 的小提琴图
  p_percent_mt <- VlnPlot(se_rna, features = "percent.mt", pt.size = 0, group.by = "orig.ident") +
    geom_hline(yintercept = max_percent_mt, linetype = "dashed", color = "red") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none") +
    NoLegend()
  
  # 使用 patchwork 组合小提琴图
  violin_plots <- p_features + p_counts + p_percent_mt + plot_layout(ncol = 3)
  
  return(violin_plots)
}


#' 生成RNA质控散点图
#'
#' @param se_rna Seurat对象,包含RNA数据
#'
#' @return patchwork对象,包含两个散点图的组合
#' @export
#'
#' @examples
#' scatter_plots <- plot_rna_qc_scatter(se_rna)
#' print(scatter_plots)
plot_rna_qc_scatter <- function(se_rna) {
  
  # 计算线粒体基因百分比(如果未计算)
  if (!"percent.mt" %in% colnames(se_rna@meta.data)) {
    se_rna[["percent.mt"]] <- PercentageFeatureSet(se_rna, pattern = "^MT-")
  }
  
  # 散点图1: nCount_RNA vs percent.mt
  plot1 <- FeatureScatter(se_rna, feature1 = "nCount_RNA", feature2 = "percent.mt")
  
  # 散点图2: nCount_RNA vs nFeature_RNA
  plot2 <- FeatureScatter(se_rna, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  
  # 使用 patchwork 组合散点图
  scatter_plots <- plot1 + plot2
  
  return(scatter_plots)
}


#' 执行RNA质控过滤
#'
#' @param se_rna Seurat对象,包含RNA数据
#' @param min_features 最低基因数阈值,默认200
#' @param max_features 最高基因数阈值,默认6000
#' @param max_percent_mt 最高线粒体基因百分比阈值,默认10
#' @param verbose 是否打印过滤信息,默认TRUE
#'
#' @return 过滤后的Seurat对象
#' @export
#'
#' @examples
#' se_rna_filtered <- filter_rna_qc(se_rna, min_features = 200, max_features = 6000, max_percent_mt = 10)
filter_rna_qc <- function(se_rna,
                          min_features = 200,
                          max_features = 6000,
                          min_counts = NULL,
                          max_counts = NULL,
                          max_percent_mt = 10,
                          verbose = TRUE) {
  
  # 计算线粒体基因百分比(如果未计算)
  if (!"percent.mt" %in% colnames(se_rna@meta.data)) {
    se_rna[["percent.mt"]] <- PercentageFeatureSet(se_rna, pattern = "^MT-")
  }
  
  # 打印过滤前细胞数量
  if (verbose) {
    print(paste0('过滤前的细胞数量为：', ncol(se_rna)))
  }
  
  keep <- rep(TRUE, ncol(se_rna))

  if (!is.null(min_features)) {
    keep <- keep & !is.na(se_rna$nFeature_RNA) & se_rna$nFeature_RNA > min_features
  }

  if (!is.null(max_features)) {
    keep <- keep & !is.na(se_rna$nFeature_RNA) & se_rna$nFeature_RNA < max_features
  }

  if (!is.null(min_counts)) {
    keep <- keep & !is.na(se_rna$nCount_RNA) & se_rna$nCount_RNA > min_counts
  }

  if (!is.null(max_counts)) {
    keep <- keep & !is.na(se_rna$nCount_RNA) & se_rna$nCount_RNA < max_counts
  }

  if (!is.null(max_percent_mt)) {
    keep <- keep & !is.na(se_rna$percent.mt) & se_rna$percent.mt < max_percent_mt
  }

  se_rna <- subset(se_rna, cells = colnames(se_rna)[keep])
  
  # 打印过滤后细胞数量
  if (verbose) {
    print(paste0('过滤后的细胞数量为：', ncol(se_rna)))
  }
  
  return(se_rna)
}

#' 计算ATAC数据的核小体信号
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @return 更新后的Seurat对象,包含核小体信号和分组信息
#' @export
calculate_nucleosome_signal <- function(se_atac) {
  se_atac <- NucleosomeSignal(object = se_atac)
  se_atac$nucleosome_group <- ifelse(se_atac$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
  return(se_atac)
}

#' 生成Fragment长度分布柱状图
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @return ggplot图形对象
#' @export
plot_fragment_histogram <- function(se_atac) {
  p <- FragmentHistogram(object = se_atac, group.by = 'nucleosome_group')
  return(p)
}

#' 计算TSS富集分数
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @return 更新后的Seurat对象,包含TSS富集分数
#' @export
calculate_tss_enrichment <- function(se_atac) {
  plan("sequential")
  se_atac <- TSSEnrichment(object = se_atac)
  return(se_atac)
}

#' 生成nCount_peaks与TSS enrichment的密度散点图
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @return ggplot图形对象
#' @export
plot_tss_density_scatter <- function(se_atac) {
  p <- DensityScatter(se_atac, x = 'nCount_peaks', y = 'TSS_fragments', 
                      log_x = TRUE, quantiles = TRUE)
  return(p)
}

#' 计算ATAC质控指标
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @param blacklist_regions 黑名单区域的GRanges对象
#' @return 更新后的Seurat对象,包含pct_reads_in_peaks和blacklist_ratio
#' @export
calculate_atac_qc_metrics <- function(se_atac, blacklist_regions = blacklist_hg38_unified) {
  # 计算 pct_reads_in_peaks比例
  se_atac$pct_reads_in_peaks <- se_atac$peak_region_fragments / se_atac$passed_filters * 100

  # 计算黑名单区域比例
  se_atac$blacklist_ratio <- FractionCountsInRegion(
    object = se_atac,
    assay = 'peaks',
    regions = blacklist_regions
  )

  return(se_atac)
}

#' 生成ATAC质控小提琴图组合
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @param ncount_min nCount_peaks最小值,默认2000
#' @param ncount_max nCount_peaks最大值,默认30000
#' @param pct_reads_min pct_reads_in_peaks最小值,默认20
#' @param tss_min TSS enrichment最小值,默认2
#' @param nucleosome_max nucleosome signal最大值,默认4
#' @param blacklist_max blacklist ratio最大值,默认0.02
#' @return 组合的ggplot图形对象
#' @export
plot_atac_qc_violins <- function(se_atac, 
                                  ncount_min = 2000,
                                  ncount_max = 30000,
                                  pct_reads_min = 20,
                                  tss_min = 2,
                                  nucleosome_max = 4,
                                  blacklist_max = 0.02) {
  
  # 生成各个小提琴图
  atsc1 <- VlnPlot(object = se_atac, features = 'nCount_peaks') +
    geom_hline(yintercept = ncount_min, linetype = "dashed", color = "red") +
    geom_hline(yintercept = ncount_max, linetype = "dashed", color = "blue")
  
  atsc2 <- VlnPlot(object = se_atac, features = 'pct_reads_in_peaks') +
    geom_hline(yintercept = pct_reads_min, linetype = "dashed", color = "red")
  
  atsc3 <- VlnPlot(object = se_atac, features = 'TSS_fragments') +
    geom_hline(yintercept = tss_min, linetype = "dashed", color = "red")
  
  atsc4 <- VlnPlot(object = se_atac, features = 'nucleosome_signal') +
    geom_hline(yintercept = nucleosome_max, linetype = "dashed", color = "red")
  
  atsc5 <- VlnPlot(object = se_atac, features = 'blacklist_ratio') +
    geom_hline(yintercept = blacklist_max, linetype = "dashed", color = "red")
  
  # 组合所有图
  combined_plot <- atsc1 + atsc2 + atsc3 + atsc4 + atsc5 +
    plot_layout(ncol = 5, guides = 'collect')
  
  return(combined_plot)
}

#' 过滤ATAC数据的低质量细胞
#'
#' @param se_atac Seurat对象,包含ATAC-seq数据
#' @param ncount_min nCount_peaks最小值,默认2000
#' @param ncount_max nCount_peaks最大值,默认30000
#' @param pct_reads_min pct_reads_in_peaks最小值,默认20
#' @param blacklist_max blacklist_ratio最大值,默认0.02
#' @param nucleosome_max nucleosome_signal最大值,默认4
#' @param tss_min TSS_fragments最小值,默认2
#' @param suffix 细胞后缀过滤,默认"-1$"
#' @param verbose 是否打印过滤信息,默认TRUE
#' @return 过滤后的Seurat对象
#' @export
filter_atac_cells <- function(se_atac, 
                               ncount_min = 2000,
                               ncount_max = 30000,
                               pct_reads_min = 20,
                               blacklist_max = 0.02,
                               nucleosome_max = 4,
                               tss_min = 2,
                               verbose = TRUE) {
  
  if (verbose) {
    print(paste0('过滤前的细胞数量为：', ncol(se_atac)))
  }
  
  # 根据质控指标过滤
  se_atac <- subset(
    x = se_atac,
    subset = nCount_peaks > ncount_min & 
      nCount_peaks < ncount_max & 
      pct_reads_in_peaks > pct_reads_min & 
      blacklist_ratio < blacklist_max & 
      nucleosome_signal < nucleosome_max & 
      TSS_fragments > tss_min
  )
  
  if (verbose) {
    print(paste0('过滤后的细胞数量为：', ncol(se_atac)))
  }
  
  return(se_atac)
}

#------------------------------------------------multi-omics integration、filtering
#' 将lineage数据映射到single-cell数据
#' 将lineage数据映射到Seurat RNA对象
#'
#' 该函数将单细胞RNA数据与lineage barcode信息关联,并按barcode聚合细胞,
#' 创建pseudo-bulk数据集用于后续分析
#'
#' @param se Seurat对象,包含单细胞RNA数据
#' @param bc 数据框,包含barcode映射信息,必须包含'cell_barcode'和'Barcode'列
#' @param dslt 目标对象,用于存储映射信息
#' @param assay 要使用的assay名称,默认为"RNA"
#' @param key 映射结果的键名前缀,将创建名为'key_map'的assay
#'
#' @return 返回一个新的Seurat对象,包含按barcode聚合后的pseudo-bulk数据,
#'         已经过SCTransform标准化
#'
#' @export
#'
#' @examples
lineage_map_seurat_rna <- function(se, bc, dslt, assay, key){
  m  <- LayerData(se, assay = assay, layer = "counts")
  md <- bc[match(colnames(m), bc$cell_barcode), , drop = FALSE]
  keep <- !is.na(md$Barcode); md <- md[keep, , drop = FALSE]; m <- m[, keep, drop = FALSE]
  rownames(md) <- md$cell_barcode
  sce <- SingleCellExperiment(list(counts = m), colData = md)
  pb  <- aggregateData(sce, assay = "counts", by = "Barcode")
  map <- transform(md, n_cells = as.integer(table(Barcode)[Barcode]))[, c("Barcode","cell_barcode","n_cells")]
  dslt[["assays"]][[paste0(key, "_map")]] <- map
  mat_rna  <- assays(pb)[[1]]
  pb_rna  <- CreateSeuratObject(counts = mat_rna)
  pb_rna <- SCTransform(pb_rna, verbose = FALSE)
  pb_rna
}

#' 将lineage数据映射到single-cell ATAC数据
#' 将lineage数据映射到Seurat ATAC对象
#'
#' 该函数将单细胞ATAC数据与lineage barcode信息关联,并按barcode聚合细胞,
#' 创建pseudo-bulk数据集用于后续分析
#'
#' @param se Seurat对象,包含单细胞ATAC数据
#' @param bc 数据框,包含barcode映射信息,必须包含'cell_barcode'和'Barcode'列
#' @param dslt 目标对象,用于存储映射信息
#' @param assay 要使用的assay名称,默认为"peaks"
#' @param key 映射结果的键名前缀,将创建名为'key_map'的assay
#'
#' @return 返回一个新的Seurat对象,包含按barcode聚合后的pseudo-bulk ATAC数据,
#'         已经过TF-IDF和SVD标准化
#'
#' @export
#'
#' @examples
lineage_map_seurat_atac <- function(se, bc, dslt, assay = "peaks", key){
  # 提取counts数据
  m  <- LayerData(se, assay = assay, layer = "counts")
  
  # 匹配barcode
  md <- bc[match(colnames(m), bc$cell_barcode), , drop = FALSE]
  
  # 过滤未匹配的细胞
  keep <- !is.na(md$Barcode)
  md <- md[keep, , drop = FALSE]
  m <- m[, keep, drop = FALSE]
  rownames(md) <- md$cell_barcode
  
  # 创建SingleCellExperiment对象并聚合
  sce <- SingleCellExperiment(list(counts = m), colData = md)
  pb  <- aggregateData(sce, assay = "counts", by = "Barcode")
  
  # 创建映射表
  map <- transform(md, n_cells = as.integer(table(Barcode)[Barcode]))[, c("Barcode","cell_barcode","n_cells")]
  dslt[["assays"]][[paste0(key, "_map")]] <- map
  
  # 提取聚合后的矩阵
  mat_atac <- assays(pb)[[1]]
  
  # 创建ChromatinAssay对象
  chrom_assay <- CreateChromatinAssay(
    counts = mat_atac,
    sep = c(":", "-"),
    genome = se@assays[[assay]]@genome,
    min.cells = 0,
    min.features = 0
  )
  
  # 创建Seurat对象
  pb_atac <- CreateSeuratObject(
    counts = chrom_assay,
    assay = "peaks"
  )

  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome(annotations) <- "hg38"
  Annotation(pb_atac) <- annotations
  rm(annotations)

  # ATAC数据标准化: TF-IDF + SVD
  pb_atac <- RunTFIDF(pb_atac, verbose = FALSE)
  pb_atac <- FindTopFeatures(pb_atac, min.cutoff = 'q0', verbose = FALSE)
  pb_atac <- RunSVD(pb_atac, verbose = FALSE)
  
  return(pb_atac)
}



#' 根据细胞barcode后缀过滤细胞
#'
#' 该函数根据指定的后缀模式过滤Seurat对象中的细胞。
#' 常用于多样本整合分析中保留特定样本的细胞。
#'
#' @param obj Seurat对象或其他包含细胞barcode的对象
#' @param suffix 字符串,用于匹配的正则表达式后缀,默认为"-1$"
#'               表示保留以"-1"结尾的细胞barcode
#'
#' @return 过滤后的对象,仅包含匹配指定后缀的细胞
#'
#' @details
#' 该函数会:
#' \itemize{
#'   \item 使用grep函数查找匹配后缀的细胞barcode
#'   \item 如果没有找到匹配的细胞,抛出错误
#'   \item 如果匹配的细胞数少于总细胞数,则进行subset操作
#'   \item 如果所有细胞都匹配,则返回原对象
#' }
#'
#' @export
#'
#' @examples
filter_cells_by_suffix <- function(obj, suffix = "-1$") {
  keep_cells <- grep(suffix, colnames(obj), value = TRUE)
  
  if (!length(keep_cells)) {
    stop(paste0("未找到以 ", gsub("\\$", "", suffix), " 结尾的细胞条形码。"))
  }
  
  if (length(keep_cells) < ncol(obj)) {
    obj <- subset(obj, cells = keep_cells)
  }
  
  return(obj)
}

map_lineage_to_single_cell <- function(dslt, mapping, assays = NULL, umi_col = "UMI", suffix = NULL) {
  lineage_assays <- dslt[["assays"]][["lineage"]]

  available_assays <- names(lineage_assays)
  if (is.null(available_assays)) {
    available_assays <- character(0)
  }

  if (is.null(assays)) {
    assays_to_process <- available_assays
  } else {
    assays_to_process <- intersect(assays, available_assays)
  }

  if (!length(assays_to_process)) {
    return(dslt)
  }
  
  mapping <- mapping[
    !is.na(mapping[["cell_barcode"]]) & mapping[["cell_barcode"]] != "" &
      !is.na(mapping[["Barcode"]]) & mapping[["Barcode"]] != "",
    ,
    drop = FALSE
  ]

  if (!is.null(suffix) && nzchar(suffix)) {
    keep_cells <- grepl(suffix, mapping[["cell_barcode"]])
    mapping <- mapping[keep_cells, , drop = FALSE]
    if (nrow(mapping) == 0) {
      empty_assays <- lapply(lineage_assays[assays_to_process], function(mat) {
        mat <- as.matrix(mat)
        empty_mat <- mat[0, , drop = FALSE]
        empty_df <- as.data.frame(empty_mat)
        rownames(empty_df) <- rownames(empty_mat)
        colnames(empty_df) <- colnames(empty_mat)
        empty_df
      })
      existing_single_cell <- dslt[["assays"]][["single_cell"]]
      if (is.null(existing_single_cell)) {
        existing_single_cell <- list()
      }
      existing_single_cell[assays_to_process] <- empty_assays
      dslt[["assays"]][["single_cell"]] <- existing_single_cell
      return(dslt)
    }
  }
  mapping[[umi_col]] <- as.character(mapping[[umi_col]])

  umi_counts <- stats::aggregate(
    mapping[[umi_col]],
    by = list(cell_barcode = mapping[["cell_barcode"]], Barcode = mapping[["Barcode"]]),
    FUN = function(x) length(unique(x))
  )
  names(umi_counts)[3] <- "umi_count"
  
  single_cell_assays <- lapply(lineage_assays[assays_to_process], function(mat) {
    mat <- as.matrix(mat)
    idx <- match(umi_counts[["Barcode"]], rownames(mat))
    keep <- !is.na(idx)
    if (!any(keep)) {
      empty_mat <- mat[0, , drop = FALSE]
      empty_df <- as.data.frame(empty_mat)
      rownames(empty_df) <- rownames(empty_mat)
      colnames(empty_df) <- colnames(empty_mat)
      return(empty_df)
    }
    counts_use <- umi_counts[keep, , drop = FALSE]
    wts <- counts_use[["umi_count"]]
    values <- mat[idx[keep], , drop = FALSE] * wts
    summed <- rowsum(values, counts_use[["cell_barcode"]])
    wsum <- rowsum(wts, counts_use[["cell_barcode"]])
    wsum_vec <- as.numeric(wsum)
    names(wsum_vec) <- rownames(wsum)
    averaged <- sweep(summed, 1, wsum_vec, "/")
    averaged_df <- as.data.frame(averaged)
    rownames(averaged_df) <- rownames(averaged)
    colnames(averaged_df) <- colnames(averaged)
    averaged_df
  })

  existing_single_cell <- dslt[["assays"]][["single_cell"]]
  if (is.null(existing_single_cell)) {
    existing_single_cell <- list()
  } else {
    existing_single_cell <- lapply(existing_single_cell, function(entry) {
      if (is.matrix(entry) || is.data.frame(entry)) {
        df <- as.data.frame(entry)
        rownames(df) <- rownames(entry)
        colnames(df) <- colnames(entry)
        return(df)
      }
      entry
    })
  }

  existing_single_cell[assays_to_process] <- single_cell_assays
  dslt[["assays"]][["single_cell"]] <- existing_single_cell

  dslt
}

#------------------------------------------------Denoizing and archetype
#' @title Initialize DatasetLT from LMM output
#' @description Reads the RDS file for the specified cell line and prepares a DatasetLT object with active samples.
#' @param cell_line_id Character scalar cell line identifier.
#' @param output_dir Character scalar path to the LMM results.
#' @return A `DatasetLT` object with active samples set.
initializeDsltFromLmm <- function(output_dir,threshold = 2e-5) {
  lmm_result_path <- output_dir
  lmm_result_object <- readRDS(lmm_result_path)
  dslt <- DatasetLT$new()
  dslt[["assays"]][["lineage"]] <- lapply(lmm_result_object$result, function(entry) {
    if (is.matrix(entry) || is.data.frame(entry)) {
      df <- as.data.frame(entry)
      rownames(df) <- rownames(entry)
      colnames(df) <- colnames(entry)
      return(df)
    }
    entry
  })
  rm(lmm_result_object)
  dslt$printLayerNames()
  fraction_assays <- dslt[["assays"]][["lineage"]]
  fraction_assay <- fraction_assays[["fraction"]]
  fraction_assay <- as.matrix(fraction_assay)
  active_sample_ids <- rownames(fraction_assay)[(matrixStats::rowMaxs(fraction_assay) > threshold)]
  dslt$setActiveSamples(active_sample_ids)
  dslt
}

#' @title Apply adaptive kernel denoising
#' @description Scales the selected assay and stores the smoothed assay obtained via adaptive kernel denoising.
#' @param dslt A `DatasetLT` instance to be updated.
#' @param assay_name Character scalar with the assay to analyse.
#' @param smoothed_assay_name Character scalar for the smoothed assay name.
#' @param center_flag Logical indicating centering.
#' @param scale_flag Logical indicating scaling.
#' @param k_neighbors_prop_value Numeric proportion of neighbours.
#' @param snn_threshold_prop_value Numeric SNN threshold proportion.
#' @param gamma_value Numeric gamma parameter.
#' @return The updated `DatasetLT` object.
applyAdaptiveKernelDenoising <- function(dslt,
                                         assay_name,
                                         center_flag = TRUE,
                                         scale_flag = TRUE,
                                         k_neighbors_prop_value = 0.025,
                                         snn_threshold_prop_value = 0.2,
                                         gamma_value = 1) {
  smoothed_assay_name <- paste0(assay_name, "_smoothed")                                        
  lineage_assays <- dslt[["assays"]][["lineage"]]
  assay_data <- lineage_assays[[assay_name]]
  scaled_matrix <- scale2(as.matrix(assay_data), center = center_flag, scale = scale_flag)
  KD_mat <- AdaptiveKernelDenoizing(scaled_matrix,
                                     use_ica = TRUE,
                                     k_neighbors_prop = k_neighbors_prop_value,
                                     snn_threshold_prop = snn_threshold_prop_value,
                                     gamma = gamma_value)
  KD_df <- as.data.frame(KD_mat)
  rownames(KD_df) <- rownames(KD_mat)
  colnames(KD_df) <- colnames(KD_mat)
  dslt[["assays"]][["lineage"]][[smoothed_assay_name]] <- KD_df
  list(dslt = dslt, smoothed_assay_name = smoothed_assay_name)
}

#' @title Perform archetype analysis and annotate DatasetLT
#' @description Runs archetype analysis on the smoothed assay and stores embeddings and metadata.
#' @param dslt A `DatasetLT` instance to be updated.
#' @param smoothed_assay_name Character scalar for the smoothed assay name.
#' @param center_flag Logical indicating centering.
#' @param scale_flag Logical indicating scaling.
#' @param kappa_value Numeric kappa parameter.
#' @param worker_count Integer number of workers.
#' @param projected_count Integer number of projected dimensions.
#' @return The updated `DatasetLT` object.
performArchetypeAndAnnotate <- function(dslt,
                                        smoothed_assay_name,
                                        levels = NULL,
                                        center_flag = TRUE,
                                        scale_flag = TRUE,
                                        kappa_value = 9,
                                        worker_count = 20,
                                        projected_count = 4) {
  if (levels == "lineage") {
    active_sample_ids <- rownames(dslt[["assays"]][["lineage"]][[smoothed_assay_name]])
    dslt$setActiveSamples(active_sample_ids)
    archetype_input <- as.data.frame(dslt[["assays"]][["lineage"]][[smoothed_assay_name]])
    archetype_input <- scale2(archetype_input, center = center_flag, scale = scale_flag)
    archetype_results <- find_params_and_perform_arch(archetype_input,
                                                      kappa = kappa_value,
                                                      nworkers = worker_count,
                                                      nprojected = projected_count)
    dslt[["embeddings"]][["lineage"]][[smoothed_assay_name]][["archetype_alpha"]] <- archetype_results$A
    dslt[["columnMetadata"]][["lineage"]][[smoothed_assay_name]][["archetypes"]] <- t(archetype_results$BY)
  } else {
    active_sample_ids <- rownames(state$seurat$sc_rna@meta.data)
    dslt$setActiveSamples(active_sample_ids)
    archetype_input <- as.data.frame(dslt[["assays"]][["single_cell"]][[smoothed_assay_name]])
    archetype_input <- scale2(archetype_input, center = center_flag, scale = scale_flag)
    archetype_results <- find_params_and_perform_arch(archetype_input,
                                                      kappa = kappa_value,
                                                      nworkers = worker_count,
                                                      nprojected = projected_count)
    dslt[["embeddings"]][["single_cell"]][[smoothed_assay_name]][["archetype_alpha"]] <- archetype_results$A
    dslt[["columnMetadata"]][["single_cell"]][[smoothed_assay_name]][["archetypes"]] <- t(archetype_results$BY)
  }
  dslt
}

#------------------------------------------------Plotting preparation, clustering

#' @title Run KNN analysis and store results
#' @description Performs KNN analysis on the smoothed assay and stores graphs and embeddings in the DatasetLT object.
#' @param dslt A `DatasetLT` instance to be updated.
#' @param smoothed_assay_name Character scalar for the smoothed assay name.
#' @param k_neighbors_prop_value Numeric proportion of neighbours.
#' @param snn_threshold_prop_value Numeric SNN threshold proportion.
#' @param gamma_value Numeric gamma parameter.
#' @param min_neighbor_count Numeric minimum neighbours (kept for compatibility).
#' @return The updated `DatasetLT` object.
runKnnAnalysisAndStore <- function(dslt,
                                   smoothed_assay_name,
                                   levels = NULL,
                                   k_neighbors_prop_value = 0.025,
                                   snn_threshold_prop_value = 0.2,
                                   gamma_value = 1,
                                   min_neighbor_count = 10) {
  if (levels == "lineage") {
    active_sample_ids <- rownames(dslt[["assays"]][["lineage"]][[smoothed_assay_name]])
    dslt$setActiveSamples(active_sample_ids)

    lineage_assay <- dslt[["assays"]][["lineage"]][[smoothed_assay_name]]
    knn_results_object <- runKnnAnalysis(as.matrix(lineage_assay),
                                        k_neighbors_prop = k_neighbors_prop_value,
                                        snn_threshold_prop = snn_threshold_prop_value,
                                        gamma = gamma_value)
    graph_components <- knn_results_object[c("similarity_full", "similarity_snn", "adjacency_snn", "distance_matrix")]
    graph_components <- lapply(graph_components, function(mat) {
      if (is.null(mat)) {
        return(NULL)
      }
      mat_df <- as.data.frame(mat)
      if (!is.null(rownames(mat))) {
        rownames(mat_df) <- rownames(mat)
      }
      if (!is.null(colnames(mat))) {
        colnames(mat_df) <- colnames(mat)
      }
      mat_df
    })
    dslt[["graphs"]][["lineage"]][[smoothed_assay_name]] <- graph_components
    dslt[["embeddings"]][["lineage"]][[smoothed_assay_name]][["louvain_clusters"]] <- data.frame(louvain_clusters = (knn_results_object[[c("louvain_clusters")]] ))
    dist_snn <- distance_form_similarity_log(as.matrix(dslt[["graphs"]][["lineage"]][[smoothed_assay_name]][["similarity_snn"]]))
    dist_snn_df <- as.data.frame(dist_snn)
    rownames(dist_snn_df) <- rownames(dist_snn)
    colnames(dist_snn_df) <- colnames(dist_snn)
    dslt[["graphs"]][["lineage"]][[smoothed_assay_name]][["dist_snn"]] <- dist_snn_df
    umap_embedding <- uwot::umap(dist_snn, n_neighbors = 150)
    umap_df <- as.data.frame(umap_embedding)
    rownames(umap_df) <- rownames(dist_snn)
    dslt[["embeddings"]][["lineage"]][[smoothed_assay_name]][["umap"]] <- umap_df
  } else {
    active_sample_ids <- rownames(state$seurat$sc_rna@meta.data)
    dslt$setActiveSamples(active_sample_ids)

    single_cell_assay <- dslt[["assays"]][["single_cell"]][[smoothed_assay_name]]
    knn_results_object <- runKnnAnalysis(as.matrix(single_cell_assay),
                                        k_neighbors_prop = k_neighbors_prop_value,
                                        snn_threshold_prop = snn_threshold_prop_value,
                                        gamma = gamma_value)
    graph_components <- knn_results_object[c("similarity_full", "similarity_snn", "adjacency_snn", "distance_matrix")]
    graph_components <- lapply(graph_components, function(mat) {
      if (is.null(mat)) {
        return(NULL)
      }
      mat_df <- as.data.frame(mat)
      if (!is.null(rownames(mat))) {
        rownames(mat_df) <- rownames(mat)
      }
      if (!is.null(colnames(mat))) {
        colnames(mat_df) <- colnames(mat)
      }
      mat_df
    })
    dslt[["graphs"]][["single_cell"]][[smoothed_assay_name]] <- graph_components
    dslt[["embeddings"]][["single_cell"]][[smoothed_assay_name]][["louvain_clusters"]] <- data.frame(louvain_clusters = (knn_results_object[[c("louvain_clusters")]] ))
    dist_snn <- distance_form_similarity_log(as.matrix(dslt[["graphs"]][["single_cell"]][[smoothed_assay_name]][["similarity_snn"]]))
    dist_snn_df <- as.data.frame(dist_snn)
    rownames(dist_snn_df) <- rownames(dist_snn)
    colnames(dist_snn_df) <- colnames(dist_snn)
    dslt[["graphs"]][["single_cell"]][[smoothed_assay_name]][["dist_snn"]] <- dist_snn_df
    umap_embedding <- uwot::umap(dist_snn, n_neighbors = 150)
    umap_df <- as.data.frame(umap_embedding)
    rownames(umap_df) <- rownames(dist_snn)
    dslt[["embeddings"]][["single_cell"]][[smoothed_assay_name]][["umap"]] <- umap_df
  }
  dslt
}

update_dslt_embedding <- function(dslt, assays, level, name, embedding) {
  if (is.null(dslt) || is.null(embedding)) {
    return(dslt)
  }

  tryCatch({
    embeddings <- dslt$embeddings
    if (is.null(embeddings)) embeddings <- list()
    if (is.null(embeddings[[assays]])) embeddings[[assays]] <- list()
    if (is.null(embeddings[[assays]][[level]])) embeddings[[assays]][[level]] <- list()
    embedding_df <- as.data.frame(embedding)
    rownames(embedding_df) <- rownames(embedding)
    colnames(embedding_df) <- colnames(embedding)
    embeddings[[assays]][[level]][[name]] <- embedding_df
    dslt$embeddings <- embeddings
    dslt
  }, error = function(e) {
    dslt
  })
}

#' 运行PCA降维
#'
#' @param seu Seurat对象
#' @param npcs PCA主成分数量,默认50
#' @param verbose 是否显示详细信息,默认FALSE
#'
#' @return 列表,包含更新后的Seurat对象和(可选)更新的dslt对象
#' @export
run_pca <- function(seu, dslt = NULL, npcs = 50, assays = "RNA", level = "lineage", verbose = FALSE) {
  seu <- RunPCA(seu, npcs = npcs)
  dslt <- update_dslt_embedding(dslt, assays, level, "pca", Embeddings(seu, "pca"))
  list(seu = seu, dslt = dslt)
}

run_svd <- function(seu, dslt = NULL, npcs = 50, assays = "ATAC", level = "lineage", verbose = FALSE) {
  seu <- RunSVD(seu, n = npcs, reduction.name = "lsi", verbose = verbose)
  dslt <- update_dslt_embedding(dslt, assays, level, "lsi", Embeddings(seu, "lsi"))
  list(seu = seu, dslt = dslt)
}

#' 运行UMAP降维并添加到dslt
#'
#' @param seu Seurat对象
#' @param dimsl PCA维度起始值,默认1
#' @param dimsh PCA维度结束值,默认30
#' @param verbose 是否显示详细信息,默认FALSE
#'
#' @return 列表,包含更新后的Seurat对象和(可选)更新的dslt对象
#' @export
run_umap <- function(seu, dslt = NULL, dimsl = 1, dimsh = 30, assays = "RNA", level = "single cell", verbose = FALSE, reduction = "pca") {
  seu <- RunUMAP(seu, dims = dimsl:dimsh, verbose = verbose, reduction = reduction)
  dslt <- update_dslt_embedding(dslt, assays, level, "umap", Embeddings(seu, "umap"))
  list(seu = seu, dslt = dslt)
}


#' 运行TSNE降维并添加到dslt
#'
#' @param seu Seurat对象
#' @param dimsl PCA维度起始值,默认1
#' @param dimsh PCA维度结束值,默认30
#' @param verbose 是否显示详细信息,默认FALSE
#'
#' @return 列表,包含更新后的Seurat对象和(可选)更新的dslt对象
#' @export
run_tsne <- function(seu, dslt = NULL, dimsl = 1, dimsh = 30, assays = "RNA", level = "single cell", verbose = FALSE, reduction = "pca") {
  seu <- RunTSNE(seu, dims = dimsl:dimsh, verbose = verbose, reduction = reduction)
  dslt <- update_dslt_embedding(dslt, assays, level, "tsne", Embeddings(seu, "tsne"))
  list(seu = seu, dslt = dslt)
}


#' 运行FindNeighbors和FindClusters
#'
#' @param seu Seurat对象
#' @param res 聚类分辨率参数(仅适用于RNA等需要设置分辨率的assay)
#' @param dimsl PCA维度起始值,默认1
#' @param dimsh PCA维度结束值,默认30
#'
#' @return 包含邻居图和聚类结果的Seurat对象
#' @export
run_neighbors_and_clusters <- function(seu, res = NULL, dimsl = 1, dimsh = 30, reduction = "pca", assays = NULL) {
  seu <- FindNeighbors(seu, dims = dimsl:dimsh, reduction = reduction)
  assay_name <- if (!is.null(assays)) {
    toupper(assays)
  } else {
    toupper(assays)
  }

  if (assay_name == "ATAC") {
    seu <- FindClusters(seu, verbose = FALSE, algorithm = 3)
  } else {
    resolution <- if (!is.null(res)) {
      res
    } else {
      0.8
    }
    seu <- FindClusters(seu, resolution = resolution)
  }
  seu
}

#' 运行kmeans聚类
#'
#' @param seu Seurat对象
#' @param knum kmeans聚类中心数量,默认5
#' @param dimsl PCA维度起始值,默认1
#' @param dimsh PCA维度结束值,默认30
#' @param nstart kmeans随机起始次数,默认20
#'
#' @return 包含kmeans聚类结果的Seurat对象
#' @export
run_kmeans_clustering <- function(seu, knum = 5, dimsl = 1, dimsh = 30, nstart = 20, reduction = "pca") {
  pc <- Embeddings(seu, reduction)[, dimsl:dimsh]
  km <- kmeans(pc, centers = knum, iter.max = 1000, nstart = nstart)
  seu[[paste0("kmeans", knum)]] <- as.factor(km$cluster)
  seu
}

#' 将聚类结果添加到dslt的columnMetadata
#'
#' @param dslt DatasetLT对象
#' @param sc_seu 单细胞Seurat对象
#' @param pb_seu pseudo-bulk Seurat对象
#' @param knum kmeans聚类中心数量,默认5
#' @param assays assay名称,默认"RNA"
#'
#' @return 更新后的dslt对象
#' @export
add_clusters_to_dslt <- function(dslt, sc_seu, knum = 5, assays = "RNA", level = "single cell") {
  colmeta <- dslt$columnMetadata
  if (is.null(colmeta)) colmeta <- list()
  if (is.null(colmeta[[assays]])) colmeta[[assays]] <- list()

  row_names_sc <- rownames(sc_seu@meta.data)
  target_cols <- c(paste0("kmeans", knum), "seurat_clusters")
  available_cols <- intersect(target_cols, colnames(sc_seu@meta.data))
  selected_cols_sc <- sc_seu@meta.data[, available_cols, drop = FALSE]
  new_df_sc <- data.frame(Barcode = row_names_sc, selected_cols_sc, stringsAsFactors = FALSE)
  if (paste0("kmeans", knum) %in% available_cols) {
    colnames(new_df_sc)[colnames(new_df_sc) == paste0("kmeans", knum)] <- "k"
  }
  if ("seurat_clusters" %in% available_cols) {
    colnames(new_df_sc)[colnames(new_df_sc) == "seurat_clusters"] <- "s"
  }

  existing <- colmeta[[assays]][[level]]
  if (!is.null(existing)) {
    if (!"Barcode" %in% colnames(existing)) {
      existing <- data.frame(Barcode = rownames(existing), existing,
                             stringsAsFactors = FALSE, check.names = FALSE)
    }
    keep_cols <- setdiff(colnames(existing), c("k", "s", "k.new", "s.new"))
    existing <- existing[, keep_cols, drop = FALSE]
    merged <- merge(existing, new_df_sc, by = "Barcode", all = TRUE, suffixes = c("", ".new"))
    if ("k.new" %in% colnames(merged)) {
      merged$k <- ifelse(is.na(merged$k.new), merged$k, merged$k.new)
      merged$k.new <- NULL
    }
    if ("s.new" %in% colnames(merged)) {
      merged$s <- ifelse(is.na(merged$s.new), merged$s, merged$s.new)
      merged$s.new <- NULL
    }
    new_df_sc <- merged
  }

  colmeta[[assays]][[level]] <- new_df_sc
  dslt$columnMetadata <- colmeta

  dslt
}

#------------------------------------------------Plotting

ggscatter_single = function(coords, values, column_name, 
                           ggObj = ggplot(), 
                           size_mult = 1,
                           colors = NULL,
                           symmQuant = NULL, 
                           color_text = "value",
                           legend.position = "bottom",
                           gg_theme = NULL) {
  
  common_rows <- intersect(rownames(coords), rownames(values))
  if (length(common_rows) == 0) {
    stop("No shared rownames between coords and values; cannot align data.")
  }
  coords_aligned <- coords[common_rows, , drop = FALSE]
  values_aligned <- values[common_rows, column_name, drop = TRUE]

  ggscatter_colored(coords = coords_aligned,
                   values = values_aligned,
                   ggObj = ggObj,
                   size_mult = size_mult,
                   colors = colors,
                   symmQuant = symmQuant,
                   gg_theme = gg_theme) +
    ggtitle(column_name) +
    labs(color = color_text) +
    theme(legend.position = legend.position)
}


plot_multi_violin_and_feature <- function(
    seu,
    features,                  
    assay = NULL,              
    layer = "data",            
    region_assay = NULL,       
    reduction = NULL,
    cluster_by = NULL          # <-- 新增：data.frame(barcode, cluster)
){
  stopifnot(inherits(seu, "Seurat"))
  suppressPackageStartupMessages({
    library(Seurat); library(Signac); library(Matrix)
    library(GenomicRanges); library(ggplot2); library(tidyr)
  })
  
  if (is.null(assay)) assay <- DefaultAssay(seu)
  if (is.null(reduction)) {
    rds <- names(Reductions(seu))
    reduction <- if ("umap" %in% rds) "umap" else if ("tsne" %in% rds) "tsne" else "pca"
  }
  
  # ====== 分簇来源：优先使用自定义映射 ======
  if (!is.null(cluster_by)) {
    clu_map <- setNames(cluster_by$cluster, cluster_by$barcode)
    clv <- unname(clu_map[colnames(seu)])
    df <- data.frame(cell = colnames(seu), cluster = clv, check.names = FALSE)
    # 让分簇水平按传入顺序固定，便于对照
    df$cluster <- factor(df$cluster, levels = unique(cluster_by$cluster))
  } else {
    df <- data.frame(cell = colnames(seu), cluster = as.factor(Idents(seu)), check.names = FALSE)
  }
  
  values_list <- list()
  
  # ====== 基因：来自指定 assay/layer ======
  gene_mat <- try(LayerData(seu, assay = assay, layer = layer), silent = TRUE)
  if (!inherits(gene_mat, "try-error")) {
    gene_feats <- intersect(features, rownames(gene_mat))
    for (g in gene_feats) values_list[[g]] <- as.numeric(gene_mat[g, ])
  }
  
  # ====== 基因组区域：在 ChromatinAssay 里按覆盖 peaks 求每细胞均值 ======
  is_region <- grepl("^chr[^:]+[:\\-]\\d+[-_]\\d+$", features, ignore.case = TRUE)
  region_vec <- features[is_region]
  if (length(region_vec)) {
    if (is.null(region_assay)) {
      cand <- intersect(c("ATAC", "peaks"), names(seu@assays))
      if (length(cand) == 0) {
        chrom_assays <- names(Filter(function(a) inherits(a, "ChromatinAssay"), seu@assays))
        if (length(chrom_assays)) cand <- chrom_assays
      }
      region_assay <- if (length(cand)) cand[1] else stop("未找到 ChromatinAssay（ATAC/peaks）。")
    }
    pm <- LayerData(seu, assay = region_assay, layer = "counts")
    rn <- rownames(pm)
    peaks_gr <- try(Signac::StringToGRanges(rn, sep = c("-", "-")), silent = TRUE)
    if (inherits(peaks_gr, "try-error")) peaks_gr <- Signac::StringToGRanges(rn, sep = c(":", "-"))
    
    for (reg in region_vec) {
      sep <- if (grepl(":", reg)) c(":", "-") else c("-", "-")
      reg_gr <- Signac::StringToGRanges(reg, sep = sep)
      hits <- GenomicRanges::findOverlaps(peaks_gr, reg_gr)
      idx  <- unique(S4Vectors::queryHits(hits))
      lab  <- paste0("region_", gsub("[:\\-]", "_", reg))
      if (length(idx)) {
        score <- Matrix::colMeans(pm[idx, , drop = FALSE])
        seu[[lab]] <- score
        values_list[[lab]] <- score
      } else {
        warning("该区域未覆盖任何 peak: ", reg)
        seu[[lab]] <- NA_real_
        values_list[[lab]] <- rep(NA_real_, ncol(seu))
      }
    }
  }
  
  # ====== 长表 -> 并排小提琴 ======
  for (nm in names(values_list)) df[[nm]] <- values_list[[nm]]
  df_long <- tidyr::pivot_longer(df, cols = all_of(names(values_list)),
                                 names_to = "feature", values_to = "value")
  df_long <- df_long[!is.na(df_long$cluster), ]  # 去掉未匹配的条形码
  
  vln <- ggplot(df_long, aes(x = cluster, y = value, fill = feature)) +
    geom_violin(position = position_dodge(width = 0.8),
                scale = "width", trim = TRUE, width = 0.7) +
    labs(x = NULL, y = paste0(layer, " (", assay, ")"), fill = "Feature") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # ====== FeaturePlot 组合图 ======
  old_assay <- DefaultAssay(seu)
  DefaultAssay(seu) <- assay
  fp <- FeaturePlot(seu, features = names(values_list), reduction = reduction, combine = TRUE)
  DefaultAssay(seu) <- old_assay
  
  list(violin = vln, feature_plot = fp)
}


plot_lineage_river <- function(dslt = dslt, embedding_name = "cGR_smoothed") {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggalluvial)
  
  rna_clust <- dslt[["columnMetadata"]][["RNA"]][["lineage"]] %>%
    select(Barcode, s) %>% rename(s_rna = s)
  
  cg_obj <- dslt[["embeddings"]][["lineage"]][[embedding_name]][["louvain_clusters"]]
  cg_clust <- if (is.data.frame(cg_obj)) cg_obj else data.frame(louvain = cg_obj)
  if (!"Barcode" %in% names(cg_clust)) cg_clust$Barcode <- rownames(cg_clust)
  names(cg_clust)[names(cg_clust) == "louvain_clusters"] <- "louvain"
  cg_clust <- cg_clust %>%
    select(Barcode, louvain) %>%
    mutate(louvain = as.factor(louvain))
  
  atac_ok <- FALSE
  atac_try <- try(dslt[["columnMetadata"]][["ATAC"]][["lineage"]], silent = TRUE)
  if (!inherits(atac_try, "try-error") && !is.null(atac_try) && NROW(atac_try) > 0) atac_ok <- TRUE
  if (atac_ok) {
    atac_clust <- atac_try %>% select(Barcode, s) %>% rename(s_atac = s)
  }
  
  if (atac_ok) {
    df <- rna_clust %>%
      inner_join(cg_clust,  by = "Barcode") %>%
      inner_join(atac_clust, by = "Barcode") %>%
      mutate(s_rna = as.factor(s_rna), s_atac = as.factor(s_atac))
    dat_long <- df %>%
      transmute(Barcode, `ATAC louvain` = s_atac,
                `Drug: louvain` = louvain, `RNA louvain` = s_rna) %>%
      pivot_longer(cols = c(`ATAC louvain`, `Drug: louvain`, `RNA louvain`),
                   names_to = "axis", values_to = "cluster") %>%
      mutate(freq = 1L,
             axis = factor(axis, levels = c("ATAC louvain", "Drug: louvain", "RNA louvain")))
  } else {
    df <- rna_clust %>%
      inner_join(cg_clust, by = "Barcode") %>%
      mutate(s_rna = as.factor(s_rna))
    dat_long <- df %>%
      transmute(Barcode, `Drug: louvain` = louvain, `RNA louvain` = s_rna) %>%
      pivot_longer(cols = c(`Drug: louvain`, `RNA louvain`),
                   names_to = "axis", values_to = "cluster") %>%
      mutate(freq = 1L,
             axis = factor(axis, levels = c("Drug: louvain", "RNA louvain")))
  }
  

  p <- ggplot(dat_long,
              aes(x = axis, stratum = cluster, alluvium = Barcode,
                  y = freq, fill = cluster, label = cluster)) +
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    geom_text(stat = "stratum", size = 4) +
    theme_minimal() +
    theme(legend.position = "none")
  
  return(p)
}


plot_rna_atac_hist <- function(dslt, seu_rna, seu_atac = NULL, binwidth = 1) {

  rna_map  <- dslt[["assays"]][["rna_map"]]
  atac_map <- dslt[["assays"]][["atac_map"]]
  
  keep_rna  <- rownames(seu_rna@meta.data)
  keep_atac <- if (!is.null(seu_atac)) rownames(seu_atac@meta.data) else keep_rna
  
  make_hist <- function(df, keep_ids, title) {
    df_f <- df %>%
      filter(Barcode %in% keep_ids) %>%
      distinct(Barcode, .keep_all = TRUE)
    
    ggplot(df_f, aes(x = n_cells)) +
      geom_histogram(binwidth = binwidth, fill = "steelblue", color = "white") +
      scale_y_continuous(trans = "log1p") +
      labs(x = "n_cells", y = "Count (log1p)", title = title) +
      theme_classic()
  }
  
  # 绘制 RNA 图
  p_rna <- make_hist(rna_map, keep_rna, "RNA map")
  
  # 判断 ATAC 是否为空
  if (is.null(atac_map) || (is.data.frame(atac_map) && nrow(atac_map) == 0)) {
    p_rna
  } else {
    p_atac <- make_hist(atac_map, keep_atac, "ATAC map")
    p_rna + p_atac
  }
}


compact_plot <- function(p) {
  if (is.null(p) || !inherits(p, "gg")) return(NULL)
  ggplotGrob(p)
}