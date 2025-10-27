library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(bslib)
library(Seurat)
library(Signac)
library(Matrix)
library(patchwork)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)
library(stringr)
library(ggplot2)
library(glue)
library(EnsDb.Hsapiens.v86)
library(SingleCellExperiment)
library(muscat)
source("fuctions/functions.R")

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

  observeEvent(input$load_supplements, {
    withProgress(message = "Loading supplement files...", value = 0, {
      tryCatch({
        incProgress(0.05, detail = "Initializing...")

        state$seurat$sc_rna_raw <- NULL
        state$seurat$sc_rna <- NULL
        state$seurat$pb_rna <- NULL
        state$seurat$sc_atac_raw <- NULL
        state$seurat$sc_atac <- NULL
        state$seurat$pb_atac <- NULL
        state$mapping <- list(rna = NULL, atac = NULL)
        state$metadata$atac <- NULL
        state$qc_applied <- FALSE

        incProgress(0.1, detail = "Checking RNA inputs...")

        if (has_uploaded_files(input$single_cell_rna_rds_file)) {
          incProgress(0.1, detail = "Reading RNA Seurat object...")
          sc_rna <- tryCatch(readRDS(input$single_cell_rna_rds_file$datapath), error = function(e) {
            showNotification(paste("Failed to read RNA RDS:", e$message), type = "error", duration = 5)
            stop(e$message)
          })
          if (!inherits(sc_rna, "Seurat")) {
            stop("The provided RNA RDS does not contain a Seurat object.")
          }
          state$seurat$sc_rna_raw <- sc_rna
          showNotification("RNA Seurat object loaded successfully", type = "message", duration = 3)
        } else if (has_uploaded_files(input$single_cell_rna_matrix_files)) {
          incProgress(0.1, detail = "Reading RNA matrix files...")
          rna_counts <- read_10x_from_upload(input$single_cell_rna_matrix_files)
          if (!is.null(rna_counts)) {
            incProgress(0.05, detail = "Creating RNA Seurat object...")
            sc_rna <- create_rna_seurat(matrix_data = rna_counts)
            state$seurat$sc_rna_raw <- sc_rna
            showNotification("RNA matrix loaded successfully", type = "message", duration = 3)
          }
        } else {
          incProgress(0.15, detail = "No RNA inputs provided.")
        }

        incProgress(0.1, detail = "Checking ATAC inputs...")

        atac_rds_provided <- has_uploaded_files(input$single_cell_atac_rds_file)
        atac_matrix_provided <- has_uploaded_files(input$single_cell_atac_matrix_files)
        atac_fragments_provided <- has_uploaded_files(input$single_cell_atac_fragments_file) &&
          any(grepl("\\.(tsv|tsv\\.gz)$", tolower(input$single_cell_atac_fragments_file$name)))

        if (atac_rds_provided) {
          incProgress(0.1, detail = "Reading ATAC Seurat object...")
          sc_atac <- tryCatch(readRDS(input$single_cell_atac_rds_file$datapath), error = function(e) {
            showNotification(paste("Failed to read ATAC RDS:", e$message), type = "error", duration = 5)
            stop(e$message)
          })
          if (!inherits(sc_atac, "Seurat")) {
            stop("The provided ATAC RDS does not contain a Seurat object.")
          }
          state$seurat$sc_atac_raw <- sc_atac
          showNotification("ATAC Seurat object loaded successfully", type = "message", duration = 3)
        } else if (atac_matrix_provided && atac_fragments_provided) {
          incProgress(0.1, detail = "Reading ATAC matrix files...")
          atac_counts <- read_10x_from_upload(input$single_cell_atac_matrix_files)

          incProgress(0.1, detail = "Processing fragments file...")
          fragments_info <- preserve_fragments_with_index(input$single_cell_atac_fragments_file)
          fragments_path <- fragments_info$fragments
          if (is.null(fragments_path)) {
            stop("Fragments file missing from upload.")
          }
          ensure_tabix_index(fragments_path)

          atac_metadata <- NULL
          if (has_uploaded_files(input$single_cell_atac_metadata_file)) {
            incProgress(0.05, detail = "Reading ATAC metadata...")
            atac_metadata <- tryCatch({
              read.csv(
                file = input$single_cell_atac_metadata_file$datapath,
                header = TRUE,
                row.names = 1,
                check.names = FALSE,
                stringsAsFactors = FALSE
              )
            }, error = function(e) {
              showNotification(paste("Failed to read ATAC metadata:", e$message), type = "error", duration = 5)
              NULL
            })
            if (!is.null(atac_metadata)) {
              atac_metadata <- as.data.frame(atac_metadata, stringsAsFactors = FALSE)
              state$metadata$atac <- atac_metadata
              showNotification("ATAC metadata loaded successfully", type = "message", duration = 3)
            }
          }

          incProgress(0.1, detail = "Creating ATAC Seurat object...")
          sc_atac <- create_atac_seurat(
            counts_data = atac_counts,
            fragments_file = fragments_path,
            metadata = state$metadata$atac
          )
          state$seurat$sc_atac_raw <- sc_atac
          showNotification("ATAC data loaded successfully", type = "message", duration = 3)
        } else if (atac_matrix_provided || atac_fragments_provided) {
          incProgress(0.2, detail = "ATAC files incomplete, skipping...")
          showNotification("Both ATAC matrix and fragments files are required for ATAC processing.", type = "warning", duration = 5)
        } else {
          incProgress(0.2, detail = "No ATAC inputs provided.")
        }

        incProgress(0.05, detail = "Reading mapping files...")

        if (has_uploaded_files(input$lineage_rna_mapping_file)) {
          state$mapping$rna <- read_mapping_file(input$lineage_rna_mapping_file$datapath)
          showNotification("RNA mapping file loaded", type = "message", duration = 2)
        }
        if (has_uploaded_files(input$single_cell_atac_mapping_file)) {
          state$mapping$atac <- read_mapping_file(input$single_cell_atac_mapping_file$datapath)
          showNotification("ATAC mapping file loaded", type = "message", duration = 2)
        }

        incProgress(0.05, detail = "Processing external drug matrix...")

        if (has_uploaded_files(input$drug_matrix_file)) {
          validate(need(!is.null(state$dslt), "Load lineage data before integrating external drug matrix."))
          state$drug_matrix <- normalize_matrix(read_table_like(input$drug_matrix_file$datapath))
          state$dslt[["assays"]][["lineage"]][["external_drug"]] <- state$drug_matrix
          state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
          drug_matrix_text("External drug matrix integrated.")
          updatePickerInput(session, "lineage_drug_select", choices = extract_drug_choices(state$lineage_drug_values))
          showNotification("External drug matrix integrated", type = "message", duration = 3)
        }

        incProgress(0.05, detail = "Finalizing...")

        if (!is.null(state$seurat$sc_rna_raw) || !is.null(state$seurat$sc_atac_raw)) {
          single_cell_status("Single-cell inputs loaded. Apply QC to generate pseudo-bulk lineage data.")
        } else {
          single_cell_status("No single-cell inputs loaded.")
        }

        output$upload_summary <- renderText(compose_upload_summary(state))
        showNotification("Supplement files loaded successfully!", type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Failed to load supplements:", e$message), type = "error", duration = 8)
      })
    })
  })

  observeEvent(input$clear_single_cell, {
    state$seurat$sc_rna_raw <- NULL
    state$seurat$sc_rna <- NULL
    state$seurat$sc_atac_raw <- NULL
    state$seurat$sc_atac <- NULL
    state$mapping <- list(rna = NULL, atac = NULL)
    state$metadata$atac <- NULL
    single_cell_status("Single-cell inputs cleared from memory.")
    output$upload_summary <- renderText(compose_upload_summary(state))
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
  needs_reduction <- !reduction %in% names(Reductions(seu))
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

# -------------------------------------------------------------------------
# Server ------------------------------------------------------------------
# -------------------------------------------------------------------------

server <- function(input, output, session) {
  state <- reactiveValues(
    dslt = NULL,
    denoised_assays = character(),
    drug_matrix = NULL,
    lineage_drug_values = NULL,
    qc = list(
      rna_violin = NULL,
      rna_scatter = NULL,
      rna_elbow = NULL,
      atac_density = NULL,
      atac_fragment = NULL,
      atac_violin = NULL,
      atac_depth = NULL,
      integration = NULL
    ),
    qc_applied = FALSE,
    export_history = data.frame(timestamp = character(), action = character(), stringsAsFactors = FALSE),
    seurat = list(
      sc_rna_raw = NULL,
      sc_rna = NULL,
      sc_atac_raw = NULL,
      sc_atac = NULL,
      pb_rna = NULL,
      pb_atac = NULL
    ),
    mapping = list(rna = NULL, atac = NULL),
    metadata = list(rna = NULL, atac = NULL),
    single_cell_counts = NULL,
    lineage_cell_counts = NULL
  )

  # -----------------------------------------------------------------------
  # Initial upload status messages ---------------------------------------
  # -----------------------------------------------------------------------
  output$lineage_upload_status <- renderText("Awaiting lineage RDS upload.")
  has_uploaded_files <- function(info) {
    !is.null(info) && nrow(info) > 0
  }

  single_cell_status <- reactiveVal("Awaiting single-cell uploads.")
  rna_matrix_status <- reactiveVal("Upload RNA matrix files (matrix.mtx, features.tsv, barcodes.tsv) or provide a preprocessed RDS object.")
  rna_rds_status <- reactiveVal("Optional: upload a preprocessed single-cell RNA Seurat object (.rds).")
  atac_matrix_status <- reactiveVal("Upload ATAC matrix files (matrix.mtx, peaks.tsv, barcodes.tsv) or provide a preprocessed RDS object.")
  atac_rds_status <- reactiveVal("Optional: upload a preprocessed single-cell ATAC Seurat object (.rds).")
  atac_fragments_status <- reactiveVal("Awaiting ATAC fragments file.")
  atac_metadata_status <- reactiveVal("Optional: upload ATAC metadata CSV with cell barcodes as row names.")
  rna_mapping_status <- reactiveVal("Optional: upload RNA barcode mapping file.")
  atac_mapping_status <- reactiveVal("Optional: upload ATAC barcode mapping file.")
  drug_matrix_text <- reactiveVal("Optional: upload external drug matrix.")

  output$single_cell_upload_status <- renderText(single_cell_status())
  output$single_cell_upload_rna_matrix_status <- renderText(rna_matrix_status())
  output$single_cell_rna_rds_status <- renderText(rna_rds_status())
  output$single_cell_upload_atac_matrix_status <- renderText(atac_matrix_status())
  output$single_cell_atac_rds_status <- renderText(atac_rds_status())
  output$single_cell_atac_fragments_status <- renderText(atac_fragments_status())
  output$single_cell_atac_metadata_status <- renderText(atac_metadata_status())
  output$lineage_rna_mapping_status <- renderText(rna_mapping_status())
  output$single_cell_atac_mapping_status <- renderText(atac_mapping_status())
  output$drug_matrix_status <- renderText(drug_matrix_text())
  output$qc_snapshot_status <- renderText("Upload a QC snapshot to restore a previous session.")

  observe({
    info <- input$single_cell_rna_matrix_files
    if (has_uploaded_files(info)) {
      rna_matrix_status("Single-cell RNA matrix files uploaded.")
    } else {
      rna_matrix_status("Upload RNA matrix files (matrix.mtx, features.tsv, barcodes.tsv) or provide a preprocessed RDS object.")
    }
  })

  observe({
    info <- input$single_cell_rna_rds_file
    if (has_uploaded_files(info)) {
      rna_rds_status("Preprocessed single-cell RNA Seurat object uploaded.")
    } else {
      rna_rds_status("Optional: upload a preprocessed single-cell RNA Seurat object (.rds).")
    }
  })

  observe({
    info <- input$single_cell_atac_matrix_files
    if (has_uploaded_files(info)) {
      atac_matrix_status("Single-cell ATAC matrix files uploaded.")
    } else {
      atac_matrix_status("Upload ATAC matrix files (matrix.mtx, peaks.tsv, barcodes.tsv) or provide a preprocessed RDS object.")
    }
  })

  observe({
    info <- input$single_cell_atac_rds_file
    if (has_uploaded_files(info)) {
      atac_rds_status("Preprocessed single-cell ATAC Seurat object uploaded.")
    } else {
      atac_rds_status("Optional: upload a preprocessed single-cell ATAC Seurat object (.rds).")
    }
  })

  observe({
    info <- input$single_cell_atac_fragments_file
    if (has_uploaded_files(info)) {
      lower_names <- tolower(info$name)
      if (any(grepl("\\.tbi(\\.gz)?$", lower_names))) {
        atac_fragments_status("ATAC fragments and index files uploaded.")
      } else {
        atac_fragments_status("ATAC fragments file uploaded. A matching .tbi index will be generated when possible.")
      }
    } else {
      atac_fragments_status("Awaiting ATAC fragments file.")
    }
  })

  observe({
    info <- input$single_cell_atac_metadata_file
    if (has_uploaded_files(info)) {
      atac_metadata_status("ATAC metadata CSV uploaded.")
    } else {
      atac_metadata_status("Optional: upload ATAC metadata CSV with cell barcodes as row names.")
    }
  })

  observe({
    info <- input$lineage_rna_mapping_file
    if (has_uploaded_files(info)) {
      rna_mapping_status("RNA barcode mapping uploaded.")
    } else {
      rna_mapping_status("Optional: upload RNA barcode mapping file.")
    }
  })

  observe({
    info <- input$single_cell_atac_mapping_file
    if (has_uploaded_files(info)) {
      atac_mapping_status("ATAC barcode mapping uploaded.")
    } else {
      atac_mapping_status("Optional: upload ATAC barcode mapping file.")
    }
  })

  observe({
    if (!is.null(state$seurat$pb_rna) || !is.null(state$seurat$pb_atac)) {
      return()
    }
    rna_ready <- has_uploaded_files(input$single_cell_rna_matrix_files) ||
      has_uploaded_files(input$single_cell_rna_rds_file)
    atac_ready <- has_uploaded_files(input$single_cell_atac_matrix_files) ||
      has_uploaded_files(input$single_cell_atac_rds_file)

    if (rna_ready && atac_ready) {
      single_cell_status("Single-cell RNA and ATAC inputs detected.")
    } else if (rna_ready) {
      single_cell_status("Single-cell RNA inputs detected. Upload ATAC inputs if needed.")
    } else if (atac_ready) {
      single_cell_status("Single-cell ATAC inputs detected. Upload RNA inputs if needed.")
    } else {
      single_cell_status("Awaiting single-cell uploads.")
    }
  })

  observe({
    if (has_uploaded_files(input$drug_matrix_file)) {
      drug_matrix_text("External drug matrix uploaded.")
    } else if (is.null(state$drug_matrix)) {
      drug_matrix_text("Optional: upload external drug matrix.")
    }
  })

  observeEvent(input$lineage_rds_file, {
    req(input$lineage_rds_file)
    output$lineage_upload_status <- renderText("Lineage data ready for loading.")
  })

  observeEvent(input$qc_snapshot_file, {
    req(input$qc_snapshot_file)
    output$qc_snapshot_status <- renderText("QC snapshot ready to load.")
  })

  # -----------------------------------------------------------------------
  # Load data buttons (split into 3 separate actions) -------------------
  # -----------------------------------------------------------------------
  
  # Button 1: Load Lineage RDS and update denoise choices
  observeEvent(input$load_lineage_rds, {
    req(input$lineage_rds_file)
    
    withProgress(message = "Loading lineage data...", value = 0, {
      tryCatch({
        incProgress(0.3, detail = "Reading RDS file")
        
        # Load dslt from LMM output
        state$dslt <- initializeDsltFromLmm(input$lineage_rds_file$datapath)
        
        incProgress(0.4, detail = "Updating UI choices")
        
        # Update denoise assay choices
        lineage_assays <- names(state$dslt[["assays"]][["lineage"]])
        updatePickerInput(session, "denoise_assays_choice", choices = lineage_assays)
        
        # Store lineage drug values
        state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
        
        # Update drug selection
        drug_choices <- extract_drug_choices(state$lineage_drug_values)
        updatePickerInput(session, "lineage_drug_select", choices = drug_choices)
        
        # Update dataset selection
        dataset_choices <- names(state$lineage_drug_values)
        if (length(dataset_choices)) {
          previous_selection <- isolate(input$lineage_rds_object_select)
          default_selection <- if (!is.null(previous_selection) && previous_selection %in% dataset_choices) {
            previous_selection
          } else {
            dataset_choices[1]
          }
          updatePickerInput(session, "lineage_rds_object_select",
            choices = dataset_choices,
            selected = default_selection
          )
        } else {
          updatePickerInput(session, "lineage_rds_object_select", choices = character(0))
        }
        
        incProgress(0.3, detail = "Complete")
        
        showNotification("Lineage RDS file loaded successfully!", type = "message")
        
        output$upload_summary <- renderText(compose_upload_summary(state))
        
      }, error = function(e) {
        showNotification(paste("Failed to load lineage RDS:", e$message), type = "error")
      })
    })
  })
  
    observeEvent(input$load_qc_snapshot, {
    req(input$qc_snapshot_file)

    withProgress(message = "Loading QC snapshot...", value = 0, {
      tryCatch({
        incProgress(0.3, detail = "Reading snapshot file")
        snapshot <- readRDS(input$qc_snapshot_file$datapath)

        validate(need(is.list(snapshot), "Snapshot file must contain a list object."))
        validate(need(!is.null(snapshot$dslt), "Snapshot is missing 'dslt' data."))

        incProgress(0.6, detail = "Restoring lineage data")
        state$dslt <- snapshot$dslt
        if (!is.null(snapshot$denoised_assays)) {
          state$denoised_assays <- snapshot$denoised_assays
        } else {
          state$denoised_assays <- character()
        }
        if (!is.null(snapshot$drug_matrix)) {
          state$drug_matrix <- snapshot$drug_matrix
          drug_matrix_text("External drug matrix integrated.")
        }
        if (!is.null(snapshot$pb_rna)) {
          state$seurat$pb_rna <- snapshot$pb_rna
        }
        if (!is.null(snapshot$pb_atac)) {
          state$seurat$pb_atac <- snapshot$pb_atac
        }
        state$seurat$sc_rna_raw <- NULL
        state$seurat$sc_rna <- NULL
        state$seurat$sc_atac_raw <- NULL
        state$seurat$sc_atac <- NULL
        if (!is.null(snapshot$mapping) && is.list(snapshot$mapping)) {
          state$mapping <- modifyList(state$mapping, snapshot$mapping)
        }
        if (!is.null(state$dslt[["assays"]][["lineage"]])) {
          state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
        } else if (!is.null(snapshot$lineage_drug_values)) {
          state$lineage_drug_values <- snapshot$lineage_drug_values
        }
        if (!is.null(snapshot$qc) && is.list(snapshot$qc)) {
          state$qc <- modifyList(state$qc, snapshot$qc)
        }

        incProgress(0.8, detail = "Refreshing selections")
        lineage_assays <- names(state$lineage_drug_values)
        updatePickerInput(session, "denoise_assays_choice",
          choices = lineage_assays,
          selected = intersect(state$denoised_assays, lineage_assays)
        )
        updatePickerInput(session, "lineage_drug_select", choices = extract_drug_choices(state$lineage_drug_values))
        updatePickerInput(session, "lineage_rds_object_select",
          choices = lineage_assays,
          selected = if (length(lineage_assays)) lineage_assays[1] else character(0)
        )

        state$qc_applied <- TRUE
        if (!is.null(state$seurat$pb_rna) || !is.null(state$seurat$pb_atac)) {
          single_cell_status("Pseudo-bulk data restored from snapshot.")
        } else {
          single_cell_status("Awaiting single-cell uploads.")
        }
        output$lineage_upload_status <- renderText("Lineage data restored from QC snapshot.")
        output$upload_summary <- renderText(compose_upload_summary(state))
        showNotification("QC snapshot loaded successfully!", type = "message", duration = 5)
      }, error = function(e) {
        output$qc_snapshot_status <- renderText(paste("Failed to load snapshot:", e$message))
        showNotification(paste("Failed to load QC snapshot:", e$message), type = "error", duration = 8)
      })
    })
  })

  # Button 3: Apply denoising to selected assays
  observeEvent(input$apply_denoise, {
    req(state$dslt)
    
    selected_assays <- input$denoise_assays_choice
    if (is.null(selected_assays) || length(selected_assays) == 0) {
      showNotification("Please select at least one assay to denoise.", type = "warning")
      return()
    }
    
    withProgress(message = "Applying denoising...", value = 0, {
      tryCatch({
        selected_assays <- intersect(selected_assays, names(state$dslt[["assays"]][["lineage"]]))
        
        if (length(selected_assays) == 0) {
          showNotification("Selected assays not found in lineage data.", type = "error")
          return()
        }
        
        denoised_names <- character()
        total <- length(selected_assays) + if (!is.null(state$drug_matrix)) 1 else 0
        progress_step <- 1 / total
        
        # Denoise each selected assay
        for (assay in selected_assays) {
          incProgress(progress_step, detail = paste("Denoising", assay))
          res <- applyAdaptiveKernelDenoising(state$dslt, assay)
          state$dslt <- res$dslt
          denoised_names <- c(denoised_names, res$smoothed_assay_name)
        }
        
        # Denoise external drug matrix if present
        if (!is.null(state$drug_matrix)) {
          incProgress(progress_step, detail = "Denoising external drug matrix")
          res <- applyAdaptiveKernelDenoising(state$dslt, "external_drug")
          state$dslt <- res$dslt
          denoised_names <- c(denoised_names, res$smoothed_assay_name)
        }
        
        state$denoised_assays <- denoised_names
        
        # Update drug values and choices
        state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
        drug_choices <- extract_drug_choices(state$lineage_drug_values)
        updatePickerInput(session, "lineage_drug_select", choices = drug_choices)
        
        # Update dataset choices
        dataset_choices <- names(state$lineage_drug_values)
        if (length(dataset_choices)) {
          default_selection <- tail(denoised_names, 1)
          if (!length(default_selection)) {
            default_selection <- isolate(input$lineage_rds_object_select)
          }
          if (!length(default_selection) || !(default_selection %in% dataset_choices)) {
            default_selection <- dataset_choices[length(dataset_choices)]
          }
          updatePickerInput(session, "lineage_rds_object_select",
            choices = dataset_choices,
            selected = default_selection
          )
        }
        
        showNotification(
          paste("Denoising completed! Smoothed assays:", paste(denoised_names, collapse = ", ")),
          type = "message",
          duration = 5
        )
        
        output$upload_summary <- renderText(compose_upload_summary(state))
        
      }, error = function(e) {
        showNotification(paste("Denoising failed:", e$message), type = "error")
      })
    })
  })

  # -----------------------------------------------------------------------
  # QC Settings handling --------------------------------------------------
  # -----------------------------------------------------------------------
  observeEvent(input$reset_qc_settings, {
    updateNumericInput(session, "rna_nfeature_min", value = 200)
    updateNumericInput(session, "rna_nfeature_max", value = 5000)
    updateNumericInput(session, "rna_ncount_min", value = 500)
    updateNumericInput(session, "rna_ncount_max", value = 30000)
    updateNumericInput(session, "rna_percent_mt_max", value = 20)
    updateNumericInput(session, "atac_peak_fragments_min", value = 2000)
    updateNumericInput(session, "atac_peak_fragments_max", value = 30000)
    updateNumericInput(session, "atac_pct_reads_peaks_min", value = 20)
    updateNumericInput(session, "atac_blacklist_ratio_max", value = 0.02)
    updateNumericInput(session, "atac_nucleosome_signal_max", value = 4)
    updateNumericInput(session, "atac_tss_enrichment_min", value = 2)
    output$qc_settings_status <- renderText("QC settings reset to default values!")
  })

  observeEvent(input$apply_qc_settings, {
    withProgress(message = "Applying QC...", value = 0, {
      validate(need(!is.null(state$dslt), "Load datasets before applying QC."))

      suffix <- input$barcode_suffix
      if (!nzchar(suffix)) suffix <- NULL

      lineage_rna_mapped <- FALSE
      lineage_atac_mapped <- FALSE
      had_single_inputs <- !is.null(state$seurat$sc_rna_raw) || !is.null(state$seurat$sc_atac_raw)

      if (!is.null(state$seurat$sc_rna_raw)) {
        incProgress(0.05, detail = "Filtering RNA cells by suffix...")
        seu <- state$seurat$sc_rna_raw
        if (!is.null(suffix)) {
          seu <- tryCatch(filter_cells_by_suffix(seu, suffix = suffix), error = function(e) seu)
        }

        incProgress(0.1, detail = "Generating RNA QC plots...")
        state$qc$rna_violin <- plot_rna_qc_violin(
          seu,
          min_features = input$rna_nfeature_min,
          max_features = input$rna_nfeature_max,
          max_percent_mt = input$rna_percent_mt_max
        )
        state$qc$rna_scatter <- plot_rna_qc_scatter(seu)

        incProgress(0.15, detail = "Filtering RNA cells by QC thresholds...")
        seu <- filter_rna_qc(
          seu,
          min_features = input$rna_nfeature_min,
          max_features = input$rna_nfeature_max,
          min_counts = input$rna_ncount_min,
          max_counts = input$rna_ncount_max,
          max_percent_mt = input$rna_percent_mt_max,
          verbose = FALSE
        )

        incProgress(0.2, detail = "Normalizing single-cell RNA data...")
        seu <- SCTransform(seu, verbose = FALSE)
        seu <- RunPCA(seu, npcs = 30, verbose = FALSE)
        state$qc$rna_elbow <- ElbowPlot(seu)
        state$seurat$sc_rna <- seu
      }

      if (!is.null(state$seurat$sc_atac_raw)) {
        incProgress(0.25, detail = "Filtering ATAC cells by suffix...")
        seu_atac <- state$seurat$sc_atac_raw
        if (!is.null(suffix)) {
          seu_atac <- tryCatch(filter_cells_by_suffix(seu_atac, suffix = suffix), error = function(e) seu_atac)
        }

        incProgress(0.3, detail = "Computing ATAC QC metrics...")
        seu_atac <- calculate_nucleosome_signal(seu_atac)
        if (!"TSS_fragments" %in% colnames(seu_atac@meta.data)) {
          seu_atac <- calculate_tss_enrichment(seu_atac)
        }
        seu_atac <- calculate_atac_qc_metrics(seu_atac)

        state$qc$atac_density <- plot_tss_density_scatter(seu_atac)
        state$qc$atac_fragment <- plot_fragment_histogram(seu_atac)
        state$qc$atac_violin <- plot_atac_qc_violins(
          seu_atac,
          ncount_min = input$atac_peak_fragments_min,
          ncount_max = input$atac_peak_fragments_max,
          pct_reads_min = input$atac_pct_reads_peaks_min,
          tss_min = input$atac_tss_enrichment_min,
          nucleosome_max = input$atac_nucleosome_signal_max,
          blacklist_max = input$atac_blacklist_ratio_max
        )

        incProgress(0.35, detail = "Filtering ATAC cells by QC thresholds...")
        seu_atac <- filter_atac_cells(
          seu_atac,
          ncount_min = input$atac_peak_fragments_min,
          ncount_max = input$atac_peak_fragments_max,
          pct_reads_min = input$atac_pct_reads_peaks_min,
          blacklist_max = input$atac_blacklist_ratio_max,
          nucleosome_max = input$atac_nucleosome_signal_max,
          tss_min = input$atac_tss_enrichment_min,
          verbose = FALSE
        )
        seu_atac <- RunTFIDF(seu_atac, verbose = FALSE)
        seu_atac <- FindTopFeatures(seu_atac, min.cutoff = "q0", verbose = FALSE)
        seu_atac <- RunSVD(seu_atac, verbose = FALSE)
        state$seurat$sc_atac <- seu_atac
      }

      incProgress(0.45, detail = "Mapping to lineage level...")

      if (!is.null(state$mapping$rna) && !is.null(state$seurat$sc_rna)) {
        pb_rna <- lineage_map_seurat_rna(
          se = state$seurat$sc_rna,
          bc = state$mapping$rna,
          dslt = state$dslt,
          assay = "RNA",
          key = "rna"
        )
        lineage_dims <- c(1, 30)
        pca_res <- ensure_pca(pb_rna, dslt = state$dslt, npcs = 30, assays = "RNA", level = "lineage")
        pb_rna <- pca_res$seu
        state$dslt <- pca_res$dslt

        kmeans_res <- ensure_kmeans_clusters(
          pb_rna,
          k = 5,
          res = 0.5,
          dims = 1:30,
          dslt = state$dslt,
          assays = "RNA",
          level = "lineage"
        )
        pb_rna <- kmeans_res$seu
        state$dslt <- kmeans_res$dslt

        umap_res <- ensure_umap(pb_rna, lineage_dims, dslt = state$dslt, assays = "RNA", level = "lineage")
        pb_rna <- umap_res$seu
        state$dslt <- umap_res$dslt

        tsne_res <- ensure_tsne(pb_rna, lineage_dims, dslt = state$dslt, assays = "RNA", level = "lineage")
        pb_rna <- tsne_res$seu
        state$dslt <- tsne_res$dslt

        state$seurat$pb_rna <- pb_rna
        lineage_rna_mapped <- TRUE
      }

      if (!is.null(state$mapping$atac) && !is.null(state$seurat$sc_atac)) {
        pb_atac <- lineage_map_seurat_atac(
          se = state$seurat$sc_atac,
          bc = state$mapping$atac,
          dslt = state$dslt,
          assay = "peaks",
          key = "atac"
        )
        atac_dims <- c(1, 30)
        pca_res <- ensure_pca(pb_atac, dslt = state$dslt, npcs = 30, assays = "ATAC", level = "lineage")
        pb_atac <- pca_res$seu
        state$dslt <- pca_res$dslt

        kmeans_res <- ensure_kmeans_clusters(
          pb_atac,
          k = 5,
          res = 0.5,
          dims = 2:30,
          dslt = state$dslt,
          assays = "ATAC",
          level = "lineage"
        )
        pb_atac <- kmeans_res$seu
        state$dslt <- kmeans_res$dslt

        umap_res <- ensure_umap(pb_atac, atac_dims, dslt = state$dslt, assays = "ATAC", level = "lineage")
        pb_atac <- umap_res$seu
        state$dslt <- umap_res$dslt

        tsne_res <- ensure_tsne(pb_atac, atac_dims, dslt = state$dslt, assays = "ATAC", level = "lineage")
        pb_atac <- tsne_res$seu
        state$dslt <- tsne_res$dslt

        state$seurat$pb_atac <- pb_atac
        lineage_atac_mapped <- TRUE
      }

      state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
      state$qc_applied <- TRUE

      removed_inputs <- c()
      if (lineage_rna_mapped) {
        state$seurat$sc_rna <- NULL
        state$seurat$sc_rna_raw <- NULL
        removed_inputs <- c(removed_inputs, "RNA")
      }
      if (lineage_atac_mapped) {
        state$seurat$sc_atac <- NULL
        state$seurat$sc_atac_raw <- NULL
        state$metadata$atac <- NULL
        removed_inputs <- c(removed_inputs, "ATAC")
      }

      if (length(removed_inputs)) {
        single_cell_status(glue::glue(
          'Pseudo-bulk data generated; single-cell {paste(removed_inputs, collapse = " and ")} objects removed from memory.'
        ))
      } else if (had_single_inputs) {
        single_cell_status("Single-cell inputs processed without mapping; provide barcode mapping to create pseudo-bulk data.")
      }

      if (length(state$denoised_assays)) {
        for (assay in state$denoised_assays) {
          if (!nzchar(assay)) next
          tryCatch({
            ensure_knn_embeddings(state, assay, level = "lineage")
            ensure_archetype_metadata(state, assay, level = "lineage")
          }, error = function(e) {
            showNotification(paste("Failed to refresh analysis for", assay, ":", e$message), type = "warning")
          })
        }
      }

      updatePickerInput(session, "lineage_drug_select", choices = extract_drug_choices(state$lineage_drug_values))
      updatePickerInput(session, "lineage_rds_object_select",
        choices = names(state$lineage_drug_values),
        selected = if (length(state$lineage_drug_values)) names(state$lineage_drug_values)[1] else character(0)
      )

      output$qc_settings_status <- renderText({
        status <- c()
        if (lineage_rna_mapped) status <- c(status, "RNA pseudo-bulk ready")
        if (lineage_atac_mapped) status <- c(status, "ATAC pseudo-bulk ready")
        if (!length(status)) {
          "Lineage QC applied without single-cell mapping."
        } else {
          paste("QC complete:", paste(status, collapse = ", "))
        }
      })

      output$upload_summary <- renderText(compose_upload_summary(state))

      tryCatch({
        snapshot_data <- list(
          dslt = state$dslt,
          denoised_assays = state$denoised_assays,
          drug_matrix = state$drug_matrix,
          lineage_drug_values = state$lineage_drug_values,
          qc = state$qc,
          pb_rna = state$seurat$pb_rna,
          pb_atac = state$seurat$pb_atac
        )
        filename <- paste0("qc_snapshot_", format(Sys.time(), "%Y-%m-%d_%H%M%S"), ".rds")
        saveRDS(snapshot_data, file = filename)
        showNotification(paste("QC snapshot saved to server at:", filename), type = "message", duration = 5)
      }, error = function(e) {
        showNotification(paste("Failed to save QC snapshot:", e$message), type = "error", duration = 10)
      })
    })
  })
  # -----------------------------------------------------------------------
  # QC plots --------------------------------------------------------------
  # -----------------------------------------------------------------------
  output$qc_rna_qc <- renderPlot({
    if (is.null(state$qc$rna_violin)) return(NULL)
    state$qc$rna_violin
  })

  output$qc_rna_pca_elbow <- renderPlot({
    if (is.null(state$qc$rna_elbow)) return(NULL)
    state$qc$rna_elbow
  })

  output$qc_rna_fea <- renderPlot({
    if (is.null(state$qc$rna_scatter)) return(NULL)
    state$qc$rna_scatter
  })

  output$qc_atac_Density <- renderPlot({
    if (is.null(state$qc$atac_density)) return(NULL)
    state$qc$atac_density
  })

  output$qc_atac_Fragment <- renderPlot({
    if (is.null(state$qc$atac_fragment)) return(NULL)
    state$qc$atac_fragment
  })

  output$qc_atac_QC <- renderPlot({
    if (is.null(state$qc$atac_violin)) return(NULL)
    state$qc$atac_violin
  })

  # -----------------------------------------------------------------------
  # Helper to obtain plotting reduction ----------------------------------
  # -----------------------------------------------------------------------

  lineage_plot_settings <- eventReactive(
    input$lineage_refresh_plots,
    {
      list(
        reduction = input$lineage_red_method,
        clustering = input$lineage_clustering_method,
        kmeans_k = suppressWarnings(as.integer(input$lineage_kmeans_input)),
        umap_pca_dims = input$lineage_umap_pca_dims,
        umap_svd_dims = input$lineage_umap_svd_dims,
        tsne_pca_dims = input$lineage_tsne_pca_dims,
        tsne_svd_dims = input$lineage_tsne_svd_dims,
        dataset = input$lineage_rds_object_select,
        selected_drugs = input$lineage_drug_select,
        bubble_mode = input$lineage_bubble_mode,
        track_mode = input$lineage_track_mode,
        gene_input = input$lineage_gene_input,
        region_input = input$lineage_atac_region_input
      )
    },
    ignoreNULL = FALSE
  )

  lineage_cluster_plot_obj <- reactive({
    req(state$qc_applied)
    settings <- lineage_plot_settings()
    seu <- state$seurat$pb_rna
    req(seu)

    reduction <- settings$reduction
    if (is.null(reduction) || !nzchar(reduction)) reduction <- "umap"
    cluster_col <- settings$clustering
    if (is.null(cluster_col) || !nzchar(cluster_col)) cluster_col <- "louvain"

    if (identical(reduction, "umap")) {
      dims <- settings$umap_pca_dims
      if (is.null(dims)) dims <- c(1, 30)
      res <- ensure_umap(seu, dims, dslt = state$dslt, assays = "RNA", level = "lineage")
      seu <- res$seu
      state$dslt <- res$dslt
      reduction <- res$reduction
    } else if (identical(reduction, "tsne")) {
      dims <- settings$tsne_pca_dims
      if (is.null(dims)) dims <- c(1, 30)
      res <- ensure_tsne(seu, dims, dslt = state$dslt, assays = "RNA", level = "lineage")
      seu <- res$seu
      state$dslt <- res$dslt
      reduction <- res$reduction
    } else {
      res <- ensure_pca(seu, dslt = state$dslt, assays = "RNA", level = "lineage")
      seu <- res$seu
      state$dslt <- res$dslt
      reduction <- res$reduction
    }

    if (identical(cluster_col, "kmeans")) {
      k <- settings$kmeans_k
      if (is.na(k) || k < 1) k <- 5L
      res <- ensure_kmeans_clusters(seu, k = k, dims = seq(1, 30), dslt = state$dslt, assays = "RNA", level = "lineage")
      seu <- res$seu
      state$dslt <- res$dslt
      cluster_col <- res$column
    } else {
      res <- ensure_clusters(seu, dslt = state$dslt, assays = "RNA", level = "lineage", dims = seq(1, 30))
      seu <- res$seu
      state$dslt <- res$dslt
      cluster_col <- "seurat_clusters"
    }

    DimPlot(
      seu,
      cols = cluster_palette,
      reduction = reduction,
      group.by = cluster_col
    ) +
      theme_minimal()
  })

  lineage_knn_plot_obj <- reactive({
    settings <- lineage_plot_settings()
    assay_name <- settings$dataset
    req(state$qc_applied, assay_name)
    res <- ensure_knn_embeddings(state, assay_name, level = "lineage")
    validate(need(!is.null(res$coords) && !is.null(res$clusters), "Run KNN analysis on the selected assay."))
    coords <- as.matrix(res$coords)
    validate(need(ncol(coords) >= 2, "Embedding requires at least two dimensions."))
    df <- as.data.frame(coords)
    names(df)[1:2] <- c("Dim1", "Dim2")
    df$cluster <- as.factor(res$clusters)
    ggplot(df, aes(x = Dim1, y = Dim2, color = cluster)) +
      geom_point(size = 1) +
      scale_color_manual(values = cluster_palette) +
      theme_umap_legend +
      coord_fixed()
  })

  output$lineage_rna_cluster_plot <- renderPlot({ req(lineage_cluster_plot_obj()); lineage_cluster_plot_obj() })
  output$lineage_knn_plot <- renderPlot({ req(lineage_knn_plot_obj()); lineage_knn_plot_obj() })

  lineage_atac_cluster_plot_obj <- reactive({
    settings <- lineage_plot_settings()
    req(state$seurat$pb_atac)
    seu <- state$seurat$pb_atac

    dims <- settings$umap_svd_dims
    if (is.null(dims)) dims <- c(2, 30)
    res <- ensure_umap(seu, dims, dslt = state$dslt, assays = "ATAC", level = "lineage")
    seu <- res$seu
    state$dslt <- res$dslt

    DimPlot(seu, reduction = res$reduction, group.by = "seurat_clusters", pt.size = 1) + theme_minimal()
  })

  output$lineage_atac_cluster_plot <- renderPlot({ req(lineage_atac_cluster_plot_obj()); lineage_atac_cluster_plot_obj() })

  lineage_drug_response_plot_obj <- reactive({
    settings <- lineage_plot_settings()
    req(state$dslt, state$lineage_drug_values)
    assay_name <- settings$dataset
    req(assay_name)
    res <- ensure_knn_embeddings(state, assay_name, level = "lineage")
    coords <- as.matrix(res$coords)
    validate(need(!is.null(coords), "Run KNN analysis on the selected assay."))
    mat <- get_dslt_assay_safe(state$dslt, "lineage", assay_name)
    req(!is.null(mat))
    mat <- as.matrix(mat)
    mat <- normalize_matrix(mat)
    selected_drugs <- intersect(colnames(mat), settings$selected_drugs)
    validate(need(length(selected_drugs) > 0, "Select at least one drug."))
    if (is.null(rownames(coords)) && !is.null(rownames(mat))) {
      rownames(coords) <- rownames(mat)
    }

    plot_list <- ggscatter_single(
      coords = coords,
      values = mat,
      column_name = selected_drugs,
      ggObj = ggplot() + coord_fixed(),
      size_mult = 0.2,
      colors = drug_palette,
      gg_theme = theme_umap,
      symmQuant = 0.95,
      legend.position = "right"
    )
    plot_list
  })

  output$lineage_drug_response_plot <- renderPlot({ req(lineage_drug_response_plot_obj()); lineage_drug_response_plot_obj() })

  lineage_bubble_plot_obj <- reactive({
    settings <- lineage_plot_settings()
    req(state$dslt, settings$dataset, settings$bubble_mode)
    assay_name <- settings$dataset
    if (!nzchar(assay_name)) return(NULL)
    if (settings$bubble_mode == "cluster") {
      res <- ensure_knn_embeddings(state, assay_name, level = "lineage")
      mat <- get_dslt_assay_safe(state$dslt, "lineage", assay_name)
      if (is.null(res$clusters) || is.null(mat)) return(NULL)
      cluster_df <- as.data.frame(res$clusters)
      cluster_vals <- if ("louvain_clusters" %in% colnames(cluster_df)) cluster_df$louvain_clusters else cluster_df[[1]]
      cluster_vals <- as.factor(cluster_vals)
      mat_df <- as.data.frame(as.matrix(mat))
      if (is.null(rownames(cluster_df)) && !is.null(rownames(mat_df))) {
        rownames(cluster_df) <- rownames(mat_df)
      }
      if (is.null(rownames(mat_df)) && !is.null(rownames(cluster_df))) {
        rownames(mat_df) <- rownames(cluster_df)
      }
      if (is.null(names(cluster_vals))) {
        names(cluster_vals) <- rownames(cluster_df)
      }
      if (is.null(names(cluster_vals)) && !is.null(rownames(mat_df))) {
        names(cluster_vals) <- rownames(mat_df)
      }
      summarized_clusters <- summarize_columns(mat_df, cluster_vals, order_rows = TRUE)
      summarized_mat <- as.matrix(summarized_clusters)
      ggshape_heatmap(
        summarized_mat,
        abs(summarized_mat),
        size_range = c(1, 7),
        theme_choice = ggplot2::theme_minimal() + theme(plot.title = element_text(size = 16, hjust = 0, vjust = 1)),
        shape_values = 21,
        value_label = "Sensitivity",
        size_label = "Effect size",
        row_label = "Cluster",
        column_label = "Treatment",
        cluster_rows = TRUE,
        colorscheme = rev(RdBl_mod3),
        symmQuant = 0.95,
        grid.pars = list(
          grid.size = 1,
          axis.text.x = element_text(size = 8, angle = 0, hjust = -1),
          grid.color = "#f4f4f4",
          grid.linetype = "solid"
        ),
        cluster_cols = TRUE,
        text.angle.x = 90
      ) + coord_fixed()
    } else {
      arch <- ensure_archetype_metadata(state, assay_name, level = "lineage")
      validate(need(!is.null(arch), "Archetype metadata not available."))
      arch_mat <- as.matrix(arch)
      values_mat <- t(tanh(arch_mat / 3))
      size_mat <- t(abs(tanh(arch_mat / 3)))
      ggshape_heatmap(
        values_mat,
        data_sizes = size_mat,
        theme_choice = ggplot2::theme_minimal() + theme(plot.title = element_text(size = 18, hjust = 0, vjust = 1)),
        shape_values = 21,
        size_range = c(0.5, 4),
        value_label = "Sensitivity",
        size_label = "Effect size",
        row_label = "Archetype",
        column_label = "Treatment",
        cluster_rows = TRUE,
        colorscheme = rev(RdBl_mod3),
        symmQuant = 0.95,
        grid.pars = list(
          grid.size = 0,
          grid.color = "#f4f4f4",
          grid.linetype = "solid"
        ),
        cluster_cols = TRUE,
        text.angle.x = 90
      ) + coord_fixed()
    }
  })

  output$lineage_bubble_plot <- renderPlot({ req(lineage_bubble_plot_obj()); lineage_bubble_plot_obj() })

  lineage_violin_plot_obj <- reactive({
    req(state$seurat$pb_rna)
    settings <- lineage_plot_settings()
    track_mode <- settings$track_mode
    if (is.null(track_mode) || !nzchar(track_mode)) track_mode <- "gene"
    features <- if (identical(track_mode, "gene")) settings$gene_input else settings$region_input
    if (is.null(features) || !nzchar(features)) return(NULL)
    feat_vec <- str_split(features, "[,;\\s]+", simplify = TRUE)
    feat_vec <- feat_vec[feat_vec != ""]
    res <- plot_multi_violin_and_feature(
      seu = state$seurat$pb_rna,
      features = feat_vec,
      assay = "SCT",
      reduction = "umap"
    )
    combine_feature_plots(res)
  })

  output$lineage_violin_plot <- renderPlot({
    plot_obj <- lineage_violin_plot_obj()
    if (is.null(plot_obj)) return(NULL)
    plot_obj
  })
  # -----------------------------------------------------------------------
  # Export handlers -------------------------------------------------------
  # -----------------------------------------------------------------------
  output$lineage_data_available <- reactive({ !is.null(state$dslt) })
  outputOptions(output, "lineage_data_available", suspendWhenHidden = FALSE)

  output$download_lineage_data <- downloadHandler(
    filename = function() paste0("lineage_data_", Sys.Date(), ".rds"),
    content = function(file) {
      saveRDS(list(
        dslt = state$dslt,
        denoised_assays = state$denoised_assays,
        drug_matrix = state$drug_matrix,
        pb_rna = state$seurat$pb_rna,
        pb_atac = state$seurat$pb_atac
      ), file)
      update_export_history(state, "Exported lineage data")
    }
  )

  output$download_qc_summary <- downloadHandler(
    filename = function() paste0("qc_summary_", Sys.Date(), ".csv"),
    content = function(file) {
      qc_data <- data.frame(
        Parameter = c(
          "nFeature_RNA_min", "nFeature_RNA_max", "nCount_RNA_min", "nCount_RNA_max",
          "percent.mt_max", "nCount_peaks_min", "nCount_peaks_max",
          "pct_reads_in_peaks_min", "blacklist_ratio_max", "nucleosome_signal_max", "TSS.enrichment_min"
        ),
        Value = c(
          input$rna_nfeature_min, input$rna_nfeature_max,
          input$rna_ncount_min, input$rna_ncount_max,
          input$rna_percent_mt_max,
          input$atac_peak_fragments_min, input$atac_peak_fragments_max,
          input$atac_pct_reads_peaks_min, input$atac_blacklist_ratio_max,
          input$atac_nucleosome_signal_max, input$atac_tss_enrichment_min
        )
      )
      readr::write_csv(qc_data, file)
      update_export_history(state, "Exported QC report")
    }
  )

  output$download_cluster_results <- downloadHandler(
    filename = function() paste0("cluster_results_", Sys.Date(), ".csv"),
    content = function(file) {
      settings <- lineage_plot_settings()
      assay_name <- settings$dataset
      req(assay_name)
      res <- ensure_knn_embeddings(state, assay_name, level = "lineage")
      clusters <- res$clusters
      df <- data.frame(
        barcode = if (!is.null(names(clusters))) names(clusters) else seq_along(clusters),
        cluster = clusters,
        stringsAsFactors = FALSE
      )
      readr::write_csv(df, file)
      update_export_history(state, "Exported clustering results")
    }
  )

  output$download_qc_plots <- downloadHandler(
    filename = function() paste0("qc_plots_", Sys.Date(), ".", input$export_qc_format),
    content = function(file) {
      device <- switch(input$export_qc_format, png = png, pdf = pdf, jpeg = jpeg)
      device(file, width = 10, height = 7)
      if (!is.null(state$qc$rna_violin)) print(state$qc$rna_violin)
      if (!is.null(state$qc$rna_scatter)) print(state$qc$rna_scatter)
      if (!is.null(state$qc$rna_elbow)) print(state$qc$rna_elbow)
      if (!is.null(state$qc$atac_violin)) print(state$qc$atac_violin)
      if (!is.null(state$qc$atac_density)) print(state$qc$atac_density)
      if (!is.null(state$qc$atac_fragment)) print(state$qc$atac_fragment)
      dev.off()
      update_export_history(state, "Exported QC plots")
    }
  )

  output$download_lineage_plots <- downloadHandler(
    filename = function() paste0("lineage_plots_", Sys.Date(), ".", input$export_lineage_format),
    content = function(file) {
      device <- switch(input$export_lineage_format, png = png, pdf = pdf, jpeg = jpeg)
      device(file, width = 10, height = 7)

      safe_eval <- function(expr) {
        tryCatch(expr, error = function(e) NULL)
      }

      plots <- list(
        safe_eval(isolate(lineage_cluster_plot_obj())),
        safe_eval(isolate(lineage_knn_plot_obj())),
        safe_eval(isolate(lineage_atac_cluster_plot_obj())),
        safe_eval(isolate(lineage_drug_response_plot_obj())),
        safe_eval(isolate(lineage_bubble_plot_obj())),
        safe_eval(isolate(lineage_violin_plot_obj()))
      )
      purrr::walk(Filter(Negate(is.null), plots), print)
      dev.off()
      update_export_history(state, "Exported lineage plots")
    }
  )

  output$download_custom_plot <- downloadHandler(
    filename = function() paste0(input$custom_plot_name, "_", Sys.Date(), ".", input$custom_plot_format),
    content = function(file) {
      device <- switch(input$custom_plot_format, png = png, pdf = pdf, jpeg = jpeg)
      device(file, width = 10, height = 7)
      plot_custom <- NULL
      safe_eval <- function(expr) {
        tryCatch(expr, error = function(e) NULL)
      }

      if (input$custom_plot_type == "rna_cluster") {
        plot_custom <- safe_eval(isolate(lineage_cluster_plot_obj()))
      } else if (input$custom_plot_type == "atac_cluster") {
        plot_custom <- safe_eval(isolate(lineage_atac_cluster_plot_obj()))
      } else if (input$custom_plot_type == "drug_embedding") {
        plot_custom <- safe_eval(isolate(lineage_drug_response_plot_obj()))
      } else if (input$custom_plot_type == "knn_embedding") {
        plot_custom <- safe_eval(isolate(lineage_knn_plot_obj()))
      } else if (input$custom_plot_type == "bubble") {
        plot_custom <- safe_eval(isolate(lineage_bubble_plot_obj()))
      } else if (input$custom_plot_type == "violin") {
        plot_custom <- safe_eval(isolate(lineage_violin_plot_obj()))
      }
      if (!is.null(plot_custom)) print(plot_custom)
      dev.off()
      update_export_history(state, glue::glue("Exported custom plot - {input$custom_plot_name}"))
    }
  )
  # -----------------------------------------------------------------------
  # Export info boxes -----------------------------------------------------
  # -----------------------------------------------------------------------
  output$lineage_export_status <- renderInfoBox({
    status <- if (is.null(state$dslt)) "yellow" else "green"
    message <- if (is.null(state$dslt)) "Lineage data not ready" else "Lineage data ready"
    infoBox("Lineage dataset", message, icon = icon("sitemap"), color = status)
  })

  output$qc_settings_info <- renderInfoBox({
    status <- if (state$qc_applied) "green" else "yellow"
    message <- if (state$qc_applied) "QC thresholds applied" else "Apply QC settings to enable exports"
    infoBox("QC settings", message, icon = icon("cog"), color = status)
  })

  output$custom_export_info <- renderInfoBox({
    available <- nrow(state$export_history) > 0
    status <- if (available) "blue" else "yellow"
    message <- if (available) "Recent exports listed below" else "No export actions yet"
    infoBox("Export log", message, icon = icon("clipboard-list"), color = status)
  })

  # -----------------------------------------------------------------------
  # Export history table --------------------------------------------------
  # -----------------------------------------------------------------------
  output$export_history_table <- renderTable({
    if (!nrow(state$export_history)) {
      return(data.frame(timestamp = "-", action = "No export records yet."))
    }
    state$export_history[order(state$export_history$timestamp, decreasing = TRUE), ]
  })

  observeEvent(input$clear_export_history, {
    showModal(modalDialog(
      title = "Clear export history",
      "This will remove all previous export entries. Continue?",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_clear_history", "Clear history", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirm_clear_history, {
    state$export_history <- state$export_history[0, ]
    removeModal()
  })
}
