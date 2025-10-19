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
    stop("Please upload matrix.mtx(.gz), features/peaks.tsv(.gz), and barcodes.tsv(.gz) files.")
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

get_dslt_embedding_safe <- function(dslt, assay_name, slot) {
  tryCatch(dslt$getEmbedding(assay_name, slot), error = function(e) NULL)
}

get_dslt_assay_safe <- function(dslt, level, assay_name) {
  tryCatch(dslt$getAssay(level, assay_name, force = TRUE), error = function(e) NULL)
}

get_dslt_column_metadata_safe <- function(dslt, assay_name, field) {
  tryCatch(dslt$getColumnMetadata(assay_name, field), error = function(e) NULL)
}

ensure_knn_embeddings <- function(state, assay_name) {
  if (is.null(state$dslt) || !nzchar(assay_name)) return(list(coords = NULL, clusters = NULL))
  assay_present <- !is.null(get_dslt_assay_safe(state$dslt, "lineage", assay_name))
  if (!assay_present) {
    return(list(coords = NULL, clusters = NULL))
  }
  coords <- get_dslt_embedding_safe(state$dslt, assay_name, "umap")
  clusters <- get_dslt_embedding_safe(state$dslt, assay_name, "louvain_clusters")
  if (is.null(coords) || is.null(clusters)) {
    dslt_updated <- tryCatch(
      runKnnAnalysisAndStore(state$dslt, assay_name),
      error = function(e) {
        showNotification(paste("KNN analysis failed:", e$message), type = "error")
        NULL
      }
    )
    if (is.null(dslt_updated)) {
      return(list(coords = NULL, clusters = NULL))
    }
    state$dslt <- dslt_updated
    coords <- get_dslt_embedding_safe(state$dslt, assay_name, "umap")
    clusters <- get_dslt_embedding_safe(state$dslt, assay_name, "louvain_clusters")
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

ensure_umap <- function(seu, dims, reduction_name = "umap") {
  dims <- seq(dims[1], dims[2])
  if (!reduction_name %in% names(Reductions(seu))) {
    seu <- RunUMAP(seu, dims = dims, reduction.name = reduction_name, verbose = FALSE)
  }
  list(seu = seu, reduction = reduction_name)
}

ensure_tsne <- function(seu, dims, reduction_name = "tsne") {
  dims <- seq(dims[1], dims[2])
  if (!reduction_name %in% names(Reductions(seu))) {
    seu <- RunTSNE(seu, dims = dims, reduction.name = reduction_name, verbose = FALSE)
  }
  list(seu = seu, reduction = reduction_name)
}

ensure_pca <- function(seu, npcs = 50) {
  if (!"pca" %in% names(Reductions(seu))) {
    seu <- RunPCA(seu, npcs = npcs, verbose = FALSE)
  }
  seu
}

ensure_clusters <- function(seu, res = 0.8, dims = 1:30) {
  if (!"RNA_snn" %in% names(seu@graphs)) {
    seu <- FindNeighbors(seu, dims = dims)
  }
  if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
    seu <- FindClusters(seu, resolution = res)
  }
  seu
}

ensure_kmeans_clusters <- function(seu, k, dims = 1:30) {
  k <- max(2, as.integer(k))
  seu <- ensure_pca(seu, max(dims))
  key <- paste0("kmeans", k)
  if (!key %in% colnames(seu@meta.data)) {
    seu <- run_kmeans_clustering(seu, knum = k, dimsl = dims[1], dimsh = dims[length(dims)])
  }
  list(seu = seu, column = key)
}

# -------------------------------------------------------------------------
# Server ------------------------------------------------------------------
# -------------------------------------------------------------------------

server <- function(input, output, session) {
  state <- reactiveValues(
    dslt = NULL,
    denoised_assays = character(),
    seurat = list(
      sc_rna_raw = NULL,
      sc_rna = NULL,
      sc_atac_raw = NULL,
      sc_atac = NULL,
      pb_rna = NULL,
      pb_atac = NULL
    ),
    mapping = list(rna = NULL, atac = NULL),
    drug_matrix = NULL,
    lineage_drug_values = NULL,
    single_drug_values = NULL,
    lineage_cell_counts = NULL,
    single_cell_counts = NULL,
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
    export_history = data.frame(timestamp = character(), action = character(), stringsAsFactors = FALSE)
  )

  # -----------------------------------------------------------------------
  # Initial upload status messages ---------------------------------------
  # -----------------------------------------------------------------------
  output$lineage_upload_status <- renderText("Awaiting lineage RDS upload.")
  output$single_cell_upload_status <- renderText({
    rna_ready <- !is.null(input$single_cell_rna_matrix_files) && nrow(input$single_cell_rna_matrix_files) >= 3
    atac_ready <- !is.null(input$single_cell_atac_matrix_files) && nrow(input$single_cell_atac_matrix_files) >= 3
    if (rna_ready && atac_ready) {
      "Single-cell RNA and ATAC matrices detected."
    } else if (rna_ready) {
      "Single-cell RNA matrix detected. Upload ATAC inputs if needed."
    } else if (atac_ready) {
      "Single-cell ATAC matrix detected. Upload RNA inputs if needed."
    } else {
      "Awaiting single-cell matrix uploads."
    }
  })
  output$single_cell_upload_rna_matrix_status <- renderText({
    if (!is.null(input$single_cell_rna_matrix_files) && nrow(input$single_cell_rna_matrix_files) >= 3) {
      "Single-cell RNA matrix files uploaded."
    } else {
      "Awaiting RNA matrix files (matrix.mtx, features.tsv, barcodes.tsv)."
    }
  })
  output$single_cell_upload_atac_matrix_status <- renderText({
    if (!is.null(input$single_cell_atac_matrix_files) && nrow(input$single_cell_atac_matrix_files) >= 3) {
      "Single-cell ATAC matrix files uploaded."
    } else {
      "Awaiting ATAC matrix files (matrix.mtx, peaks.tsv, barcodes.tsv)."
    }
  })
  output$single_cell_atac_metadata_status <- renderText({
    if (!is.null(input$single_cell_atac_metadata_file)) {
      "ATAC metadata file uploaded."
    } else {
      "Awaiting ATAC metadata file."
    }
  })
  output$single_cell_atac_fragments_status <- renderText({
    if (!is.null(input$single_cell_atac_fragments_file)) {
      "ATAC fragments file uploaded."
    } else {
      "Awaiting ATAC fragments file."
    }
  })
  output$lineage_rna_mapping_status <- renderText("Optional: upload RNA barcode mapping file.")
  output$single_cell_atac_mapping_status <- renderText("Optional: upload ATAC barcode mapping file.")

  observeEvent(input$lineage_rds_file, {
    req(input$lineage_rds_file)
    output$lineage_upload_status <- renderText("Lineage data ready for loading.")
  })

  observeEvent(input$lineage_rna_mapping_file, {
    req(input$lineage_rna_mapping_file)
    output$lineage_rna_mapping_status <- renderText("RNA barcode mapping uploaded.")
  })

  observeEvent(input$single_cell_atac_mapping_file, {
    req(input$single_cell_atac_mapping_file)
    output$single_cell_atac_mapping_status <- renderText("ATAC barcode mapping uploaded.")
  })

  # -----------------------------------------------------------------------
  # Load data button ------------------------------------------------------
  # -----------------------------------------------------------------------
  observeEvent(input$load_all_data, {
    req(input$lineage_rds_file)

    tryCatch({
      # Reset derived state
      state$qc_applied <- FALSE
      state$seurat <- modifyList(state$seurat, list(
        sc_rna = NULL,
        sc_atac = NULL,
        pb_rna = NULL,
        pb_atac = NULL
      ))
      state$qc <- modifyList(state$qc, list(
        rna_violin = NULL,
        rna_scatter = NULL,
        rna_elbow = NULL,
        atac_density = NULL,
        atac_fragment = NULL,
        atac_violin = NULL,
        atac_depth = NULL,
        integration = NULL
      ))

      # Load dslt from LMM output ---------------------------------------
      state$dslt <- initializeDsltFromLmm(input$lineage_rds_file$datapath)

      lineage_assays <- names(state$dslt[["assays"]][["lineage"]])
      updatePickerInput(session, "dslt_denoise_assays", choices = lineage_assays)
      updatePickerInput(session, "denoise_assays_choice", choices = lineage_assays)

      state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
      state$single_drug_values <- state$dslt[["assays"]][["single_cell"]]

      drug_choices <- extract_drug_choices(state$lineage_drug_values)
      updatePickerInput(session, "lineage_drug_select", choices = drug_choices)
      updatePickerInput(session, "single_drug_select", choices = drug_choices)

      # Optional RNA matrices ------------------------------------------
      if (!is.null(input$single_cell_rna_matrix_files) && nrow(input$single_cell_rna_matrix_files) > 0) {
        rna_counts <- read_10x_from_upload(input$single_cell_rna_matrix_files)
        if (!is.null(rna_counts)) {
          sc_rna <- create_rna_seurat(matrix_data = rna_counts)
          state$seurat$sc_rna_raw <- sc_rna
        }
      }

      # Optional ATAC matrices -----------------------------------------
      atac_matrix_provided <- !is.null(input$single_cell_atac_matrix_files) && nrow(input$single_cell_atac_matrix_files) > 0
      atac_metadata_provided <- !is.null(input$single_cell_atac_metadata_file)
      atac_fragments_provided <- !is.null(input$single_cell_atac_fragments_file)
      if (atac_matrix_provided || atac_metadata_provided || atac_fragments_provided) {
        if (!(atac_matrix_provided && atac_metadata_provided && atac_fragments_provided)) {
          stop("Please upload ATAC matrix files, metadata, and fragments together.")
        }
        atac_counts <- read_10x_from_upload(input$single_cell_atac_matrix_files)
        metadata_path <- preserve_uploaded_file(input$single_cell_atac_metadata_file)
        fragments_path <- preserve_uploaded_file(input$single_cell_atac_fragments_file)
        sc_atac <- create_atac_seurat(
          counts_data = atac_counts,
          metadata_file = metadata_path,
          fragments_file = fragments_path
        )
        state$seurat$sc_atac_raw <- sc_atac
      }

      # Mapping files ---------------------------------------------------
      if (!is.null(input$lineage_rna_mapping_file)) {
        state$mapping$rna <- read_mapping_file(input$lineage_rna_mapping_file$datapath)
      }
      if (!is.null(input$single_cell_atac_mapping_file)) {
        state$mapping$atac <- read_mapping_file(input$single_cell_atac_mapping_file$datapath)
      }

      # Additional drug matrix ----------------------------------------
      if (!is.null(input$drug_matrix_file)) {
        state$drug_matrix <- normalize_matrix(read_table_like(input$drug_matrix_file$datapath))
        state$dslt[["assays"]][["lineage"]][["external_drug"]] <- state$drug_matrix
        state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
      }

      # Apply denoising if requested -----------------------------------
      selected_assays <- unique(c(input$denoise_assays_choice, input$dslt_denoise_assays))
      selected_assays <- intersect(selected_assays, names(state$dslt[["assays"]][["lineage"]]))
      denoised_names <- character()
      for (assay in selected_assays) {
        res <- applyAdaptiveKernelDenoising(state$dslt, assay)
        state$dslt <- res$dslt
        denoised_names <- c(denoised_names, res$smoothed_assay_name)
      }
      if (!is.null(state$drug_matrix) && isTRUE(input$denoise_drug_matrix)) {
        res <- applyAdaptiveKernelDenoising(state$dslt, "external_drug")
        state$dslt <- res$dslt
        denoised_names <- c(denoised_names, res$smoothed_assay_name)
      }
      state$denoised_assays <- denoised_names

      # Update drug choices after denoising -----------------------------
      state$lineage_drug_values <- state$dslt[["assays"]][["lineage"]]
      state$single_drug_values <- state$dslt[["assays"]][["single_cell"]]
      drug_choices <- extract_drug_choices(state$lineage_drug_values)
      updatePickerInput(session, "lineage_drug_select", choices = drug_choices)
      updatePickerInput(session, "single_drug_select", choices = drug_choices)
      dataset_choices <- names(state$lineage_drug_values)
      if (length(dataset_choices)) {
        updateSelectInput(session, "lineage_rds_object_select", choices = dataset_choices, selected = dataset_choices[1])
      }
      single_dataset_choices <- names(state$single_drug_values)
      if (length(single_dataset_choices)) {
        updateSelectInput(session, "single_rds_object_select", choices = single_dataset_choices, selected = single_dataset_choices[1])
      }

      output$upload_summary <- renderText({
        lineage_status <- if (is.null(state$dslt)) "Not loaded" else "Loaded"
        single_cell_rna_status <- if (is.null(state$seurat$sc_rna_raw)) "Not uploaded" else "Uploaded"
        single_cell_atac_status <- if (is.null(state$seurat$sc_atac_raw)) "Not uploaded" else "Uploaded"
        mapping_status <- if (is.null(state$mapping$rna) && is.null(state$mapping$atac)) "Not provided" else "Provided"
        paste(
          "Lineage data:", lineage_status,
          "\nSingle-cell RNA:", single_cell_rna_status,
          "\nSingle-cell ATAC:", single_cell_atac_status,
          "\nMapping:", mapping_status,
          "\nDenoised assays:", ifelse(length(state$denoised_assays), paste(state$denoised_assays, collapse = ", "), "None")
        )
      })
    }, error = function(e) {
      showNotification(paste("Failed to load data:", e$message), type = "error")
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
    updateNumericInput(session, "atac_peak_fragments_min", value = 1000)
    updateNumericInput(session, "atac_peak_fragments_max", value = 100000)
    updateNumericInput(session, "atac_pct_reads_peaks_min", value = 40)
    updateNumericInput(session, "atac_blacklist_ratio_max", value = 5)
    updateNumericInput(session, "atac_nucleosome_signal_max", value = 2)
    updateNumericInput(session, "atac_tss_enrichment_min", value = 2)
    output$qc_settings_status <- renderText("QC settings reset to default values!")
  })

  observeEvent(input$apply_qc_settings, {
    validate(need(!is.null(state$dslt), "Load datasets before applying QC."))

    suffix <- input$barcode_suffix
    if (!nzchar(suffix)) suffix <- NULL

    # RNA QC -------------------------------------------------------------
    if (!is.null(state$seurat$sc_rna_raw)) {
      seu <- state$seurat$sc_rna_raw
      if (!is.null(suffix)) {
        seu <- tryCatch(filter_cells_by_suffix(seu, suffix = suffix), error = function(e) seu)
      }
      seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
      seu <- subset(seu,
        subset = nFeature_RNA > input$rna_nfeature_min &
          nFeature_RNA < input$rna_nfeature_max &
          nCount_RNA > input$rna_ncount_min &
          nCount_RNA < input$rna_ncount_max &
          percent.mt < input$rna_percent_mt_max)

      state$qc$rna_violin <- plot_rna_qc_violin(
        seu,
        min_features = input$rna_nfeature_min,
        max_features = input$rna_nfeature_max,
        max_percent_mt = input$rna_percent_mt_max
      )
      state$qc$rna_scatter <- plot_rna_qc_scatter(seu)
      seu <- ensure_pca(seu)
      state$qc$rna_elbow <- ElbowPlot(seu)

      seu <- ensure_clusters(seu)
      seu <- ensure_umap(seu, input$single_umap_pca_dims)$seu
      seu <- ensure_tsne(seu, input$single_tsne_pca_dims)$seu
      state$seurat$sc_rna <- seu
    }

    # ATAC QC ------------------------------------------------------------
    if (!is.null(state$seurat$sc_atac_raw)) {
      seu_atac <- state$seurat$sc_atac_raw
      if (!is.null(suffix)) {
        seu_atac <- tryCatch(filter_cells_by_suffix(seu_atac, suffix = suffix), error = function(e) seu_atac)
      }
      seu_atac <- calculate_nucleosome_signal(seu_atac)
      seu_atac <- calculate_tss_enrichment(seu_atac)
      seu_atac <- calculate_atac_qc_metrics(seu_atac)
      seu_atac <- filter_atac_cells(
        seu_atac,
        ncount_min = input$atac_peak_fragments_min,
        ncount_max = input$atac_peak_fragments_max,
        pct_reads_min = input$atac_pct_reads_peaks_min,
        blacklist_max = input$atac_blacklist_ratio_max,
        nucleosome_max = input$atac_nucleosome_signal_max,
        tss_min = input$atac_tss_enrichment_min,
        suffix = if (is.null(suffix)) "-1$" else suffix
      )
      state$qc$atac_density <- plot_tss_density_scatter(seu_atac)
      state$qc$atac_fragment <- plot_fragment_histogram(seu_atac)
      state$qc$atac_violin <- plot_atac_qc_violins(
        seu_atac,
        ncount_min = input$atac_peak_fragments_min,
        ncount_max = input$atac_peak_fragments_max,
        pct_reads_min = input$atac_pct_reads_peaks_min,
        tss_min = input$atac_tss_enrichment_min,
        nucleosome_max = input$atac_nucleosome_signal_max,
        blacklist_max = input$atac_blacklist_ratio_max / 100
      )
      state$seurat$sc_atac <- seu_atac
    }

    # Mapping to lineage level -----------------------------------------
    if (!is.null(state$mapping$rna) && !is.null(state$seurat$sc_rna)) {
      pb_rna <- lineage_map_seurat_rna(
        se = state$seurat$sc_rna,
        bc = state$mapping$rna,
        dslt = state$dslt,
        assay = "RNA",
        key = "rna"
      )
      pb_rna <- ensure_pca(pb_rna)
      pb_rna <- ensure_clusters(pb_rna)
      pb_rna <- ensure_umap(pb_rna, input$lineage_umap_pca_dims)$seu
      pb_rna <- ensure_tsne(pb_rna, input$lineage_tsne_pca_dims)$seu
      state$seurat$pb_rna <- pb_rna
      state$lineage_cell_counts <- state$dslt[["assays"]][["rna_map"]]
    }

    if (!is.null(state$mapping$atac) && !is.null(state$seurat$sc_atac)) {
      atac_counts <- LayerData(state$seurat$sc_atac, assay = "ATAC", layer = "counts")
      pb_atac <- CreateSeuratObject(counts = atac_counts)
      pb_atac <- ensure_pca(pb_atac)
      pb_atac <- ensure_clusters(pb_atac)
      pb_atac <- ensure_umap(pb_atac, input$lineage_umap_pca_dims)$seu
      pb_atac <- ensure_tsne(pb_atac, input$lineage_tsne_pca_dims)$seu
      state$seurat$pb_atac <- pb_atac
    }

    # Store single cell drug values if mapping available ----------------
    if (!is.null(state$mapping$rna)) {
      try({
        state$dslt <- map_lineage_to_single_cell(state$dslt, state$mapping$rna)
        state$single_drug_values <- state$dslt[["assays"]][["single_cell"]]
      }, silent = TRUE)
    }

    state$qc_applied <- TRUE
    output$qc_settings_status <- renderText("QC settings applied successfully!")
  })

  # -----------------------------------------------------------------------
  # QC plots --------------------------------------------------------------
  # -----------------------------------------------------------------------
  output$qc_rna_qc <- renderPlot({
    req(state$qc$rna_violin)
    state$qc$rna_violin
  })

  output$qc_rna_pca_elbow <- renderPlot({
    req(state$qc$rna_elbow)
    state$qc$rna_elbow
  })

  output$qc_rna_fea <- renderPlot({
    req(state$qc$rna_scatter)
    state$qc$rna_scatter
  })

  output$qc_atac_Density <- renderPlot({
    req(state$qc$atac_density)
    state$qc$atac_density
  })

  output$qc_atac_Fragment <- renderPlot({
    req(state$qc$atac_fragment)
    state$qc$atac_fragment
  })

  output$qc_atac_QC <- renderPlot({
    req(state$qc$atac_violin)
    state$qc$atac_violin
  })

  output$qc_atac_depth_correlation <- renderPlot({
    req(state$seurat$sc_atac)
    seu <- state$seurat$sc_atac
    ggplot(seu@meta.data, aes(x = nCount_peaks, y = peak_region_fragments)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "lm", se = FALSE, colour = "red") +
      theme_minimal() +
      labs(x = "nCount_peaks", y = "Peak region fragments")
  })

  output$integration_plot <- renderPlot({
    req(state$seurat$sc_rna, state$seurat$sc_atac)
    umap_rna <- Embeddings(state$seurat$sc_rna, "umap")
    umap_atac <- Embeddings(state$seurat$sc_atac, "umap")
    df_rna <- data.frame(umap_rna[, 1:2], assay = "RNA")
    df_atac <- data.frame(umap_atac[, 1:2], assay = "ATAC")
    colnames(df_rna) <- colnames(df_atac) <- c("dim1", "dim2", "assay")
    df <- rbind(df_rna, df_atac)
    ggplot(df, aes(x = dim1, y = dim2, colour = assay)) +
      geom_point(alpha = 0.4, size = 0.8) +
      theme_minimal() +
      coord_equal() +
      labs(x = "UMAP1", y = "UMAP2", colour = "Assay")
  })

  # -----------------------------------------------------------------------
  # Helper to obtain plotting reduction ----------------------------------
  # -----------------------------------------------------------------------
  lineage_plot_data <- reactive({
    req(state$qc_applied)
    seu <- state$seurat$pb_rna
    if (is.null(seu)) seu <- state$seurat$sc_rna
    req(seu)

    method <- input$lineage_clustering_method
    reduction <- "pca"
    cluster_col <- "seurat_clusters"
    dims <- seq(1, 30)

    if (method == "umap") {
      res <- ensure_umap(seu, input$lineage_umap_pca_dims)
      seu <- res$seu
      reduction <- res$reduction
    } else if (method == "tsne") {
      res <- ensure_tsne(seu, input$lineage_tsne_pca_dims)
      seu <- res$seu
      reduction <- res$reduction
    } else if (method == "kmeans") {
      res <- ensure_kmeans_clusters(seu, input$lineage_kmeans_input, dims = dims)
      seu <- res$seu
      cluster_col <- res$column
    } else {
      seu <- ensure_pca(seu)
      seu <- ensure_clusters(seu)
      reduction <- "pca"
    }

    state$seurat$pb_rna <- seu
    list(seu = seu, reduction = reduction, cluster_col = cluster_col)
  })

  single_plot_data <- reactive({
    req(state$qc_applied)
    seu <- state$seurat$sc_rna
    req(seu)

    method <- input$single_clustering_method
    reduction <- "pca"
    cluster_col <- "seurat_clusters"
    dims <- seq(1, 30)

    if (method == "umap") {
      res <- ensure_umap(seu, input$single_umap_pca_dims)
      seu <- res$seu
      reduction <- res$reduction
    } else if (method == "tsne") {
      res <- ensure_tsne(seu, input$single_tsne_pca_dims)
      seu <- res$seu
      reduction <- res$reduction
    } else if (method == "kmeans") {
      res <- ensure_kmeans_clusters(seu, input$single_kmeans_input, dims = dims)
      seu <- res$seu
      cluster_col <- res$column
    } else {
      seu <- ensure_pca(seu)
      seu <- ensure_clusters(seu)
      reduction <- "pca"
    }

    state$seurat$sc_rna <- seu
    list(seu = seu, reduction = reduction, cluster_col = cluster_col)
  })

  # -----------------------------------------------------------------------
  # Lineage plots --------------------------------------------------------
  # -----------------------------------------------------------------------
  lineage_cluster_plot_obj <- reactive({
    pd <- lineage_plot_data()
    seu <- pd$seu
    reduction <- pd$reduction
    cluster_col <- pd$cluster_col

    if (input$lineage_color_by == "drug_response") {
      scores <- get_drug_scores(
        state$lineage_drug_values,
        input$lineage_rds_object_select,
        input$lineage_drug_select
      )
      if (!is.null(scores)) {
        meta <- seu@meta.data
        meta <- dplyr::left_join(meta %>% tibble::rownames_to_column("barcode"), scores, by = "barcode")
        rownames(meta) <- meta$barcode
        seu@meta.data <- meta
        cluster_col <- "score"
      }
    }

    DimPlot(seu, reduction = reduction, group.by = cluster_col, pt.size = 1) +
      theme_minimal()
  })

  lineage_knn_plot_obj <- reactive({
    req(input$lineage_rds_object_select)
    res <- ensure_knn_embeddings(state, input$lineage_rds_object_select)
    validate(need(!is.null(res$coords) && !is.null(res$clusters), "Run KNN analysis on the selected assay."))
    coords <- as.matrix(res$coords)
    cluster_df <- as.data.frame(res$clusters)
    cluster_vals <- if ("louvain_clusters" %in% colnames(cluster_df)) cluster_df$louvain_clusters else cluster_df[[1]]
    cluster_vals <- as.factor(cluster_vals)
    ggscatter_colored(
      coords,
      cluster_vals,
      ggObj = ggplot(),
      dimnamesXYZ = c("UMAP1", "UMAP2", "Cluster"),
      size = 0.5,
      gg_theme = theme_umap_legend,
      colors = cluster_palette
    ) +
      guides(color = guide_legend(override.aes = list(size = 3))) +
      ggtitle("Phenomics clusters") +
      coord_fixed()
  })

  output$lineage_rna_cluster_plot <- renderPlot({ req(lineage_cluster_plot_obj()); lineage_cluster_plot_obj() })
  output$lineage_knn_plot <- renderPlot({ req(lineage_knn_plot_obj()); lineage_knn_plot_obj() })

  output$lineage_atac_cluster_plot <- renderPlot({
    req(state$seurat$pb_atac)
    seu <- state$seurat$pb_atac
    seu <- ensure_umap(seu, input$lineage_umap_pca_dims)$seu
    DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 1) + theme_minimal()
  })

  lineage_drug_response_plot_obj <- reactive({
    req(state$dslt, state$lineage_drug_values)
    assay_name <- input$lineage_rds_object_select
    req(assay_name)
    res <- ensure_knn_embeddings(state, assay_name)
    coords <- as.matrix(res$coords)
    validate(need(!is.null(coords), "Run KNN analysis on the selected assay."))
    mat <- get_dslt_assay_safe(state$dslt, "lineage", assay_name)
    req(!is.null(mat))
    mat <- as.matrix(mat)
    mat <- normalize_matrix(mat)
    selected_drugs <- intersect(colnames(mat), input$lineage_drug_select)
    validate(need(length(selected_drugs) > 0, "Select at least one drug."))
    if (is.null(rownames(coords)) && !is.null(rownames(mat))) {
      rownames(coords) <- rownames(mat)
    }
    plot_list <- lapply(selected_drugs, function(drug) {
      ggscatter_single(
        coords = coords,
        values = mat,
        column_name = drug,
        ggObj = ggplot() + coord_fixed(),
        size_mult = 0.2,
        colors = drug_palette,
        gg_theme = theme_umap,
        symmQuant = 0.95,
        legend.position = "right"
      )
    })
    patchwork::wrap_plots(plotlist = plot_list)
  })

  output$lineage_drug_response_plot <- renderPlot({ req(lineage_drug_response_plot_obj()); lineage_drug_response_plot_obj() })

  lineage_bubble_plot_obj <- reactive({
    req(state$dslt, input$lineage_rds_object_select, input$lineage_bubble_mode)
    assay_name <- input$lineage_rds_object_select
    if (!nzchar(assay_name)) return(NULL)
    cell_line <- assay_name
    if (input$lineage_bubble_mode == "cluster") {
      res <- ensure_knn_embeddings(state, assay_name)
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
        title = paste0(cell_line, " Cluster Analysis of cGR scores"),
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
      arch <- get_dslt_column_metadata_safe(state$dslt, assay_name, "archetypes")
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
        title = paste0(cell_line, " Archetypal Analysis of cGR scores"),
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

  output$lineage_violin_plot <- renderPlot({
    req(state$seurat$pb_rna)
    features <- if (input$lineage_track_mode == "gene") input$lineage_gene_input else input$lineage_atac_region_input
    if (!nzchar(features)) return(NULL)
    feat_vec <- str_split(features, "[,;\\s]+", simplify = TRUE)
    feat_vec <- feat_vec[feat_vec != ""]
    res <- plot_multi_violin_and_feature(
      seu = state$seurat$pb_rna,
      features = feat_vec,
      assay = "RNA",
      reduction = "umap"
    )
    combine_feature_plots(res)
  })

  output$lineage_coverage_plot <- renderPlot({
    req(state$seurat$pb_atac)
    features <- if (input$lineage_track_mode == "region") input$lineage_atac_region_input else input$lineage_gene_input
    if (!nzchar(features)) return(NULL)
    feat_vec <- str_split(features, "[,;\\s]+", simplify = TRUE)
    feat_vec <- feat_vec[feat_vec != ""]
    res <- plot_multi_violin_and_feature(
      seu = state$seurat$pb_atac,
      features = feat_vec,
      assay = "ATAC",
      reduction = "umap"
    )
    combine_feature_plots(res)
  })

  # -----------------------------------------------------------------------
  # Single-cell plots ----------------------------------------------------
  # -----------------------------------------------------------------------
  single_cluster_plot_obj <- reactive({
    pd <- single_plot_data()
    seu <- pd$seu
    reduction <- pd$reduction
    cluster_col <- pd$cluster_col

    if (input$single_color_by == "drug_response") {
      scores <- get_drug_scores(
        state$single_drug_values,
        input$single_rds_object_select,
        input$single_drug_select
      )
      if (!is.null(scores)) {
        meta <- seu@meta.data
        meta <- dplyr::left_join(meta %>% tibble::rownames_to_column("barcode"), scores, by = "barcode")
        rownames(meta) <- meta$barcode
        seu@meta.data <- meta
        cluster_col <- "score"
      }
    }

    DimPlot(seu, reduction = reduction, group.by = cluster_col, pt.size = 0.6) +
      theme_minimal()
  })

  single_knn_plot_obj <- reactive({
    assay_name <- input$single_rds_object_select
    res <- if (!is.null(assay_name)) ensure_knn_embeddings(state, assay_name) else list(coords = NULL, clusters = NULL)
    if (!is.null(res$coords) && !is.null(res$clusters)) {
      coords <- as.matrix(res$coords)
      cluster_df <- as.data.frame(res$clusters)
      cluster_vals <- if ("louvain_clusters" %in% colnames(cluster_df)) cluster_df$louvain_clusters else cluster_df[[1]]
      cluster_vals <- as.factor(cluster_vals)
      ggscatter_colored(
        coords,
        cluster_vals,
        ggObj = ggplot(),
        dimnamesXYZ = c("UMAP1", "UMAP2", "Cluster"),
        size = 0.4,
        gg_theme = theme_umap_legend,
        colors = cluster_palette
      ) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggtitle("Single-cell phenomics clusters") +
        coord_fixed()
    } else {
      pd <- single_plot_data()
      seu <- pd$seu
      reduction <- pd$reduction
      cluster_col <- pd$cluster_col
      coords <- Embeddings(seu, reduction)
      coords_df <- as.data.frame(coords)
      cluster_vals <- if (!is.null(cluster_col) && cluster_col %in% colnames(seu@meta.data)) {
        seu@meta.data[[cluster_col]]
      } else {
        Idents(seu)
      }
      ggscatter_colored(
        coords_df,
        as.factor(cluster_vals),
        ggObj = ggplot(),
        dimnamesXYZ = c(paste0(toupper(reduction), "1"), paste0(toupper(reduction), "2"), "Cluster"),
        size = 0.4,
        gg_theme = theme_umap_legend,
        colors = cluster_palette
      ) +
        guides(color = guide_legend(override.aes = list(size = 3))) +
        ggtitle("Single-cell clustering") +
        coord_fixed()
    }
  })

  output$single_rna_cluster_plot <- renderPlot({ req(single_cluster_plot_obj()); single_cluster_plot_obj() })
  output$single_knn_plot <- renderPlot({ req(single_knn_plot_obj()); single_knn_plot_obj() })

  output$single_atac_cluster_plot <- renderPlot({
    req(state$seurat$sc_atac)
    seu <- ensure_umap(state$seurat$sc_atac, input$single_umap_pca_dims)$seu
    DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", pt.size = 0.6) + theme_minimal()
  })

  single_drug_response_plot_obj <- reactive({
    req(input$single_rds_object_select, state$single_drug_values)
    assay_name <- input$single_rds_object_select
    coords <- NULL
    mat <- NULL
    if (!is.null(state$dslt)) {
      res <- ensure_knn_embeddings(state, assay_name)
      coords <- res$coords
      mat <- get_dslt_assay_safe(state$dslt, "single_cell", assay_name)
      if (is.null(mat)) {
        mat <- state$single_drug_values[[assay_name]]
      }
    }
    if (is.null(coords)) {
      pd <- single_plot_data()
      coords <- Embeddings(pd$seu, pd$reduction)
    }
    if (is.null(mat)) {
      mat <- state$single_drug_values[[assay_name]]
    }
    req(!is.null(coords), !is.null(mat))
    mat <- as.matrix(mat)
    mat <- normalize_matrix(mat)
    selected <- intersect(colnames(mat), input$single_drug_select)
    validate(need(length(selected) > 0, "Select at least one drug."))
    coords <- as.matrix(coords)
    if (is.null(rownames(coords)) && !is.null(rownames(mat))) {
      rownames(coords) <- rownames(mat)
    }
    plot_list <- lapply(selected, function(drug) {
      ggscatter_single(
        coords = coords,
        values = mat,
        column_name = drug,
        ggObj = ggplot() + coord_fixed(),
        size_mult = 0.2,
        colors = drug_palette,
        gg_theme = theme_umap,
        symmQuant = 0.95,
        legend.position = "right"
      )
    })
    patchwork::wrap_plots(plotlist = plot_list)
  })

  output$single_drug_response_plot <- renderPlot({ req(single_drug_response_plot_obj()); single_drug_response_plot_obj() })

  output$single_violin_plot <- renderPlot({
    req(state$seurat$sc_rna)
    features <- if (input$single_track_mode == "gene") input$single_gene_input else input$single_atac_region_input
    if (!nzchar(features)) return(NULL)
    feat_vec <- str_split(features, "[,;\\s]+", simplify = TRUE)
    feat_vec <- feat_vec[feat_vec != ""]
    pd <- single_plot_data()
    res <- plot_multi_violin_and_feature(
      seu = pd$seu,
      features = feat_vec,
      assay = "RNA",
      reduction = pd$reduction
    )
    combine_feature_plots(res)
  })

  output$single_coverage_plot <- renderPlot({
    req(state$seurat$sc_atac)
    features <- if (input$single_track_mode == "region") input$single_atac_region_input else input$single_gene_input
    if (!nzchar(features)) return(NULL)
    feat_vec <- str_split(features, "[,;\\s]+", simplify = TRUE)
    feat_vec <- feat_vec[feat_vec != ""]
    pd <- single_plot_data()
    res <- plot_multi_violin_and_feature(
      seu = state$seurat$sc_atac,
      features = feat_vec,
      assay = "ATAC",
      reduction = "umap"
    )
    combine_feature_plots(res)
  })

  # -----------------------------------------------------------------------
  # Export handlers -------------------------------------------------------
  # -----------------------------------------------------------------------
  output$lineage_data_available <- reactive({ !is.null(state$seurat$pb_rna) })
  output$single_cell_data_available <- reactive({ !is.null(state$seurat$sc_rna) })
  outputOptions(output, "lineage_data_available", suspendWhenHidden = FALSE)
  outputOptions(output, "single_cell_data_available", suspendWhenHidden = FALSE)

  output$download_lineage_data <- downloadHandler(
    filename = function() paste0("lineage_data_", Sys.Date(), ".rds"),
    content = function(file) {
      saveRDS(list(dslt = state$dslt, seurat = state$seurat$pb_rna), file)
      update_export_history(state, "Exported lineage data")
    }
  )

  output$download_single_cell_data <- downloadHandler(
    filename = function() paste0("single_cell_data_", Sys.Date(), ".rds"),
    content = function(file) {
      saveRDS(list(seurat = state$seurat$sc_rna, atac = state$seurat$sc_atac), file)
      update_export_history(state, "Exported single-cell data")
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
      req(state$seurat$sc_rna)
      df <- state$seurat$sc_rna@meta.data
      df$barcode <- rownames(df)
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
      plots <- list(
        isolate(lineage_cluster_plot_obj()),
        isolate(lineage_knn_plot_obj()),
        isolate(lineage_drug_response_plot_obj()),
        isolate(lineage_bubble_plot_obj())
      )
      purrr::walk(Filter(Negate(is.null), plots), print)
      dev.off()
      update_export_history(state, "Exported lineage plots")
    }
  )

  output$download_single_plots <- downloadHandler(
    filename = function() paste0("single_plots_", Sys.Date(), ".", input$export_single_format),
    content = function(file) {
      device <- switch(input$export_single_format, png = png, pdf = pdf, jpeg = jpeg)
      device(file, width = 10, height = 7)
      plots <- list(
        isolate(single_cluster_plot_obj()),
        isolate(single_knn_plot_obj()),
        isolate(single_drug_response_plot_obj())
      )
      purrr::walk(Filter(Negate(is.null), plots), print)
      dev.off()
      update_export_history(state, "Exported single-cell plots")
    }
  )

  output$download_custom_plot <- downloadHandler(
    filename = function() paste0(input$custom_plot_name, "_", Sys.Date(), ".", input$custom_plot_format),
    content = function(file) {
      device <- switch(input$custom_plot_format, png = png, pdf = pdf, jpeg = jpeg)
      device(file, width = 10, height = 7)
      plot_custom <- NULL
      plot_custom <- NULL
      if (input$custom_plot_type == "rna_cluster") {
        plot_custom <- if (input$custom_plot_level == "lineage") isolate(lineage_cluster_plot_obj()) else isolate(single_cluster_plot_obj())
      } else if (input$custom_plot_type == "drug_embedding") {
        plot_custom <- if (input$custom_plot_level == "lineage") isolate(lineage_drug_response_plot_obj()) else isolate(single_drug_response_plot_obj())
      } else if (input$custom_plot_type == "knn_embedding") {
        plot_custom <- if (input$custom_plot_level == "lineage") isolate(lineage_knn_plot_obj()) else isolate(single_knn_plot_obj())
      } else if (input$custom_plot_type == "bubble" && input$custom_plot_level == "lineage") {
        plot_custom <- isolate(lineage_bubble_plot_obj())
      } else if (input$custom_plot_type == "bubble") {
        showNotification("Bubble plots are only available for lineage level.", type = "warning")
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
    status <- if (is.null(state$seurat$pb_rna)) "yellow" else "green"
    message <- if (is.null(state$seurat$pb_rna)) "Lineage data not ready" else "Lineage data ready"
    infoBox("Lineage dataset", message, icon = icon("sitemap"), color = status)
  })

  output$single_export_status <- renderInfoBox({
    status <- if (is.null(state$seurat$sc_rna)) "yellow" else "green"
    message <- if (is.null(state$seurat$sc_rna)) "Single-cell data not ready" else "Single-cell data ready"
    infoBox("Single-cell dataset", message, icon = icon("braille"), color = status)
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
