library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(bslib)

# UI Definition
ui <- dashboardPage(
  skin = "blue",

  # 1. Dashboard Header
  dashboardHeader(title = "Multi-omic Response Explorer"),

  # 2. Sidebar Menu
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("QC Settings", tabName = "qc_settings", icon = icon("cog")),
      menuItem("QC", tabName = "qc", icon = icon("check-circle")),
      # -- Main Menu --
      menuItem("Main", tabName = "main", icon = icon("dna"), startExpanded = TRUE,
        menuSubItem("Lineage Level", tabName = "lineage_level", icon = icon("sitemap")),
        menuSubItem("Single Cell Level", tabName = "single_cell_level", icon = icon("dot-circle"))
      ),
      menuItem("File Export", tabName = "export", icon = icon("download"))
    )
  ),

  # 3. Dashboard Body Content
  dashboardBody(
    tags$head(
      tags$style(HTML(
        "
        .control-spacer {
          margin-top: 15px;
        }

        .action-row .btn {
          margin-right: 10px;
        }

        .export-history-table {
          margin-top: 10px;
        }

        .upload-section {
          border-left: 3px solid #3c8dbc;
          padding-left: 15px;
          margin-bottom: 20px;
        }

        .upload-section h4 {
          color: #3c8dbc;
          font-weight: 600;
          margin-bottom: 10px;
        }

        .status-text {
          padding: 8px 12px;
          background-color: #f4f4f4;
          border-radius: 4px;
          margin-top: 8px;
          font-size: 13px;
        }

        .action-buttons {
          display: flex;
          gap: 10px;
          margin-top: 15px;
        }

        .action-buttons .btn {
          flex: 1;
        }

        .upload-summary-box {
          background-color: #f9f9f9;
          border: 1px solid #ddd;
          border-radius: 4px;
          padding: 15px;
          margin-top: 20px;
        }

        .upload-summary-box h4 {
          margin-top: 0;
          color: #333;
          font-size: 16px;
        }
      "
      ))
    ),
    tabItems(
      # =================================================================
      # 3.0 Data Upload Tab
      # =================================================================
      tabItem(tabName = "upload",
        h2(icon("upload"), "Data Upload"),
        p("Upload and configure your multi-omic experiment data. Required files must be loaded first, followed by optional supplements."),
        
        fluidRow(
          # Left column: Required Files
          column(width = 6,
            box(
              title = tagList(icon("file-alt"), "Required Files"),
              width = NULL,
              solidHeader = TRUE,
              status = "primary",
              
              div(class = "upload-section",
                h4(icon("dna"), "Lineage Level Data"),
                fileInput(
                  "lineage_rds_file",
                  "Lineage drug-response RDS file",
                  accept = c(".rds", ".RDS"),
                  buttonLabel = "Browse...",
                  placeholder = "No file selected"
                ),
                helpText("Upload the processed lineage-level Seurat/SingleCellExperiment object containing drug response data."),
                div(class = "status-text", textOutput("lineage_upload_status"))
              ),
              
              hr(),
              
              div(class = "upload-section",
                h4(icon("sliders-h"), "Processing Options"),
                textInput(
                  "barcode_suffix", 
                  tagList(icon("filter"), "Barcode suffix filter"),
                  value = "-1$",
                  placeholder = "Regex pattern (e.g., -1$)"
                ),
                helpText("Filter cell barcodes using regex pattern. Default removes '-1' suffix."),
                
                pickerInput(
                  "denoise_assays_choice",
                  tagList(icon("broom"), "Select assays to denoise"),
                  choices = NULL,
                  multiple = TRUE,
                  options = pickerOptions(
                    actionsBox = TRUE, 
                    noneSelectedText = "Load lineage data first",
                    selectAllText = "Select All",
                    deselectAllText = "Deselect All"
                  )
                ),
                helpText("Apply adaptive kernel denoising to selected assays to reduce noise.")
              ),
              
              hr(),
              
              div(class = "action-buttons",
                actionButton("load_lineage_rds", tagList(icon("file-import"), "Load Lineage"), class = "btn-primary"),
                actionButton("apply_denoise", tagList(icon("magic"), "Apply Denoising"), class = "btn-warning")
              )
            )
          ),
          
          # Right column: Optional Supplements
          column(width = 6,
            box(
              title = tagList(icon("plus-circle"), "Optional Supplements"),
              width = NULL,
              solidHeader = TRUE,
              status = "info",
              
              div(class = "upload-section",
                h4(icon("dna"), "Single-cell RNA Matrix"),
                fileInput(
                  "single_cell_rna_matrix_files",
                  "RNA matrix files (10x format)",
                  accept = c(".mtx", ".mtx.gz", ".tsv", ".tsv.gz", ".h5", ".hdf5"),
                  multiple = TRUE,
                  buttonLabel = "Browse...",
                  placeholder = "No files selected"
                ),
                helpText("Upload matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz), or a single 10x HDF5 file (.h5/.hdf5)."),
                div(class = "status-text", textOutput("single_cell_upload_rna_matrix_status"))
              ),

              hr(),

              div(class = "upload-section",
                h4(icon("database"), "Preprocessed RNA Object"),
                fileInput(
                  "single_cell_rna_rds_file",
                  "RNA Seurat object (.rds)",
                  accept = c(".rds", ".RDS"),
                  buttonLabel = "Browse...",
                  placeholder = "No file selected"
                ),
                helpText("Upload a preprocessed single-cell RNA Seurat object saved as .rds."),
                div(class = "status-text", textOutput("single_cell_rna_rds_status"))
              ),

              hr(),

              div(class = "upload-section",
                h4(icon("chart-area"), "Single-cell ATAC Matrix"),
                fileInput(
                  "single_cell_atac_matrix_files",
                  "ATAC matrix files (10x format)",
                  accept = c(".mtx", ".mtx.gz", ".tsv", ".tsv.gz", ".bed", ".bed.gz", ".h5", ".hdf5"),
                  multiple = TRUE,
                  buttonLabel = "Browse...",
                  placeholder = "No files selected"
                ),
                helpText("Upload matrix.mtx(.gz), peaks.tsv(.gz), barcodes.tsv(.gz), or a single 10x HDF5 file."),
                div(class = "status-text", textOutput("single_cell_upload_atac_matrix_status"))
              ),

              hr(),

              div(class = "upload-section",
                h4(icon("layer-group"), "Preprocessed ATAC Object"),
                fileInput(
                  "single_cell_atac_rds_file",
                  "ATAC Seurat object (.rds)",
                  accept = c(".rds", ".RDS"),
                  buttonLabel = "Browse...",
                  placeholder = "No file selected"
                ),
                helpText("Upload a preprocessed single-cell ATAC Seurat object saved as .rds."),
                div(class = "status-text", textOutput("single_cell_atac_rds_status"))
              ),

              hr(),

              div(class = "action-buttons",
                actionButton("load_supplements", tagList(icon("download"), "Load Supplements"), class = "btn-info btn-block")
              )
            )
          )
        ),
        
        # Additional Files Section
        fluidRow(
          box(
            title = tagList(icon("folder-open"), "Additional Files"),
            width = 12,
            solidHeader = TRUE,
            status = "success",
            collapsible = TRUE,
            collapsed = TRUE,
            
            fluidRow(
              column(width = 4,
                div(class = "upload-section",
                  h4(icon("paperclip"), "ATAC Fragments"),
                  fileInput(
                    "single_cell_atac_fragments_file",
                    "Fragments file + index",
                    accept = c(".tsv", ".tsv.gz", ".tbi", ".tbi.gz"),
                    multiple = TRUE,
                    buttonLabel = "Browse...",
                    placeholder = "No files selected"
                  ),
                  helpText("Upload fragments.tsv(.gz) and its .tbi index file."),
                  div(class = "status-text", textOutput("single_cell_atac_fragments_status"))
                )
              ),

              column(width = 4,
                div(class = "upload-section",
                  h4(icon("link"), "RNA Mapping"),
                  fileInput(
                    "lineage_rna_mapping_file",
                    "Barcode mapping file",
                    accept = c(".rds", ".csv", ".tsv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Map lineage barcodes to single-cell RNA identifiers."),
                  div(class = "status-text", textOutput("lineage_rna_mapping_status"))
                )
              ),

              column(width = 4,
                div(class = "upload-section",
                  h4(icon("link"), "ATAC Mapping"),
                  fileInput(
                    "single_cell_atac_mapping_file",
                    "Barcode mapping file",
                    accept = c(".rds", ".csv", ".tsv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Optional mapping for ATAC barcodes."),
                  div(class = "status-text", textOutput("single_cell_atac_mapping_status"))
                )
              )
            ),

            fluidRow(
              column(width = 4,
                div(class = "upload-section",
                  h4(icon("table"), "ATAC Metadata"),
                  fileInput(
                    "single_cell_atac_metadata_file",
                    "Metadata CSV",
                    accept = c(".csv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Optional metadata with peak_region_fragments and other metrics."),
                  div(class = "status-text", textOutput("single_cell_atac_metadata_status"))
                )
              )
            ),

            hr(),
            
            fluidRow(
              column(width = 6,
                div(class = "upload-section",
                  h4(icon("capsules"), "External Drug Matrix"),
                  fileInput(
                    "drug_matrix_file",
                    "Drug response matrix",
                    accept = c(".rds", ".RDS", ".csv", ".tsv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Optional external drug-response matrix to integrate with existing data.")
                )
              )
            )
          )
        ),
        
        # Upload Summary Section
        fluidRow(
          box(
            title = tagList(icon("info-circle"), "Upload Summary"),
            width = 12,
            solidHeader = TRUE,
            status = "warning",
            
            div(class = "upload-summary-box",
              verbatimTextOutput("upload_summary")
            ),
            
            hr(),
            
            p(
              icon("exclamation-triangle"), 
              strong("Next Steps:"),
              br(),
              "1. After loading lineage data, proceed to", tags$b("QC Settings"), "to configure quality control parameters.",
              br(),
              "2. Apply denoising if needed to reduce noise in drug response data.",
              br(),
              "3. Load supplements (RNA/ATAC matrices and mapping files) to enable single-cell analysis.",
              br(),
              "4. Apply QC settings with RNA mapping to generate single-cell drug response data."
            )
          )
        )
      ),

      # =================================================================
      # 3.1 QC Settings Tab
      # =================================================================
      tabItem(tabName = "qc_settings",
        h2("QC Parameter Settings"),
        fluidRow(
          box(
            title = "Quality Control Parameters",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            fluidRow(
              column(
                width = 6,
                h4("scRNA-seq QC thresholds"),
                icon = icon("dna"),
                numericInput("rna_nfeature_min", "nFeature_RNA > ", value = 200, min = 0),
                helpText("Cells below this threshold will be flagged for removal."),
                numericInput("rna_nfeature_max", "nFeature_RNA < ", value = 5000, min = 1),
                helpText("Use a generous upper bound to keep heterogeneous cell states."),
                numericInput("rna_ncount_min", "nCount_RNA > ", value = 500, min = 0),
                numericInput("rna_ncount_max", "nCount_RNA < ", value = 30000, min = 1),
                numericInput("rna_percent_mt_max", "percent.mt < ", value = 20, min = 0, max = 100),
                helpText("High mitochondrial content indicates stressed or dying cells."),
                class = "control-spacer"
              ),
              column(
                width = 6,
                h4("scATAC-seq QC thresholds"),
                icon = icon("chart-area"),
                numericInput("atac_peak_fragments_min", "nCount_peaks > ", value = 2000, min = 0),
                numericInput("atac_peak_fragments_max", "nCount_peaks < ", value = 30000, min = 1),
                numericInput("atac_pct_reads_peaks_min", "pct_reads_in_peaks > ", value = 20, min = 0, max = 100),
                numericInput("atac_blacklist_ratio_max", "blacklist_ratio < ", value = 0.02, min = 0, max = 1, step = 0.01),
                numericInput("atac_nucleosome_signal_max", "nucleosome_signal < ", value = 4, min = 0),
                numericInput("atac_tss_enrichment_min", "TSS.enrichment > ", value = 2, min = 0),
                helpText("Adjust thresholds to balance cell retention with signal quality."),
                class = "control-spacer"
              )
            )
          )
        ),

        fluidRow(
          box(
            title = "Apply or reset settings",
            width = 12,
            solidHeader = TRUE,
            status = "success",
            div(
              class = "action-row",
              actionButton("apply_qc_settings", "Apply QC settings", class = "btn-primary"),
              actionButton("reset_qc_settings", "Reset to defaults", class = "btn-warning")
            ),
            div(class = "control-spacer", textOutput("qc_settings_status"))
          )
        )
      ),

      # =================================================================
      # 3.2 QC Tab
      # =================================================================
      tabItem(tabName = "qc",
        h2("Quality Control (QC)"),
        fluidRow(
          tabBox(
            id = "qc_tabset",
            width = 12,
            tabPanel("scRNA-seq QC",
              h3("scRNA-seq Quality Control Plots"),
              fluidRow(
                box(title = "QC metrics as a violin plot", width = 6, plotOutput("qc_rna_qc")),
                box(title = "PCA Elbow Plot", width = 6, plotOutput("qc_rna_pca_elbow"))
              ),
              fluidRow(
                box(title = "FeatureScatter", width = 12, plotOutput("qc_rna_fea"))
              )
            ),
            tabPanel("scATAC-seq QC",
              h3("scATAC-seq Quality Control Plots"),
              fluidRow(
                box(title = "DensityScatter", width = 6, plotOutput("qc_atac_Density")),
                box(title = "FragmentHistogram", width = 6, plotOutput("qc_atac_Fragment"))
              ),
              fluidRow(
                box(title = "QC metrics as a violin plot", width = 6, plotOutput("qc_atac_QC")),
                box(title = "Sequencing Depth Correlation", width = 6, plotOutput("qc_atac_depth_correlation"))
              )
            ),
            tabPanel("Data Integration",
              h3("scRNA and scATAC Integration"),
              plotOutput("integration_plot")
            )
          )
        )
      ),

      # =================================================================
      # 3.2 Lineage Level Tab
      # =================================================================
      tabItem(tabName = "lineage_level",

        fluidRow(
          # -- Left Control Panel --
          column(width = 3,
            box(
              title = "Visualization Control",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",

              # Clustering Method Selection
              selectInput("lineage_clustering_method", "Select Clustering Method:",
                choices = c("UMAP" = "umap", "PCA" = "pca", "t-SNE" = "tsne", "K-Means" = "kmeans")),

              conditionalPanel(
                condition = "input.lineage_clustering_method == 'kmeans'",
                textInput("lineage_kmeans_input", "K:", "e.g., 3")
              ),

              conditionalPanel(
                condition = "input.lineage_clustering_method == 'umap'",
                numericRangeInput("lineage_umap_pca_dims", "PCA dimensions", value = c(1, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.lineage_clustering_method == 'umap'",
                numericRangeInput("lineage_umap_svd_dims", "SVD dimensions", value = c(2, 30), min = 1, max = 50)
              ),

              conditionalPanel(
                condition = "input.lineage_clustering_method == 'tsne'",
                numericRangeInput("lineage_tsne_pca_dims", "PCA dimensions", value = c(1, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.lineage_clustering_method == 'tsne'",
                numericRangeInput("lineage_tsne_svd_dims", "SVD dimensions", value = c(2, 30), min = 1, max = 50)
              ),

              hr(),

              # UMAP Coloring Selector
              h4("UMAP Plot Coloring"),
              selectInput("lineage_color_by", "Select Coloring Variable Type:",
                choices = c(
                  "Cell Clusters" = "cluster",
                  "Drug Response" = "drug_response"
                )),

              # -- Conditional UI: Display different options based on coloring selection --
              # Drug Response/Growth Rate Selection
              conditionalPanel(
                condition = "input.lineage_color_by == 'drug_response'",
                selectInput("lineage_rds_object_select", "Select dataset",
                  choices = c("Long-term growth rate" = "cGR", "Drug resistance score" = "AUC")),
                pickerInput(
                  "lineage_drug_select",
                  "Select drugs",
                  choices = NULL,
                  multiple = TRUE,
                  options = pickerOptions(
                    actionsBox = TRUE,
                    liveSearch = TRUE,
                    noneSelectedText = "Awaiting uploaded metadata"
                  )
                )
              ),

              radioButtons(
                "lineage_bubble_mode",
                "Bubble plot mode",
                choices = c("Cluster" = "cluster", "Archetype" = "archetype"),
                inline = TRUE
              ),

              hr(),
              h4("Coloring/combining genomic tracks"),
              radioButtons(
                "lineage_track_mode",
                "Choose annotation mode",
                choices = c("Gene" = "gene", "Chromatin region" = "region"),
                inline = TRUE
              ),
              conditionalPanel(
                condition = "input.lineage_track_mode == 'gene'",
                textInput("lineage_gene_input", "Enter gene symbol", "e.g., SOX2")
              ),
              conditionalPanel(
                condition = "input.lineage_track_mode == 'region'",
                textInput("lineage_atac_region_input", "Enter chromatin region", "e.g., chr1:12345-23456")
              )

            )
          ),

          # -- Right Chart Display --
          column(width = 9,
            # -- Clustering Plots --
            fluidRow(
              box(title = "RNA-seq Clustering", width = 6, plotOutput("lineage_rna_cluster_plot")),
              box(title = "KNN Embedding", width = 6, plotOutput("lineage_knn_plot"))
            ),
            fluidRow(
              box(title = "ATAC-seq Clustering", width = 6, plotOutput("lineage_atac_cluster_plot")),
              box(title = "Drug Response Embedding", width = 6, plotOutput("lineage_drug_response_plot"))
            ),
            fluidRow(
              box(title = "Cluster/Archetype Bubble", width = 12, plotOutput("lineage_bubble_plot"))
            ),

            # -- Expression/Accessibility Analysis --
            fluidRow(
              box(title = "Combining genomic tracks", width = 12,
                tabBox(width = 12,
                  tabPanel("Cross-cluster violin", plotOutput("lineage_violin_plot")),
                  tabPanel("Coverage", plotOutput("lineage_coverage_plot"))
                )
              )
            )
          )
        )
      ),

      # =================================================================
      # 3.3 Single Cell Level Tab
      # =================================================================
      tabItem(tabName = "single_cell_level",
        fluidRow(
          # -- Left Control Panel --
          column(width = 3,
            box(
              title = "Visualization Control",
              width = NULL,
              solidHeader = TRUE,
              status = "primary",

              # Clustering Method Selection
              selectInput("single_clustering_method", "Select Clustering Method:",
                choices = c("UMAP" = "umap", "PCA" = "pca", "t-SNE" = "tsne", "K-Means" = "kmeans")),

              conditionalPanel(
                condition = "input.single_clustering_method == 'kmeans'",
                textInput("single_kmeans_input", "K:", "e.g., 3")
              ),

              conditionalPanel(
                condition = "input.single_clustering_method == 'umap'",
                numericRangeInput("single_umap_pca_dims", "PCA dimensions", value = c(1, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.single_clustering_method == 'umap'",
                numericRangeInput("single_umap_svd_dims", "SVD dimensions", value = c(2, 30), min = 1, max = 50)
              ),

              conditionalPanel(
                condition = "input.single_clustering_method == 'tsne'",
                numericRangeInput("single_tsne_pca_dims", "PCA dimensions", value = c(1, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.single_clustering_method == 'tsne'",
                numericRangeInput("single_tsne_svd_dims", "SVD dimensions", value = c(2, 30), min = 1, max = 50)
              ),

              hr(),

              # UMAP Coloring Selector
              h4("UMAP Plot Coloring"),
              selectInput("single_color_by", "Select Coloring Variable Type:",
                choices = c(
                  "Cell Clusters" = "cluster",
                  "Drug Response" = "drug_response"
                )),

              # -- Conditional UI: Display different options based on coloring selection --
              # Drug Response/Growth Rate Selection
              conditionalPanel(
                condition = "input.single_color_by == 'drug_response'",
                selectInput("single_rds_object_select", "Select dataset",
                  choices = c("Long-term growth rate" = "cGR", "Drug resistance score" = "AUC")),
                pickerInput(
                  "single_drug_select",
                  "Select drugs",
                  choices = NULL,
                  multiple = TRUE,
                  options = pickerOptions(
                    actionsBox = TRUE,
                    liveSearch = TRUE,
                    noneSelectedText = "Drug data requires RNA mapping - apply QC settings first"
                  )
                )
              ),

              hr(),
              h4("Coloring/combining genomic tracks"),
              radioButtons(
                "single_track_mode",
                "Choose annotation mode",
                choices = c("Gene" = "gene", "Chromatin region" = "region"),
                inline = TRUE
              ),
              conditionalPanel(
                condition = "input.single_track_mode == 'gene'",
                textInput("single_gene_input", "Enter gene symbol", "e.g., SOX2")
              ),
              conditionalPanel(
                condition = "input.single_track_mode == 'region'",
                textInput("single_atac_region_input", "Enter chromatin region", "e.g., chr1:12345-23456")
              )

            )
          ),

          # -- Right Chart Display --
          column(width = 9,
            # -- Clustering Plots --
            fluidRow(
              box(title = "RNA-seq Clustering", width = 6, plotOutput("single_rna_cluster_plot")),
              box(title = "KNN Embedding", width = 6, plotOutput("single_knn_plot"))
            ),
            fluidRow(
              box(title = "ATAC-seq Clustering", width = 6, plotOutput("single_atac_cluster_plot")),
              box(title = "Drug Response Embedding", width = 6, plotOutput("single_drug_response_plot"))
            ),

            # -- Expression/Accessibility Analysis --
            fluidRow(
              box(title = "Combining genomic tracks", width = 12,
                tabBox(width = 12,
                  tabPanel("Cross-cluster violin", plotOutput("single_violin_plot")),
                  tabPanel("Coverage", plotOutput("single_coverage_plot"))
                )
              )
            )
          )
        )
      ),

      # =================================================================
      # 3.4 Export Tab
      # =================================================================
      tabItem(tabName = "export",
        fluidRow(
          box(
            title = "Data export",
            width = 6,
            solidHeader = TRUE,
            status = "primary",

            p("Check availability before attempting to export."),
            infoBoxOutput("lineage_export_status", width = 12),
            conditionalPanel(
              condition = "output.lineage_data_available",
              downloadButton("download_lineage_data", "Export lineage data (.rds)", class = "btn-info")
            ),
            br(),
            infoBoxOutput("single_export_status", width = 12),
            conditionalPanel(
              condition = "output.single_cell_data_available",
              downloadButton("download_single_cell_data", "Export single-cell data (.rds)", class = "btn-info")
            ),
            hr(),
            h4("Analysis reports"),
            downloadButton("download_qc_summary", "Export QC report (.csv)", class = "btn-success"),
            br(),
            downloadButton("download_cluster_results", "Export clustering results (.csv)", class = "btn-success")
          ),

          # Plot Export
          box(
            title = "Plot export",
            width = 6,
            solidHeader = TRUE,
            status = "info",

            h4("QC plots"),
            selectInput("export_qc_format", "Image format", choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
            downloadButton("download_qc_plots", "Export QC plots", class = "btn-warning"),

            hr(),

            h4("Clustering plots"),
            selectInput("export_lineage_format", "Lineage format", choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
            downloadButton("download_lineage_plots", "Export lineage plots", class = "btn-warning"),
            br(),
            selectInput("export_single_format", "Single-cell format", choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
            downloadButton("download_single_plots", "Export single-cell plots", class = "btn-warning"),

            hr(),

            h4("Custom export"),
            textInput("custom_plot_name", "Plot name", "custom_plot"),
            selectInput("custom_plot_type", "Available plots", choices = c(
              "RNA clustering" = "rna_cluster",
              "Drug embedding" = "drug_embedding",
              "KNN embedding" = "knn_embedding",
              "Bubble" = "bubble"
            )),
            selectInput("custom_plot_level", "Analysis level", choices = c("Lineage" = "lineage", "Single-cell" = "single")),
            selectInput("custom_plot_format", "Image format", choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
            downloadButton("download_custom_plot", "Export custom plot", class = "btn-danger")
          )
        ),

        fluidRow(
          box(
            title = "Export Status",
            width = 12,
            solidHeader = TRUE,
            status = "success",

            h4("Export readiness"),
            fluidRow(
              column(6, infoBoxOutput("qc_settings_info", width = 12)),
              column(6, infoBoxOutput("custom_export_info", width = 12))
            ),

            hr(),

            h4("Export history"),
            div(class = "export-history-table", tableOutput("export_history_table")),

            hr(),

            actionButton("clear_export_history", "Clear history", class = "btn-secondary")
          )
        )
      )
    )
  )
)
