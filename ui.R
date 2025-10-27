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
        menuSubItem("Lineage Level", tabName = "lineage_level", icon = icon("sitemap"))
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
          # Lineage data upload and processing
          column(width = 12,
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
              column(
                width = 6,
                div(class = "upload-section",
                  h4(icon("dna"), "Single-cell RNA Inputs"),
                  fileInput(
                    "single_cell_rna_matrix_files",
                    "RNA matrix files (10x format)",
                    accept = c(".mtx", ".mtx.gz", ".tsv", ".tsv.gz", ".h5", ".hdf5"),
                    multiple = TRUE,
                    buttonLabel = "Browse...",
                    placeholder = "No files selected"
                  ),
                  helpText("Upload matrix.mtx(.gz), features.tsv(.gz), barcodes.tsv(.gz), or a single 10x HDF5 file."),
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
                  helpText("Optional: upload a preprocessed single-cell RNA Seurat object."),
                  div(class = "status-text", textOutput("single_cell_rna_rds_status"))
                ),

                hr(),

                div(class = "upload-section",
                  h4(icon("project-diagram"), "RNA Barcode Mapping"),
                  fileInput(
                    "lineage_rna_mapping_file",
                    "Barcode mapping file",
                    accept = c(".rds", ".RDS", ".csv", ".tsv", ".txt"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Map single-cell barcodes to lineage barcodes to enable pseudo-bulk aggregation."),
                  div(class = "status-text", textOutput("lineage_rna_mapping_status"))
                )
              ),

              column(
                width = 6,
                div(class = "upload-section",
                  h4(icon("chart-area"), "Single-cell ATAC Inputs"),
                  fileInput(
                    "single_cell_atac_matrix_files",
                    "ATAC matrix files (10x format)",
                    accept = c(".mtx", ".mtx.gz", ".tsv", ".tsv.gz", ".h5", ".hdf5"),
                    multiple = TRUE,
                    buttonLabel = "Browse...",
                    placeholder = "No files selected"
                  ),
                  helpText("Upload matrix.mtx(.gz), peaks.tsv(.gz), barcodes.tsv(.gz), or a single 10x HDF5 file."),
                  div(class = "status-text", textOutput("single_cell_upload_atac_matrix_status"))
                ),

                hr(),

                div(class = "upload-section",
                  h4(icon("stream"), "ATAC Fragments & Metadata"),
                  fileInput(
                    "single_cell_atac_fragments_file",
                    "Fragments file (.tsv/.tsv.gz)",
                    accept = c(".tsv", ".tsv.gz"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  fileInput(
                    "single_cell_atac_metadata_file",
                    "Metadata (.csv)",
                    accept = c(".csv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Provide the fragments.tsv(.gz) file and optional metadata for ATAC processing."),
                  div(class = "status-text", textOutput("single_cell_atac_fragments_status")),
                  div(class = "status-text", textOutput("single_cell_atac_metadata_status"))
                ),

                hr(),

                div(class = "upload-section",
                  h4(icon("database"), "Preprocessed ATAC Object"),
                  fileInput(
                    "single_cell_atac_rds_file",
                    "ATAC Seurat object (.rds)",
                    accept = c(".rds", ".RDS"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Optional: upload a preprocessed single-cell ATAC Seurat object."),
                  div(class = "status-text", textOutput("single_cell_atac_rds_status"))
                ),

                hr(),

                div(class = "upload-section",
                  h4(icon("project-diagram"), "ATAC Barcode Mapping"),
                  fileInput(
                    "single_cell_atac_mapping_file",
                    "Barcode mapping file",
                    accept = c(".rds", ".RDS", ".csv", ".tsv", ".txt"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Map ATAC single-cell barcodes to lineage barcodes."),
                  div(class = "status-text", textOutput("single_cell_atac_mapping_status"))
                )
              )
            ),

            hr(),

            fluidRow(
              column(
                width = 6,
                div(class = "upload-section",
                  h4(icon("capsules"), "External Drug Matrix"),
                  fileInput(
                    "drug_matrix_file",
                    "Drug response matrix",
                    accept = c(".rds", ".RDS", ".csv", ".tsv"),
                    buttonLabel = "Browse...",
                    placeholder = "No file selected"
                  ),
                  helpText("Optional external lineage-level drug-response matrix to merge with the dataset."),
                  div(class = "status-text", textOutput("drug_matrix_status"))
                )
              ),

              column(
                width = 6,
                div(class = "upload-section",
                  h4(icon("tasks"), "Actions"),
                  helpText("Load supplements to process single-cell inputs into pseudo-bulk lineage data."),
                  div(
                    class = "action-buttons",
                    actionButton("load_supplements", tagList(icon("download"), "Load Supplements"), class = "btn-info"),
                    actionButton("clear_single_cell", tagList(icon("trash"), "Clear Single-cell"), class = "btn-danger")
                  ),
                  div(class = "status-text", textOutput("single_cell_upload_status"))
                )
              )
            )
          )
        ),

        fluidRow(
          column(
            width = 12,
            box(
              title = tagList(icon("save"), "QC Snapshot"),
              width = NULL,
              solidHeader = TRUE,
              status = "success",
              div(
                class = "upload-section",
                h4(icon("history"), "Restore Previous Session"),
                fileInput(
                  "qc_snapshot_file",
                  "QC snapshot (.rds)",
                  accept = c(".rds", ".RDS"),
                  buttonLabel = "Browse...",
                  placeholder = "No file selected"
                ),
                helpText("Load a previously saved QC snapshot to restore processed data without re-running uploads or QC."),
                div(class = "status-text", textOutput("qc_snapshot_status")),
                div(
                  class = "action-buttons",
                  actionButton("load_qc_snapshot", tagList(icon("redo"), "Load Snapshot"), class = "btn-success")
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
              "1. After loading lineage data, proceed to", tags$b("QC Settings"), "to configure quality control parameters if required.",
              br(),
              "2. Apply denoising if needed to reduce noise in drug response data.",
              br(),
              "3. Explore lineage-level visualizations in the Main tab."
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
                box(title = "Sequencing Depth Correlation", width = 6, plotOutput("qc_atac_depth_correlation"))
              ),
              fluidRow(
                box(title = "QC metrics as a violin plot", width = 12, plotOutput("qc_atac_QC"))
              )
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

              h4("Reduction and clustering"),
              selectInput("lineage_red_method", "Select reduction method",
                choices = c("UMAP" = "umap", "PCA/SVD" = "pca", "t-SNE" = "tsne"),
                selected = "umap"
              ),
              selectInput("lineage_clustering_method", "Select clustering method",
                choices = c("Louvain" = "louvain", "K-means" = "kmeans"),
                selected = "louvain"
              ),
              conditionalPanel(
                condition = "input.lineage_clustering_method == 'kmeans'",
                textInput("lineage_kmeans_input", "K for K-means", value = 5)
              ),
              conditionalPanel(
                condition = "input.lineage_red_method == 'umap'",
                numericRangeInput("lineage_umap_pca_dims", "PCA dimensions", value = c(1, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.lineage_red_method == 'umap'",
                numericRangeInput("lineage_umap_svd_dims", "SVD dimensions", value = c(2, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.lineage_red_method == 'tsne'",
                numericRangeInput("lineage_tsne_pca_dims", "PCA dimensions", value = c(1, 30), min = 1, max = 50)
              ),
              conditionalPanel(
                condition = "input.lineage_red_method == 'tsne'",
                numericRangeInput("lineage_tsne_svd_dims", "SVD dimensions", value = c(2, 30), min = 1, max = 50)
              ),

              hr(),

              h4("Dataset selection"),
              pickerInput(
                "lineage_rds_object_select",
                "Select dataset",
                choices = NULL,
                multiple = FALSE,
                options = pickerOptions(
                  liveSearch = TRUE,
                  noneSelectedText = "Awaiting lineage assays"
                )
              ),

              h4("Drug selection"),
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
              ),

              radioButtons(
                "lineage_bubble_mode",
                "Bubble plot mode",
                choices = c("Cluster" = "cluster", "Archetype" = "archetype"),
                inline = TRUE
              ),

              hr(),

              h4("Genomic tracks"),
              radioButtons(
                "lineage_track_mode",
                "Annotation mode",
                choices = c("Gene" = "gene", "Chromatin region" = "region"),
                inline = TRUE
              ),
              conditionalPanel(
                condition = "input.lineage_track_mode == 'gene'",
                textInput("lineage_gene_input", "Enter gene symbol", value = "MYC")
              ),
              conditionalPanel(
                condition = "input.lineage_track_mode == 'region'",
                textInput("lineage_atac_region_input", "Enter chromatin region", value = "chr1 183993-184842")
              ),

              div(class = "control-spacer",
                actionButton(
                  "lineage_refresh_plots",
                  tagList(icon("sync"), "Refresh"),
                  class = "btn-primary btn-block"
                )
              )

            )
          ),

          # -- Right Chart Display --
          column(width = 9,
            fluidRow(
              box(title = "RNA-seq clustering", width = 6, plotOutput("lineage_rna_cluster_plot")),
              box(title = "KNN embedding", width = 6, plotOutput("lineage_knn_plot"))
            ),
            fluidRow(
              box(title = "ATAC-seq clustering", width = 6, plotOutput("lineage_atac_cluster_plot")),
              box(title = "Drug response embedding", width = 6, plotOutput("lineage_drug_response_plot"))
            ),
            fluidRow(
              box(title = "Cluster/Archetype bubble", width = 12, plotOutput("lineage_bubble_plot"))
            ),
            fluidRow(
              box(title = "Genomic track comparison", width = 12,
                tabBox(width = 12,
                  tabPanel("Cross-cluster violin", plotOutput("lineage_violin_plot"))
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

            hr(),

            h4("Custom export"),
            textInput("custom_plot_name", "Plot name", "custom_plot"),
            selectInput("custom_plot_type", "Available plots", choices = c(
              "RNA clustering" = "rna_cluster",
              "ATAC clustering" = "atac_cluster",
              "Drug embedding" = "drug_embedding",
              "KNN embedding" = "knn_embedding",
              "Bubble" = "bubble",
              "Violin" = "violin"
            )),
            selectInput("custom_plot_level", "Analysis level", choices = c("Lineage" = "lineage")),
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
