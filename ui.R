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
        "        .control-spacer {
          margin-top: 15px;
        }

        .action-row .btn {
          margin-right: 10px;
        }

        .export-history-table {
          margin-top: 10px;
        }
      "
      ))
    ),
    tabItems(
      # =================================================================
      # 3.0 Data Upload Tab
      # =================================================================
      tabItem(tabName = "upload",
        h2("Data Upload"),
        fluidRow(
          box(
            title = "Upload Experiment Outputs",
            width = 12,
            solidHeader = TRUE,
            status = "primary",

            tabBox(
              width = 12,
              id = "upload_groups",
              tabPanel(
                title = "Required Files",
                fluidRow(
                  column(
                    width = 6,
                    h4("Lineage Level"),
                    fileInput(
                      "lineage_rds_file",
                      "Select lineage drug-response RDS file",
                      accept = c(".rds", ".RDS"),
                      buttonLabel = "Browse"
                    ),
                    helpText("Upload the processed lineage-level Seurat/SingleCellExperiment object."),
                    textOutput("lineage_upload_status")
                  ),
                  column(
                    width = 6,
                    h4("Single-cell Level"),
                    helpText("Single-cell response objects are derived from uploaded matrices."),
                    textOutput("single_cell_upload_status")
                  )
                )
              ),
              tabPanel(
                title = "Optional Supplements",
                fluidRow(
                  column(
                    width = 6,
                    fileInput(
                      "single_cell_rna_matrix_files",
                      "Select single-cell RNA matrix files",
                      accept = c(".mtx", ".mtx.gz", ".tsv", ".tsv.gz"),
                      multiple = TRUE,
                      buttonLabel = "Browse"
                    ),
                    helpText("Upload the matrix.mtx(.gz), features.tsv(.gz), and barcodes.tsv(.gz) files."),
                    textOutput("single_cell_upload_rna_matrix_status")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    fileInput(
                      "single_cell_atac_matrix_files",
                      "Select single-cell ATAC matrix files",
                      accept = c(".mtx", ".mtx.gz", ".tsv", ".tsv.gz", ".bed", ".bed.gz"),
                      multiple = TRUE,
                      buttonLabel = "Browse"
                    ),
                    helpText("Upload the matrix.mtx(.gz), peaks/features.tsv(.gz), and barcodes.tsv(.gz) files."),
                    textOutput("single_cell_upload_atac_matrix_status")
                  ),
                  column(
                    width = 6,
                    fileInput(
                      "single_cell_atac_metadata_file",
                      "Select single-cell ATAC metadata (.csv/.tsv)",
                      accept = c(".csv", ".tsv", ".txt"),
                      buttonLabel = "Browse"
                    ),
                    helpText("Provide the metadata file matching the ATAC matrix barcodes."),
                    textOutput("single_cell_atac_metadata_status")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    fileInput(
                      "single_cell_atac_fragments_file",
                      "Select single-cell ATAC fragments file",
                      accept = c(".tsv", ".tsv.gz"),
                      buttonLabel = "Browse"
                    ),
                    helpText("Upload the fragments.tsv(.gz) file for ATAC integration."),
                    textOutput("single_cell_atac_fragments_status")
                  ),
                  column(
                    width = 6,
                    fileInput(
                      "lineage_rna_mapping_file",
                      "Select lineage/single-cell RNA barcode mapping file",
                      accept = c(".rds", ".csv", ".tsv"),
                      buttonLabel = "Browse"
                    ),
                    helpText("Map lineage-level barcodes to single-cell identifiers."),
                    textOutput("lineage_rna_mapping_status")
                  ),
                  column(
                    width = 6,
                    fileInput(
                      "single_cell_atac_mapping_file",
                      "Select lineage/single-cell ATAC barcode mapping file",
                      accept = c(".rds", ".csv", ".tsv"),
                      buttonLabel = "Browse"
                    ),
                    helpText("Optional mapping to reconcile ATAC barcodes across assays."),
                    textOutput("single_cell_atac_mapping_status")
                  )
                ),
                fluidRow(
                  column(
                    width = 6,
                    fileInput(
                      "drug_matrix_file",
                      "Optional drug-response matrix",
                      accept = c(".rds", ".RDS", ".csv", ".tsv"),
                      buttonLabel = "Browse"
                    ),
                    checkboxInput("denoise_drug_matrix", "Apply denoising to uploaded drug matrix", FALSE)
                  ),
                  column(
                    width = 6,
                    pickerInput(
                      "dslt_denoise_assays",
                      "Select assays to denoise",
                      choices = NULL,
                      multiple = TRUE,
                      options = pickerOptions(actionsBox = TRUE, noneSelectedText = "Awaiting lineage upload")
                    )
                  )
                )
              )
            ),

            hr(),

            fluidRow(
              column(
                width = 4,
                textInput("barcode_suffix", "Barcode suffix filter", value = "-1$",
                  placeholder = "Regex suffix, e.g. -1$")
              ),
              column(
                width = 4,
                pickerInput(
                  "denoise_assays_choice",
                  "Denoise lineage assays",
                  choices = NULL,
                  multiple = TRUE,
                  options = pickerOptions(actionsBox = TRUE, noneSelectedText = "Awaiting lineage upload")
                )
              ),
              column(
                width = 4,
                actionButton("load_all_data", "Load selected data", class = "btn-success btn-block")
              )
            ),

            fluidRow(
              column(
                width = 12,
                h4("Upload Summary"),
                verbatimTextOutput("upload_summary")
              )
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
                numericInput("rna_nfeature_min", "Minimum detected genes", value = 200, min = 0),
                helpText("Cells below this threshold will be flagged for removal."),
                numericInput("rna_nfeature_max", "Maximum detected genes", value = 5000, min = 1),
                helpText("Use a generous upper bound to keep heterogeneous cell states."),
                numericInput("rna_ncount_min", "Minimum UMI counts", value = 500, min = 0),
                numericInput("rna_ncount_max", "Maximum UMI counts", value = 30000, min = 1),
                numericInput("rna_percent_mt_max", "Maximum mitochondrial percentage (%)", value = 20, min = 0, max = 100),
                helpText("High mitochondrial content indicates stressed or dying cells."),
                class = "control-spacer"
              ),
              column(
                width = 6,
                h4("scATAC-seq QC thresholds"),
                icon = icon("chart-area"),
                numericInput("atac_peak_fragments_min", "Minimum peak fragments", value = 1000, min = 0),
                numericInput("atac_peak_fragments_max", "Maximum peak fragments", value = 100000, min = 1),
                numericInput("atac_pct_reads_peaks_min", "Minimum reads in peaks (%)", value = 40, min = 0, max = 100),
                numericInput("atac_blacklist_ratio_max", "Maximum blacklist ratio (%)", value = 5, min = 0, max = 100),
                numericInput("atac_nucleosome_signal_max", "Maximum nucleosome signal", value = 2, min = 0),
                numericInput("atac_tss_enrichment_min", "Minimum TSS enrichment", value = 2, min = 0),
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
                    noneSelectedText = "Awaiting uploaded metadata"
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
