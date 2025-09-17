library(shiny)
library(shinydashboard)
library(shinyWidgets)

# UI Definition
ui <- dashboardPage(
  skin = "blue",
  
  # 1. Dashboard Header
  dashboardHeader(title = "Name?"),
  
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
    tags$style(HTML("
      .box .box-body {
        text-align: center;
      }
    "))
  ),
    tabItems(
      # =================================================================
      # 3.0 Data Upload Tab
      # =================================================================
      tabItem(tabName = "upload",
        h2("Data Upload"),
        fluidRow(
          box(
            title = "RDS File Upload",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            
            fluidRow(
              column(width = 6,
                h4("Lineage Level Data"),
                fileInput("lineage_rds_file", 
                         "Select Lineage RDS File:",
                         accept = c(".rds", ".RDS")),
                textOutput("lineage_upload_status")
              ),
              column(width = 6,
                h4("Single Cell Level Data"),
                fileInput("single_cell_rds_file", 
                         "Select Single Cell RDS File:",
                         accept = c(".rds", ".RDS")),
                textOutput("single_cell_upload_status")
              )
            ),
            
            hr(),
            
            fluidRow(
              column(width = 12,
                h4("Upload Status"),
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
          # RNA-seq QC Parameters
          box(
            title = "RNA-seq QC Parameters",
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Gene Expression Related"),
            numericInput("rna_nfeature_min", "nFeature_RNA Minimum:", value = 200, min = 0),
            numericInput("rna_nfeature_max", "nFeature_RNA Maximum:", value = 5000, min = 1),
            
            h4("UMI Count Related"),
            numericInput("rna_ncount_min", "nCount_RNA Minimum:", value = 500, min = 0),
            numericInput("rna_ncount_max", "nCount_RNA Maximum:", value = 30000, min = 1),
            
            h4("Mitochondrial Gene Percentage"),
            numericInput("rna_percent_mt_max", "percent.mt Maximum (%):", value = 20, min = 0, max = 100)
          ),
          
          # ATAC-seq QC Parameters
          box(
            title = "ATAC-seq QC Parameters",
            width = 6,
            solidHeader = TRUE,
            status = "info",
            
            h4("Peak Region Fragments"),
            numericInput("atac_peak_fragments_min", "peak_region_fragments Minimum:", value = 1000, min = 0),
            numericInput("atac_peak_fragments_max", "peak_region_fragments Maximum:", value = 100000, min = 1),
            
            h4("Percentage Reads in Peaks"),
            numericInput("atac_pct_reads_peaks_min", "pct_reads_in_peaks Minimum (%):", value = 40, min = 0, max = 100),
            
            h4("Blacklist Region Ratio"),
            numericInput("atac_blacklist_ratio_max", "blacklist_ratio Maximum (%):", value = 5, min = 0, max = 100),
            
            h4("Nucleosome Signal"),
            numericInput("atac_nucleosome_signal_max", "nucleosome_signal Maximum:", value = 2, min = 0),
            
            h4("TSS Enrichment"),
            numericInput("atac_tss_enrichment_min", "TSS.enrichment Minimum:", value = 2, min = 0)
          )
        ),
        
        fluidRow(
          box(
            title = "Operations",
            width = 12,
            solidHeader = TRUE,
            status = "success",
            
            actionButton("apply_qc_settings", "Apply QC Settings", class = "btn-primary"),
            br(), br(),
            actionButton("reset_qc_settings", "Reset to Default Values", class = "btn-warning"),
            br(), br(),
            textOutput("qc_settings_status")
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
            title = "QC Report",
            id = "qc_tabset",
            width = 12,
            tabPanel("scRNA-seq QC", 
                     h3("scRNA-seq Quality Control Plots"),
                     fluidRow(box(title = "QC metrics as a violin plot", width=12, plotOutput("qc_rna_qc"))),
                     fluidRow(box(title = "FeatureScatter", width=12, plotOutput("qc_rna_fea")))
            ),
            tabPanel("scATAC-seq QC", 
                     h3("scATAC-seq Quality Control Plots"),
                     fluidRow(
                       box(title = "DensityScatter", width=6, plotOutput("qc_atac_Density")),
                       box(title = "FragmentHistogram", width=6, plotOutput("qc_atac_Fragment"))
                     ),
                     fluidRow(box(title = "QC metrics as a violin plot", width=12, plotOutput("qc_atac_QC")))
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
                selectInput("lineage_rds_object_select", "Select Dataset (From R list object):",
                            choices = c("Long-term Growth Rate" = "cGR", "Drug Resistance Score" = "AUC")),
                selectInput("lineage_drug_select", "Select Drug:",
                            choices = c("Drug A", "Drug B", "Drug C")) # Mock data
              ),

              hr(),
              h4("Coloring/Combining genomic tracks"),
              textInput("lineage_gene_input", "Enter Gene Name:", "e.g., SOX2"),
              h4("OR"),
              textInput("lineage_atac_region_input", "Enter Chromatin Region:", "e.g., chr1:12345-23456")

            )
          ),
          
          # -- Right Chart Display --
          column(width = 9,
            # -- Clustering Plots --
            fluidRow(
              box(title = "RNA-seq Clustering", width = 6, plotOutput("lineage_rna_cluster_plot", height = "350px")),
              box(title = "Cluster Density Plot", width = 6, plotOutput("lineage_density_plot", height = "350px"))
            ),
            fluidRow(
              box(title = "ATAC-seq Clustering", width = 6, plotOutput("lineage_atac_cluster_plot", height = "350px")),
              box(title = "Drug Response Heatmap", width = 6, plotOutput("lineage_drug_response_heatmap", height = "350px"))
            ),
            
            # -- Expression/Accessibility Analysis --
            fluidRow(
              box(title = "Combining genomic tracks ", width = 12,
                  fluidRow(
                    box(title = "Cross-Cluster Expression/Accessibility Violin Plot", width = 12, solidHeader = TRUE, plotOutput("lineage_violin_plot", height = "350px"))
                  ),
                  fluidRow(
                    box(title = "Coverage Plot", width = 12, solidHeader = TRUE, plotOutput("lineage_coverage_plot", height = "300px"))
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
                selectInput("single_rds_object_select", "Select Dataset (From R list object):",
                            choices = c("Long-term Growth Rate" = "cGR", "Drug Resistance Score" = "AUC")),
                selectInput("single_drug_select", "Select Drug:",
                            choices = c("Drug A", "Drug B", "Drug C")) # Mock data
              ),

              hr(),
              h4("Coloring/Combining genomic tracks"),
              textInput("single_gene_input", "Enter Gene Name:", "e.g., SOX2"),
              h4("OR"),
              textInput("single_atac_region_input", "Enter Chromatin Region:", "e.g., chr1:12345-23456")

            )
          ),
          
          # -- Right Chart Display --
          column(width = 9,
            # -- Clustering Plots --
            fluidRow(
              box(title = "scRNA-seq Clustering", width = 6, plotOutput("single_rna_cluster_plot", height = "350px")),
              box(title = "Cluster Density Plot", width = 6, plotOutput("single_density_plot", height = "350px"))
            ),
            fluidRow(
              box(title = "scATAC-seq Clustering", width = 6, plotOutput("single_atac_cluster_plot", height = "350px")),
              box(title = "Drug Response Heatmap", width = 6, plotOutput("single_drug_response_heatmap", height = "350px"))
            ),
            
            # -- Expression/Accessibility Analysis --
            fluidRow(
              box(title = "Combining genomic tracks ", width = 12,
                  fluidRow(
                    box(title = "Coverage Plot", width = 12, solidHeader = TRUE, plotOutput("single_coverage_plot", height = "900px"))
                  )
              )
            )
          )
        )
      ),
      
      # =================================================================
      # 3.4 File Export Tab
      # =================================================================
      tabItem(tabName = "export",
        h2("File Export"),
        fluidRow(
          # Data Export
          box(
            title = "Data Export",
            width = 6,
            solidHeader = TRUE,
            status = "primary",
            
            h4("Processed Data"),
            p("Export QC-processed and analyzed data files"),
            
            fluidRow(
              column(width = 12,
                h5("Lineage Level Data"),
                conditionalPanel(
                  condition = "output.lineage_data_available",
                  downloadButton("download_lineage_data", "Export Lineage Data (.rds)", 
                                class = "btn-info", style = "margin-bottom: 10px;")
                ),
                conditionalPanel(
                  condition = "!output.lineage_data_available",
                  p("Lineage data not ready", style = "color: #999;")
                )
              )
            ),
            
            fluidRow(
              column(width = 12,
                h5("Single Cell Level Data"),
                conditionalPanel(
                  condition = "output.single_cell_data_available",
                  downloadButton("download_single_cell_data", "Export Single Cell Data (.rds)", 
                                class = "btn-info", style = "margin-bottom: 10px;")
                ),
                conditionalPanel(
                  condition = "!output.single_cell_data_available",
                  p("Single Cell data not ready", style = "color: #999;")
                )
              )
            ),
            
            hr(),
            
            h4("Analysis Results"),
            downloadButton("download_qc_summary", "Export QC Report (.csv)", 
                          class = "btn-success", style = "margin-bottom: 10px;"),
            br(),
            downloadButton("download_cluster_results", "Export Clustering Results (.csv)", 
                          class = "btn-success")
          ),
          
          # Plot Export
          box(
            title = "Plot Export",
            width = 6,
            solidHeader = TRUE,
            status = "info",
            
            h4("QC Plots"),
            selectInput("export_qc_format", "Select Image Format:",
                       choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
            downloadButton("download_qc_plots", "Export QC Plots", 
                          class = "btn-warning", style = "margin-bottom: 10px;"),
            
            hr(),
            
            h4("Clustering Analysis Plots"),
            fluidRow(
              column(width = 6,
                h5("Lineage Level"),
                selectInput("export_lineage_format", "Image Format:",
                           choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
                downloadButton("download_lineage_plots", "Export Lineage Plots", 
                              class = "btn-warning", style = "margin-bottom: 5px;")
              ),
              column(width = 6,
                h5("Single Cell Level"),
                selectInput("export_single_format", "Image Format:",
                           choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
                downloadButton("download_single_plots", "Export Single Cell Plots", 
                              class = "btn-warning", style = "margin-bottom: 5px;")
              )
            ),
            
            hr(),
            
            h4("Custom Plot Export"),
            textInput("custom_plot_name", "Plot Name:", "custom_plot"),
            selectInput("custom_plot_type", "Select Plot to Export:",
                       choices = c(
                         "RNA Clustering Plot" = "rna_cluster",
                         "ATAC Clustering Plot" = "atac_cluster", 
                         "Density Plot" = "density",
                         "Drug Response Heatmap" = "drug_heatmap",
                         "Violin Plot" = "violin",
                         "Coverage Plot" = "coverage"
                       )),
            selectInput("custom_plot_level", "Select Analysis Level:",
                       choices = c("Lineage Level" = "lineage", "Single Cell Level" = "single")),
            selectInput("custom_plot_format", "Image Format:",
                       choices = c("PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg")),
            downloadButton("download_custom_plot", "Export Custom Plot", 
                          class = "btn-danger")
          )
        ),
        
        fluidRow(
          box(
            title = "Export Status",
            width = 12,
            solidHeader = TRUE,
            status = "success",
            
            h4("Export History"),
            verbatimTextOutput("export_history"),
            
            hr(),
            
            actionButton("clear_export_history", "Clear History", class = "btn-secondary")
          )
        )
      )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # --- Data Storage ---
  values <- reactiveValues(
    lineage_data = NULL,
    single_cell_data = NULL,
    qc_applied = FALSE,
    export_history = c()
  )
  
  # --- File Upload Processing ---
  observeEvent(input$lineage_rds_file, {
    req(input$lineage_rds_file)
    tryCatch({
      values$lineage_data <- readRDS(input$lineage_rds_file$datapath)
      output$lineage_upload_status <- renderText("Lineage data uploaded successfully!")
    }, error = function(e) {
      output$lineage_upload_status <- renderText(paste("Upload failed:", e$message))
    })
  })
  
  observeEvent(input$single_cell_rds_file, {
    req(input$single_cell_rds_file)
    tryCatch({
      values$single_cell_data <- readRDS(input$single_cell_rds_file$datapath)
      output$single_cell_upload_status <- renderText("Single Cell data uploaded successfully!")
    }, error = function(e) {
      output$single_cell_upload_status <- renderText(paste("Upload failed:", e$message))
    })
  })
  
  # --- Upload Status Summary ---
  output$upload_summary <- renderText({
    lineage_status <- if(is.null(values$lineage_data)) "Not uploaded" else "Uploaded"
    single_cell_status <- if(is.null(values$single_cell_data)) "Not uploaded" else "Uploaded"
    
    paste("Lineage Data:", lineage_status, "\n",
          "Single Cell Data:", single_cell_status, "\n",
          "Status: Data", if(is.null(values$lineage_data) && is.null(values$single_cell_data)) "not ready" else "ready")
  })
  
  # --- QC Settings Application ---
  observeEvent(input$apply_qc_settings, {
    values$qc_applied <- TRUE
    output$qc_settings_status <- renderText("QC settings applied successfully!")
  })
  
  # --- QC Settings Reset ---
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
  
  # --- QC 图表 (修改为显示本地图片) ---
  output$qc_atac_Density <- renderImage({
    list(src = "images/qc_atac_density.png",
         contentType = 'image/png',
         width = 500,
         height = 400,
         alt = "ATAC DensityScatter")
  }, deleteFile = FALSE)
  
  output$qc_atac_Fragment <- renderImage({
    list(src = "images/qc_atac_fragment.png",
         contentType = 'image/png',
         width = 500,
         height = 400,
         alt = "Fragment Length Distribution")
  }, deleteFile = FALSE)

  output$qc_atac_QC <- renderImage({
    list(src = "images/qc_atac_QC.png",
         contentType = 'image/png',
         width = 900,
         height = 300,
         alt = "ATAC Quality Control")
  }, deleteFile = FALSE)
  
  # --- Lineage Level 图表 (修改为显示本地图片) ---
  output$lineage_rna_cluster_plot <- renderImage({
    list(src = "images/lineage_rna_cluster.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Lineage RNA UMAP")
  }, deleteFile = FALSE)
  
  output$lineage_atac_cluster_plot <- renderImage({
    list(src = "images/lineage_atac_cluster.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Lineage ATAC UMAP")
  }, deleteFile = FALSE)
  
  output$lineage_density_plot <- renderImage({
    list(src = "images/lineage_density.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Lineage 密度图")
  }, deleteFile = FALSE)
  
  output$lineage_drug_response_heatmap <- renderImage({
    list(src = "images/lineage_drug_heatmap.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Lineage 药物反应热图")
  }, deleteFile = FALSE)
  
  output$lineage_violin_plot <- renderImage({
    list(src = "images/lineage_violin.png",
         contentType = 'image/png',
         width = 600,
         height = 350,
         alt = "Lineage 小提琴图")
  }, deleteFile = FALSE)
  
  output$lineage_coverage_plot <- renderImage({
    list(src = "images/lineage_coverage.png",
         contentType = 'image/png',
         width = 600,
         height = 300,
         alt = "Lineage Coverage Plot")
  }, deleteFile = FALSE)
  
  # --- Single Cell Level 图表 (修改为显示本地图片) ---
  output$single_rna_cluster_plot <- renderImage({
    list(src = "images/single_rna_cluster.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Single Cell RNA UMAP")
  }, deleteFile = FALSE)
  
  output$single_atac_cluster_plot <- renderImage({
    list(src = "images/single_atac_cluster.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Single Cell ATAC UMAP")
  }, deleteFile = FALSE)
  
  output$single_density_plot <- renderImage({
    list(src = "images/single_density.png",
         contentType = 'image/png',
         width = 350,
         height = 350,
         alt = "Single Cell 密度图")
  }, deleteFile = FALSE)
  
  output$single_drug_response_heatmap <- renderImage({
    list(src = "images/single_drug_heatmap.png",
         contentType = 'image/png',
         width = 350,
         height = 250,
         alt = "Single Cell 药物反应热图")
  }, deleteFile = FALSE)

  output$single_coverage_plot <- renderImage({
    list(src = "images/single_coverage.png",
         contentType = 'image/png',
         width = 600,
         height = 800,
         alt = "Single Cell Coverage Plot")
  }, deleteFile = FALSE)
  
  # --- QC 图表 (修改为显示本地图片) ---
  output$qc_rna_qc <- renderImage({
    list(src = "images/qc_rna_qc.png",
         contentType = 'image/png',
         width = 800,
         height = 400,
         alt = "UMI 计数")
  }, deleteFile = FALSE)
  
  output$qc_rna_fea <- renderImage({
    list(src = "images/qc_rna_features.png",
         contentType = 'image/png',
         width = 800,
         height = 400,
         alt = "基因计数")
  }, deleteFile = FALSE)
  
  output$integration_plot <- renderImage({
    list(src = "images/integration_plot.png",
         contentType = 'image/png',
         width = 600,
         height = 400,
         alt = "整合图")
  }, deleteFile = FALSE)
  
  # --- 数据可用性检查 ---
  output$lineage_data_available <- reactive({
    !is.null(values$lineage_data)
  })
  
  output$single_cell_data_available <- reactive({
    !is.null(values$single_cell_data)
  })
  
  outputOptions(output, "lineage_data_available", suspendWhenHidden = FALSE)
  outputOptions(output, "single_cell_data_available", suspendWhenHidden = FALSE)
  
  # --- 数据导出处理 ---
  output$download_lineage_data <- downloadHandler(
    filename = function() {
      paste("lineage_data_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(values$lineage_data, file)
      values$export_history <- c(values$export_history, 
                                paste(Sys.time(), ": Exported Lineage data"))
    }
  )
  
  output$download_single_cell_data <- downloadHandler(
    filename = function() {
      paste("single_cell_data_", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(values$single_cell_data, file)
      values$export_history <- c(values$export_history, 
                                paste(Sys.time(), ": Exported Single Cell data"))
    }
  )
  
  output$download_qc_summary <- downloadHandler(
    filename = function() {
      paste("qc_summary_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Create mock QC summary data
      qc_data <- data.frame(
        Parameter = c("nFeature_RNA_min", "nFeature_RNA_max", "nCount_RNA_min", 
                     "nCount_RNA_max", "percent.mt_max", "TSS.enrichment_min"),
        Value = c(input$rna_nfeature_min, input$rna_nfeature_max, 
                 input$rna_ncount_min, input$rna_ncount_max, 
                 input$rna_percent_mt_max, input$atac_tss_enrichment_min)
      )
      write.csv(qc_data, file, row.names = FALSE)
      values$export_history <- c(values$export_history, 
                                paste(Sys.time(), ": Exported QC report"))
    }
  )
  
  output$download_cluster_results <- downloadHandler(
    filename = function() {
      paste("cluster_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Create mock clustering results
      cluster_data <- data.frame(
        Cell_ID = paste("Cell", 1:100, sep = "_"),
        Cluster = sample(1:5, 100, replace = TRUE),
        UMAP_1 = rnorm(100),
        UMAP_2 = rnorm(100)
      )
      write.csv(cluster_data, file, row.names = FALSE)
      values$export_history <- c(values$export_history, 
                                paste(Sys.time(), ": Exported clustering results"))
    }
  )
  
  # --- 图表导出处理 ---
  output$download_custom_plot <- downloadHandler(
    filename = function() {
      paste(input$custom_plot_name, "_", input$custom_plot_level, "_", 
            Sys.Date(), ".", input$custom_plot_format, sep = "")
    },
    content = function(file) {
      # 根据选择生成对应的图表
      if(input$custom_plot_format == "png") {
        png(file, width = 800, height = 600)
      } else if(input$custom_plot_format == "pdf") {
        pdf(file, width = 8, height = 6)
      } else if(input$custom_plot_format == "jpeg") {
        jpeg(file, width = 800, height = 600)
      }
      
      # 根据图表类型和级别生成图表
      if(input$custom_plot_type == "rna_cluster") {
        plot(rnorm(100), rnorm(100), main = paste("RNA Cluster -", input$custom_plot_level))
      } else if(input$custom_plot_type == "density") {
        plot(density(rnorm(100)), main = paste("Density Plot -", input$custom_plot_level))
      }
      # 可以根据需要添加更多图表类型
      
      dev.off()
      values$export_history <- c(values$export_history, 
                                paste(Sys.time(), ": Exported custom plot -", input$custom_plot_name))
    }
  )
  
  # --- 导出历史显示 ---
  output$export_history <- renderText({
    if(length(values$export_history) == 0) {
      "No export records"
    } else {
      paste(rev(values$export_history), collapse = "\n")
    }
  })
  
  # --- 清除导出历史 ---
  observeEvent(input$clear_export_history, {
    values$export_history <- c()
  })
}

# 运行App
shinyApp(ui, server)
