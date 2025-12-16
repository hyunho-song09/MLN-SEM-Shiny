################################################################################
# MLN-SEM Shiny Application
# 
# Complete workflow for Multi-Layer Network Structural Equation Modeling
# with Extended Visualizations
################################################################################

library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(ggplot2)
library(plotly)
library(tidyverse)
library(AER)
library(glmnet)
library(ivreg)
library(sandwich)
library(lmtest)
library(lavaan)
library(igraph)

# Source required framework files (will be loaded when app starts)
source("MLN_SEM_v2_1_FINAL.R")
source("MLN_SEM_Visualization_FIXED.R")
source("MLN_SEM_Extended_Network_Visualization_v2.R")
source("MLN_SEM_Module_Assignments_Simple.R")

################################################################################
# UI
################################################################################

ui <- dashboardPage(
  skin = "blue",
  
  # Header
  dashboardHeader(
    title = "MLN-SEM Analysis Platform",
    titleWidth = 300
  ),
  
  # Sidebar
  dashboardSidebar(
    width = 300,
    sidebarMenu(
      menuItem("📁 Data Upload", tabName = "upload", icon = icon("upload")),
      menuItem("⚙️ Parameters", tabName = "parameters", icon = icon("sliders-h")),
      menuItem("🚀 Run Analysis", tabName = "analysis", icon = icon("play-circle")),
      menuItem("📊 Results", tabName = "results", icon = icon("chart-bar")),
      menuItem("🎨 Visualizations", tabName = "visualizations", icon = icon("image")),
      menuItem("📥 Downloads", tabName = "downloads", icon = icon("download")),
      menuItem("ℹ️ Help", tabName = "help", icon = icon("question-circle"))
    )
  ),
  
  # Body
  dashboardBody(
    # Custom CSS
    tags$head(
      tags$style(HTML("
        .box-title { font-weight: bold; }
        .progress-bar { background-color: #3c8dbc; }
        .status-box { padding: 15px; border-radius: 5px; margin-bottom: 10px; }
        .success-box { background-color: #d4edda; border-left: 5px solid #28a745; }
        .warning-box { background-color: #fff3cd; border-left: 5px solid #ffc107; }
        .info-box { background-color: #d1ecf1; border-left: 5px solid #17a2b8; }
      "))
    ),
    
    tabItems(
      # Tab 1: Data Upload
      tabItem(
        tabName = "upload",
        fluidRow(
          box(
            title = "📁 Upload Data Files",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            p("Upload three CSV files for MLN-SEM analysis:"),
            
            fluidRow(
              column(4,
                fileInput("file_clinical", "Clinical Data (01_cp.csv)",
                         accept = c(".csv"), buttonLabel = "Browse...",
                         placeholder = "No file selected")
              ),
              column(4,
                fileInput("file_metabolome", "Metabolome Data (02_met.csv)",
                         accept = c(".csv"), buttonLabel = "Browse...",
                         placeholder = "No file selected")
              ),
              column(4,
                fileInput("file_microbiome", "Microbiome Data (03_micro.csv)",
                         accept = c(".csv"), buttonLabel = "Browse...",
                         placeholder = "No file selected")
              )
            ),
            
            hr(),
            
            h4("📦 Example Data"),
            p("Don't have data yet? Try our example dataset (120 samples, liver disease study):"),
            
            fluidRow(
              column(3,
                div(class = "info-box", style = "padding: 10px;",
                  h5("🔽 Download Examples"),
                  downloadButton("download_example_clinical", "Clinical CSV", 
                               class = "btn-info btn-sm", style = "width: 100%; margin-bottom: 5px;"),
                  downloadButton("download_example_metabolome", "Metabolome CSV",
                               class = "btn-info btn-sm", style = "width: 100%; margin-bottom: 5px;"),
                  downloadButton("download_example_microbiome", "Microbiome CSV",
                               class = "btn-info btn-sm", style = "width: 100%;")
                )
              ),
              column(3,
                div(class = "success-box", style = "padding: 10px;",
                  h5("⚡ Quick Load"),
                  actionButton("load_example_data", "Load Example Data", 
                             icon = icon("bolt"),
                             class = "btn-success btn-lg",
                             style = "width: 100%; height: 80px; font-size: 16px;"),
                  br(), br(),
                  helpText("Instantly load all three files")
                )
              ),
              column(6,
                div(class = "info-box", style = "padding: 10px;",
                  h5("📊 Example Dataset Info"),
                  tags$ul(
                    tags$li(strong("Samples:"), " 120 patients"),
                    tags$li(strong("Clinical:"), " 18 variables (AST, ALT, BMI, etc.)"),
                    tags$li(strong("Metabolome:"), " ~200 metabolites"),
                    tags$li(strong("Microbiome:"), " ~150 microbial species"),
                    tags$li(strong("Study:"), " Liver disease (MASLD/ALD)")
                  )
                )
              )
            ),
            
            hr(),
            
            h4("Data Preview"),
            tabsetPanel(
              id = "data_preview_tabs",
              tabPanel("Clinical", DTOutput("preview_clinical")),
              tabPanel("Metabolome", DTOutput("preview_metabolome")),
              tabPanel("Microbiome", DTOutput("preview_microbiome"))
            )
          )
        ),
        
        fluidRow(
          box(
            title = "📋 Data Summary",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("data_summary")
          )
        )
      ),
      
      # Tab 2: Parameters
      tabItem(
        tabName = "parameters",
        fluidRow(
          box(
            title = "⚙️ Analysis Parameters",
            status = "warning",
            solidHeader = TRUE,
            width = 6,
            
            h4("Outcomes & Covariates"),
            helpText("Select clinical variables to analyze"),
            
            uiOutput("outcome_selector"),
            uiOutput("covariate_selector"),
            
            hr(),
            
            h4("Module Construction"),
            sliderInput("metab_k", "Metabolome Modules (k):",
                       min = 2, max = 20, value = 8, step = 1),
            sliderInput("micro_k", "Microbiome Modules (k):",
                       min = 2, max = 20, value = 8, step = 1),
            sliderInput("min_size", "Minimum Module Size:",
                       min = 2, max = 10, value = 3, step = 1),
            
            hr(),
            
            h4("Quality Control"),
            sliderInput("max_inst", "Max Instruments per Module:",
                       min = 5, max = 30, value = 20, step = 1),
            sliderInput("weak_f", "Weak-F Threshold:",
                       min = 5, max = 20, value = 10, step = 1),
            sliderInput("p_threshold", "P-value Threshold:",
                       min = 0.01, max = 0.10, value = 0.05, step = 0.01),
            
            checkboxInput("use_clr", "Use CLR Transformation", value = TRUE),
            sliderInput("alpha", "Alpha (Ridge penalty):",
                       min = 0, max = 1, value = 1, step = 0.1)
          ),
          
          box(
            title = "💡 Parameter Suggestions",
            status = "success",
            solidHeader = TRUE,
            width = 6,
            
            actionButton("suggest_params", "Get Suggested Parameters",
                        icon = icon("lightbulb"), class = "btn-success btn-lg",
                        style = "width: 100%; margin-bottom: 20px;"),
            
            hr(),
            
            h4("Suggested Values"),
            verbatimTextOutput("suggested_params"),
            
            hr(),
            
            actionButton("apply_suggestions", "Apply Suggestions",
                        icon = icon("check"), class = "btn-info",
                        style = "width: 100%;"),
            
            hr(),
            
            h4("Parameter Guidelines"),
            div(class = "info-box",
              tags$ul(
                tags$li("Metabolome modules: 6-12 recommended"),
                tags$li("Microbiome modules: 6-12 recommended"),
                tags$li("Weak-F > 10 indicates strong instruments"),
                tags$li("Use CLR for compositional data"),
                tags$li("Alpha = 1 for LASSO, 0 for Ridge")
              )
            )
          )
        )
      ),
      
      # Tab 3: Run Analysis
      tabItem(
        tabName = "analysis",
        fluidRow(
          box(
            title = "🚀 Execute MLN-SEM Analysis",
            status = "danger",
            solidHeader = TRUE,
            width = 12,
            
            fluidRow(
              column(6,
                actionButton("run_analysis", "▶ Run Complete Analysis",
                           icon = icon("play-circle"),
                           class = "btn-danger btn-lg",
                           style = "width: 100%; height: 60px; font-size: 18px;")
              ),
              column(6,
                actionButton("stop_analysis", "⏹ Stop Analysis",
                           icon = icon("stop-circle"),
                           class = "btn-warning btn-lg",
                           style = "width: 100%; height: 60px; font-size: 18px;")
              )
            ),
            
            hr(),
            
            h4("Analysis Progress"),
            progressBar(id = "analysis_progress", value = 0, 
                       title = "Waiting to start...",
                       display_pct = TRUE, striped = TRUE, status = "info"),
            
            hr(),
            
            h4("Analysis Log"),
            verbatimTextOutput("analysis_log", placeholder = TRUE),
            
            hr(),
            
            h4("Status"),
            uiOutput("analysis_status")
          )
        )
      ),
      
      # Tab 4: Results
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            title = "📊 Analysis Results Summary",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            
            h4("Module Construction"),
            verbatimTextOutput("module_summary"),
            
            hr(),
            
            h4("IV Analysis Summary"),
            verbatimTextOutput("iv_summary"),
            
            hr(),
            
            h4("Top Validated Pathways"),
            DTOutput("top_pathways")
          )
        ),
        
        fluidRow(
          box(
            title = "📈 All IV Results",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            DTOutput("all_iv_results")
          )
        ),
        
        fluidRow(
          box(
            title = "🎯 Module Assignments",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            collapsible = TRUE,
            collapsed = TRUE,
            
            tabsetPanel(
              tabPanel("Microbiome", DTOutput("microbiome_modules")),
              tabPanel("Metabolome", DTOutput("metabolome_modules"))
            )
          )
        )
      ),
      
      # Tab 5: Visualizations
      tabItem(
        tabName = "visualizations",
        fluidRow(
          box(
            title = "🎨 Network Visualizations",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            
            tabsetPanel(
              id = "viz_tabs",
              
              tabPanel("Validated Pathways Network",
                      br(),
                      plotOutput("plot_validated_network", height = "700px"),
                      hr(),
                      downloadButton("download_validated_pdf", "Download PDF")
              ),
              
              tabPanel("Network Comparison",
                      br(),
                      plotOutput("plot_network_comparison", height = "700px"),
                      hr(),
                      downloadButton("download_comparison_pdf", "Download PDF")
              ),
              
              tabPanel("Module Heatmap",
                      br(),
                      plotOutput("plot_module_heatmap", height = "600px"),
                      hr(),
                      downloadButton("download_heatmap_pdf", "Download PDF")
              ),
              
              tabPanel("Effect Sizes",
                      br(),
                      plotOutput("plot_effect_sizes", height = "600px"),
                      hr(),
                      downloadButton("download_effects_pdf", "Download PDF")
              ),
              
              tabPanel("Weak-F Statistics",
                      br(),
                      plotOutput("plot_weak_f", height = "600px"),
                      hr(),
                      downloadButton("download_weakf_pdf", "Download PDF")
              )
            )
          )
        )
      ),
      
      # Tab 6: Downloads
      tabItem(
        tabName = "downloads",
        fluidRow(
          box(
            title = "📥 Download Results",
            status = "success",
            solidHeader = TRUE,
            width = 12,
            
            h4("Result Files"),
            p("Download analysis results in various formats:"),
            
            fluidRow(
              column(3,
                div(class = "success-box",
                  h4("Complete Results"),
                  downloadButton("download_rds", "Download RDS", 
                               class = "btn-success", style = "width: 100%;"),
                  br(), br(),
                  helpText("Contains all analysis objects")
                )
              ),
              column(3,
                div(class = "success-box",
                  h4("Passed Pathways"),
                  downloadButton("download_pathways_csv", "Download CSV",
                               class = "btn-success", style = "width: 100%;"),
                  br(), br(),
                  helpText("Validated pathways table")
                )
              ),
              column(3,
                div(class = "success-box",
                  h4("All IV Results"),
                  downloadButton("download_iv_csv", "Download CSV",
                               class = "btn-success", style = "width: 100%;"),
                  br(), br(),
                  helpText("Complete IV analysis results")
                )
              ),
              column(3,
                div(class = "success-box",
                  h4("Module Assignments"),
                  downloadButton("download_modules_xlsx", "Download XLSX",
                               class = "btn-success", style = "width: 100%;"),
                  br(), br(),
                  helpText("Feature-to-module mapping")
                )
              )
            ),
            
            hr(),
            
            h4("All Visualizations"),
            fluidRow(
              column(6,
                downloadButton("download_all_plots", "Download All Plots (ZIP)",
                             class = "btn-primary btn-lg", 
                             style = "width: 100%; height: 50px; font-size: 16px;")
              ),
              column(6,
                downloadButton("download_report", "Generate HTML Report",
                             class = "btn-info btn-lg",
                             style = "width: 100%; height: 50px; font-size: 16px;")
              )
            )
          )
        )
      ),
      
      # Tab 7: Help
      tabItem(
        tabName = "help",
        fluidRow(
          box(
            title = "ℹ️ User Guide",
            status = "info",
            solidHeader = TRUE,
            width = 12,
            
            h3("MLN-SEM Analysis Workflow"),
            
            div(class = "info-box",
              h4("1. Data Upload"),
              p("Upload three CSV files:"),
              tags$ul(
                tags$li(strong("Clinical data (01_cp.csv):"), " Patient phenotype and clinical parameters"),
                tags$li(strong("Metabolome data (02_met.csv):"), " Metabolite abundance data"),
                tags$li(strong("Microbiome data (03_micro.csv):"), " Microbial abundance data")
              ),
              p("All files must have an 'Index' column for sample matching.")
            ),
            
            hr(),
            
            div(class = "info-box",
              h4("2. Parameter Configuration"),
              p("Configure analysis parameters or use suggested values:"),
              tags$ul(
                tags$li(strong("Outcomes:"), " Clinical variables to predict"),
                tags$li(strong("Covariates:"), " Variables to adjust for (Age, Gender, BMI)"),
                tags$li(strong("Module counts:"), " Number of modules for dimension reduction"),
                tags$li(strong("Quality thresholds:"), " Weak-F and p-value cutoffs")
              )
            ),
            
            hr(),
            
            div(class = "info-box",
              h4("3. Run Analysis"),
              p("Execute the complete MLN-SEM pipeline:"),
              tags$ul(
                tags$li("Data preprocessing and CLR transformation"),
                tags$li("Module construction using eigengene-based clustering"),
                tags$li("Instrumental variable analysis"),
                tags$li("Pathway validation with Sargan tests"),
                tags$li("Network visualization generation")
              ),
              p(tags$em("Note: Analysis may take several minutes depending on data size."))
            ),
            
            hr(),
            
            div(class = "success-box",
              h4("4. Explore Results"),
              p("View and download:"),
              tags$ul(
                tags$li("Module summaries and assignments"),
                tags$li("Validated causal pathways"),
                tags$li("Network visualizations"),
                tags$li("Statistical tables and effect sizes")
              )
            ),
            
            hr(),
            
            h4("Reference"),
            p("For detailed methodology, see:"),
            tags$ul(
              tags$li("Multi-Layer Network SEM framework documentation"),
              tags$li("WGCNA module construction methods"),
              tags$li("Instrumental variable regression theory")
            ),
            
            hr(),
            
            h4("Contact"),
            p("For questions or issues, contact: HHSONG @ Seoul National University"),
            p("Metabolomics Lab, Prof. Do Yup Lee")
          )
        )
      )
    )
  )
)

################################################################################
# SERVER
################################################################################

server <- function(input, output, session) {
  
  # Reactive values to store data and results
  rv <- reactiveValues(
    clinical = NULL,
    metabolome = NULL,
    microbiome = NULL,
    results = NULL,
    module_assignments = NULL,
    analysis_running = FALSE,
    log_text = ""
  )
  
  # Helper function to append log messages
  add_log <- function(message) {
    timestamp <- format(Sys.time(), "[%H:%M:%S]")
    rv$log_text <- paste0(rv$log_text, timestamp, " ", message, "\n")
  }
  
  ############################################################################
  # DATA UPLOAD
  ############################################################################
  
  # Load clinical data
  observeEvent(input$file_clinical, {
    req(input$file_clinical)
    tryCatch({
      rv$clinical <- read.csv(input$file_clinical$datapath, stringsAsFactors = FALSE)
      add_log("✓ Clinical data loaded successfully")
      showNotification("Clinical data loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading clinical data:", e$message), type = "error")
    })
  })
  
  # Load metabolome data
  observeEvent(input$file_metabolome, {
    req(input$file_metabolome)
    tryCatch({
      rv$metabolome <- read.csv(input$file_metabolome$datapath, stringsAsFactors = FALSE)
      add_log("✓ Metabolome data loaded successfully")
      showNotification("Metabolome data loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading metabolome data:", e$message), type = "error")
    })
  })
  
  # Load microbiome data
  observeEvent(input$file_microbiome, {
    req(input$file_microbiome)
    tryCatch({
      rv$microbiome <- read.csv(input$file_microbiome$datapath, stringsAsFactors = FALSE)
      add_log("✓ Microbiome data loaded successfully")
      showNotification("Microbiome data loaded!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading microbiome data:", e$message), type = "error")
    })
  })
  
  # Data previews
  output$preview_clinical <- renderDT({
    req(rv$clinical)
    datatable(head(rv$clinical, 100), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$preview_metabolome <- renderDT({
    req(rv$metabolome)
    datatable(head(rv$metabolome, 100), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$preview_microbiome <- renderDT({
    req(rv$microbiome)
    datatable(head(rv$microbiome, 100), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  ############################################################################
  # EXAMPLE DATA DOWNLOAD
  ############################################################################
  
  # Download example clinical data
  output$download_example_clinical <- downloadHandler(
    filename = function() {
      "example_01_cp.csv"
    },
    content = function(file) {
      # Use embedded CSV file in app directory
      file.copy("example_clinical.csv", file)
    }
  )
  
  # Download example metabolome data
  output$download_example_metabolome <- downloadHandler(
    filename = function() {
      "example_02_met.csv"
    },
    content = function(file) {
      file.copy("example_metabolome.csv", file)
    }
  )
  
  # Download example microbiome data
  output$download_example_microbiome <- downloadHandler(
    filename = function() {
      "example_03_micro.csv"
    },
    content = function(file) {
      file.copy("example_microbiome.csv", file)
    }
  )
  
  # Load example data with one click
  observeEvent(input$load_example_data, {
    
    withProgress(message = "Loading example data...", value = 0, {
      
      tryCatch({
        # Load from embedded CSV files in app directory
        
        # Load clinical data
        incProgress(0.2, detail = "Loading clinical data...")
        rv$clinical <- read.csv("example_clinical.csv", stringsAsFactors = FALSE)
        add_log("✓ Example clinical data loaded (120 samples, 18 variables)")
        
        # Load metabolome data
        incProgress(0.4, detail = "Loading metabolome data...")
        rv$metabolome <- read.csv("example_metabolome.csv", stringsAsFactors = FALSE)
        metab_count <- length(grep("^M_", colnames(rv$metabolome)))
        add_log(sprintf("✓ Example metabolome data loaded (120 samples, %d metabolites)", metab_count))
        
        # Load microbiome data
        incProgress(0.6, detail = "Loading microbiome data...")
        rv$microbiome <- read.csv("example_microbiome.csv", stringsAsFactors = FALSE)
        micro_count <- length(grep("^s_", colnames(rv$microbiome)))
        add_log(sprintf("✓ Example microbiome data loaded (120 samples, %d microbes)", micro_count))
        
        incProgress(1.0, detail = "Complete!")
        
        showNotification(
          "Example data loaded successfully! You can now proceed to Parameters tab.",
          type = "message",
          duration = 5
        )
        
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        add_log("📦 Example dataset ready for analysis")
        add_log("   Next: Go to Parameters tab → Get Suggested Parameters")
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        
      }, error = function(e) {
        showNotification(
          paste("Error loading example data:", e$message, 
                "\nMake sure example CSV files are in the app directory."),
          type = "error",
          duration = NULL
        )
        add_log(sprintf("✗ Error loading example data: %s", e$message))
        add_log("  → Check if example_clinical.csv, example_metabolome.csv, example_microbiome.csv exist")
      })
    })
  })
  
  # Data summary
  output$data_summary <- renderText({
    if (is.null(rv$clinical) && is.null(rv$metabolome) && is.null(rv$microbiome)) {
      return("No data loaded yet. Please upload CSV files.")
    }
    
    summary_text <- ""
    
    if (!is.null(rv$clinical)) {
      summary_text <- paste0(summary_text, 
        sprintf("Clinical Data:\n  - Samples: %d\n  - Variables: %d\n\n",
                nrow(rv$clinical), ncol(rv$clinical)))
    }
    
    if (!is.null(rv$metabolome)) {
      metab_features <- length(grep("^M_", colnames(rv$metabolome)))
      summary_text <- paste0(summary_text,
        sprintf("Metabolome Data:\n  - Samples: %d\n  - Total columns: %d\n  - Metabolite features: %d\n\n",
                nrow(rv$metabolome), ncol(rv$metabolome), metab_features))
    }
    
    if (!is.null(rv$microbiome)) {
      micro_features <- length(grep("^s_", colnames(rv$microbiome)))
      summary_text <- paste0(summary_text,
        sprintf("Microbiome Data:\n  - Samples: %d\n  - Total columns: %d\n  - Microbial features: %d\n\n",
                nrow(rv$microbiome), ncol(rv$microbiome), micro_features))
    }
    
    return(summary_text)
  })
  
  ############################################################################
  # PARAMETERS
  ############################################################################
  
  # Dynamic outcome selector
  output$outcome_selector <- renderUI({
    req(rv$clinical)
    
    # Find potential outcome columns (after Age, Gender, BMI, etc.)
    all_cols <- colnames(rv$clinical)
    basic_cols <- c("Index", "Sample", "Group", "Gender", "Age", "BMI")
    outcome_cols <- setdiff(all_cols, basic_cols)
    
    checkboxGroupInput("outcomes", "Select Outcome Variables:",
                      choices = outcome_cols,
                      selected = head(outcome_cols, min(5, length(outcome_cols))))
  })
  
  # Dynamic covariate selector
  output$covariate_selector <- renderUI({
    req(rv$clinical)
    
    covariate_options <- c("Age", "Gender", "BMI")
    available_covariates <- covariate_options[covariate_options %in% colnames(rv$clinical)]
    
    checkboxGroupInput("covariates", "Select Covariates:",
                      choices = available_covariates,
                      selected = available_covariates)
  })
  
  # Suggest parameters
  observeEvent(input$suggest_params, {
    req(rv$clinical, rv$metabolome, rv$microbiome)
    
    withProgress(message = "Calculating suggested parameters...", value = 0.5, {
      
      params <- suggest_parameters(
        n = nrow(rv$clinical),
        p_metab = length(grep("^M_", colnames(rv$metabolome))),
        p_micro = length(grep("^s_", colnames(rv$microbiome)))
      )
      
      output$suggested_params <- renderText({
        paste0(
          sprintf("Based on your data dimensions:\n\n"),
          sprintf("Metabolome modules (k): %d\n", params$metab_k),
          sprintf("Microbiome modules (k): %d\n", params$micro_k),
          sprintf("Minimum module size: %d\n", params$min_size),
          sprintf("Max instruments: %d\n", params$max_inst),
          sprintf("Weak-F threshold: %d\n\n", params$weak_f_threshold),
          sprintf("These values are optimized for:\n"),
          sprintf("  - %d samples\n", nrow(rv$clinical)),
          sprintf("  - %d metabolites\n", length(grep("^M_", colnames(rv$metabolome)))),
          sprintf("  - %d microbes\n", length(grep("^s_", colnames(rv$microbiome))))
        )
      })
      
      rv$suggested_params <- params
    })
  })
  
  # Apply suggestions
  observeEvent(input$apply_suggestions, {
    req(rv$suggested_params)
    
    updateSliderInput(session, "metab_k", value = rv$suggested_params$metab_k)
    updateSliderInput(session, "micro_k", value = rv$suggested_params$micro_k)
    updateSliderInput(session, "min_size", value = rv$suggested_params$min_size)
    updateSliderInput(session, "max_inst", value = rv$suggested_params$max_inst)
    updateSliderInput(session, "weak_f", value = rv$suggested_params$weak_f_threshold)
    
    showNotification("Parameters applied!", type = "message")
  })
  
  ############################################################################
  # RUN ANALYSIS
  ############################################################################
  
  observeEvent(input$run_analysis, {
    req(rv$clinical, rv$metabolome, rv$microbiome)
    req(input$outcomes, input$covariates)
    
    # Prevent multiple runs
    if (rv$analysis_running) {
      showNotification("Analysis is already running!", type = "warning")
      return()
    }
    
    rv$analysis_running <- TRUE
    rv$log_text <- ""
    
    withProgress(message = "Running MLN-SEM Analysis", value = 0, {
      
      tryCatch({
        # Step 1: Preprocessing
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        add_log("STEP 1: Preprocessing data...")
        updateProgressBar(session, "analysis_progress", value = 10, title = "Preprocessing...")
        
        pheno_liver <- rv$clinical
        metab_liver <- rv$metabolome
        micro_liver <- rv$microbiome
        
        # Clean column names
        colnames(pheno_liver) <- make.names(colnames(pheno_liver))
        
        # Fix AST/ALT if present
        if("AST.ALT" %in% colnames(pheno_liver)) {
          colnames(pheno_liver)[colnames(pheno_liver) == "AST.ALT"] <- "AST_ALT"
        }
        
        # Handle Gender
        if(length(unique(pheno_liver$Gender)) > 1) {
          pheno_liver$Gender <- factor(pheno_liver$Gender)
          add_log("  ✓ Gender converted to factor")
        } else {
          pheno_liver$Gender <- as.numeric(pheno_liver$Gender)
          add_log("  ✓ Gender converted to numeric")
        }
        
        add_log(sprintf("  ✓ Samples: %d", nrow(pheno_liver)))
        add_log(sprintf("  ✓ Metabolites: %d", length(grep("^M_", colnames(metab_liver)))))
        add_log(sprintf("  ✓ Microbes: %d", length(grep("^s_", colnames(micro_liver)))))
        
        # Step 2: Run Analysis
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        add_log("STEP 2: Running MLN-SEM analysis...")
        add_log("  (This may take several minutes)")
        updateProgressBar(session, "analysis_progress", value = 20, title = "Running MLN-SEM...")
        
        results <- run_mln_sem(
          pheno = pheno_liver,
          metabolome = metab_liver,
          microbiome = micro_liver,
          outcomes = input$outcomes,
          covariates = input$covariates,
          metab_k = input$metab_k,
          micro_k = input$micro_k,
          min_size = input$min_size,
          max_inst = input$max_inst,
          weak_f_threshold = input$weak_f,
          p_threshold = input$p_threshold,
          sargan_threshold = 0.05,
          use_clr = input$use_clr,
          alpha = input$alpha
        )
        
        add_log("  ✓ MLN-SEM analysis complete!")
        updateProgressBar(session, "analysis_progress", value = 60, title = "Generating visualizations...")
        
        # Step 3: Generate visualizations
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        add_log("STEP 3: Generating visualizations...")
        
        # Create temp directory for plots
        temp_dir <- tempdir()
        
        extended_plots <- visualize_extended_networks(
          results = results,
          pheno_data = pheno_liver,
          output_dir = temp_dir,
          prefix = "analysis"
        )
        
        add_log("  ✓ Extended network visualizations created")
        updateProgressBar(session, "analysis_progress", value = 80, title = "Extracting modules...")
        
        # Step 4: Extract module assignments
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        add_log("STEP 4: Extracting module assignments...")
        
        module_assignments <- extract_and_save_simple_modules(
          results = results,
          microbiome = micro_liver,
          metabolome = metab_liver,
          output_dir = temp_dir,
          prefix = "analysis"
        )
        
        add_log("  ✓ Module assignments extracted")
        
        # Store results
        rv$results <- results
        rv$module_assignments <- module_assignments
        rv$pheno_data <- pheno_liver
        
        updateProgressBar(session, "analysis_progress", value = 100, 
                         title = "Analysis Complete!", status = "success")
        
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        add_log("✓ ANALYSIS COMPLETE!")
        add_log(sprintf("  • Metabolome modules: %d", ncol(results$metabolome_modules$MEs)))
        add_log(sprintf("  • Microbiome modules: %d", ncol(results$microbiome_modules$MEs)))
        add_log(sprintf("  • Pathways tested: %d", nrow(results$iv_table)))
        add_log(sprintf("  • Pathways validated: %d", nrow(results$passed_pathways)))
        add_log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
        
        showNotification("Analysis completed successfully!", type = "message", duration = 5)
        
        # Switch to results tab
        updateTabItems(session, "tabs", "results")
        
      }, error = function(e) {
        add_log(sprintf("✗ ERROR: %s", e$message))
        updateProgressBar(session, "analysis_progress", value = 0, 
                         title = "Error occurred", status = "danger")
        showNotification(paste("Analysis failed:", e$message), type = "error", duration = NULL)
      }, finally = {
        rv$analysis_running <- FALSE
      })
    })
  })
  
  # Stop analysis
  observeEvent(input$stop_analysis, {
    if (rv$analysis_running) {
      rv$analysis_running <- FALSE
      add_log("⚠ Analysis stopped by user")
      updateProgressBar(session, "analysis_progress", value = 0, 
                       title = "Stopped", status = "warning")
      showNotification("Analysis stopped", type = "warning")
    }
  })
  
  # Analysis log output
  output$analysis_log <- renderText({
    rv$log_text
  })
  
  # Analysis status
  output$analysis_status <- renderUI({
    if (is.null(rv$results)) {
      div(class = "warning-box",
        h4("⚠ No Results Available"),
        p("Please run the analysis first.")
      )
    } else {
      div(class = "success-box",
        h4("✓ Analysis Complete"),
        p(sprintf("Validated %d pathways from %d tested",
                 nrow(rv$results$passed_pathways),
                 nrow(rv$results$iv_table)))
      )
    }
  })
  
  ############################################################################
  # RESULTS
  ############################################################################
  
  # Module summary
  output$module_summary <- renderText({
    req(rv$results)
    
    paste0(
      sprintf("Metabolome Modules: %d\n", ncol(rv$results$metabolome_modules$MEs)),
      sprintf("  Modules: %s\n\n", 
              paste(colnames(rv$results$metabolome_modules$MEs), collapse = ", ")),
      sprintf("Microbiome Modules: %d\n", ncol(rv$results$microbiome_modules$MEs)),
      sprintf("  Modules: %s\n",
              paste(colnames(rv$results$microbiome_modules$MEs), collapse = ", "))
    )
  })
  
  # IV summary
  output$iv_summary <- renderText({
    req(rv$results)
    
    iv_table <- rv$results$iv_table
    
    if (nrow(iv_table) == 0) {
      return("No IV analysis results available.")
    }
    
    weakF_vals <- iv_table$weak_F[!is.na(iv_table$weak_F)]
    sig_pathways <- sum(iv_table$p_value < 0.05, na.rm = TRUE)
    
    paste0(
      sprintf("Total Pathways Tested: %d\n\n", nrow(iv_table)),
      sprintf("Weak-F Statistics:\n"),
      sprintf("  Range: %.2f - %.2f\n", min(weakF_vals), max(weakF_vals)),
      sprintf("  Strong instruments (F>10): %d\n\n", sum(weakF_vals > 10)),
      sprintf("Significant Pathways (p<0.05): %d\n", sig_pathways)
    )
  })
  
  # Top pathways table
  output$top_pathways <- renderDT({
    req(rv$results)
    
    if (nrow(rv$results$passed_pathways) == 0) {
      return(data.frame(Message = "No validated pathways"))
    }
    
    top_pw <- rv$results$passed_pathways %>%
      select(module, outcome, beta, se, p_value, weak_F, sargan_p) %>%
      arrange(p_value) %>%
      head(20)
    
    datatable(top_pw, 
             options = list(scrollX = TRUE, pageLength = 10),
             rownames = FALSE) %>%
      formatRound(columns = c('beta', 'se', 'weak_F'), digits = 3) %>%
      formatSignif(columns = c('p_value', 'sargan_p'), digits = 3)
  })
  
  # All IV results
  output$all_iv_results <- renderDT({
    req(rv$results)
    
    datatable(rv$results$iv_table,
             options = list(scrollX = TRUE, pageLength = 25),
             rownames = FALSE,
             filter = 'top') %>%
      formatRound(columns = c('beta', 'se', 'weak_F'), digits = 3) %>%
      formatSignif(columns = c('p_value', 'sargan_p'), digits = 3)
  })
  
  # Module assignments
  output$microbiome_modules <- renderDT({
    req(rv$module_assignments)
    
    datatable(rv$module_assignments$microbiome,
             options = list(scrollX = TRUE, pageLength = 25),
             rownames = FALSE,
             filter = 'top')
  })
  
  output$metabolome_modules <- renderDT({
    req(rv$module_assignments)
    
    datatable(rv$module_assignments$metabolome,
             options = list(scrollX = TRUE, pageLength = 25),
             rownames = FALSE,
             filter = 'top')
  })
  
  ############################################################################
  # VISUALIZATIONS
  ############################################################################
  
  # Validated pathways network
  output$plot_validated_network <- renderPlot({
    req(rv$results, rv$pheno_data)
    
    plot_validated_pathways_network(
      results = rv$results,
      pheno_data = rv$pheno_data,
      layout = "fr",
      show_labels = TRUE,
      output_file = NULL
    )
  })
  
  # Network comparison
  output$plot_network_comparison <- renderPlot({
    req(rv$results, rv$pheno_data)
    
    plot_network_comparison(
      results = rv$results,
      pheno_data = rv$pheno_data,
      cor_threshold = 0.3,
      correlation_method = "spearman",
      layout = "fr",
      show_labels = TRUE,
      output_file = NULL
    )
  })
  
  # Module heatmap (if visualize_all creates one)
  output$plot_module_heatmap <- renderPlot({
    req(rv$results)
    
    # Create a simple heatmap of module eigengenes
    MEs <- rv$results$metabolome_modules$MEs
    
    # Correlation heatmap
    cor_matrix <- cor(MEs, use = "pairwise.complete.obs")
    
    library(pheatmap)
    pheatmap::pheatmap(cor_matrix,
                      main = "Metabolome Module Correlation",
                      color = colorRampPalette(c("blue", "white", "red"))(100),
                      breaks = seq(-1, 1, length.out = 101),
                      display_numbers = TRUE,
                      number_format = "%.2f",
                      fontsize = 10)
  })
  
  # Effect sizes plot
  output$plot_effect_sizes <- renderPlot({
    req(rv$results)
    
    if (nrow(rv$results$passed_pathways) == 0) {
      plot.new()
      text(0.5, 0.5, "No validated pathways to display", cex = 1.5)
      return()
    }
    
    pw_data <- rv$results$passed_pathways %>%
      arrange(desc(abs(beta))) %>%
      head(20) %>%
      mutate(pathway = paste(module, "→", outcome),
             pathway = factor(pathway, levels = pathway))
    
    ggplot(pw_data, aes(x = beta, y = pathway)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      geom_errorbarh(aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se),
                     height = 0.3, color = "gray40") +
      geom_point(aes(size = weak_F, color = p_value < 0.01), alpha = 0.8) +
      scale_color_manual(values = c("TRUE" = "red", "FALSE" = "orange"),
                        name = "p < 0.01") +
      scale_size_continuous(name = "Weak-F", range = c(3, 10)) +
      labs(title = "Top 20 Pathway Effect Sizes",
           subtitle = "Error bars show 95% CI",
           x = "Effect Size (β)", y = "") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "right")
  })
  
  # Weak-F statistics plot
  output$plot_weak_f <- renderPlot({
    req(rv$results)
    
    if (nrow(rv$results$passed_pathways) == 0) {
      plot.new()
      text(0.5, 0.5, "No validated pathways to display", cex = 1.5)
      return()
    }
    
    pw_data <- rv$results$passed_pathways %>%
      arrange(desc(weak_F)) %>%
      head(30) %>%
      mutate(pathway = paste(module, "→", outcome),
             pathway = factor(pathway, levels = rev(pathway)),
             strength = case_when(
               weak_F > 20 ~ "Very Strong (F > 20)",
               weak_F > 10 ~ "Strong (F > 10)",
               TRUE ~ "Moderate (F < 10)"
             ))
    
    ggplot(pw_data, aes(x = weak_F, y = pathway, fill = strength)) +
      geom_col(alpha = 0.8) +
      geom_vline(xintercept = 10, linetype = "dashed", color = "red", size = 1) +
      scale_fill_manual(values = c("Very Strong (F > 20)" = "#2E7D32",
                                   "Strong (F > 10)" = "#66BB6A",
                                   "Moderate (F < 10)" = "#FFA726")) +
      labs(title = "Instrumental Variable Strength (Weak-F Statistics)",
           subtitle = "Dashed line indicates F = 10 threshold",
           x = "Weak-F Statistic", y = "",
           fill = "Instrument Strength") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")
  })
  
  ############################################################################
  # DOWNLOADS
  ############################################################################
  
  # Download RDS
  output$download_rds <- downloadHandler(
    filename = function() {
      paste0("MLN_SEM_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds")
    },
    content = function(file) {
      req(rv$results)
      saveRDS(rv$results, file)
    }
  )
  
  # Download passed pathways CSV
  output$download_pathways_csv <- downloadHandler(
    filename = function() {
      paste0("Passed_Pathways_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(rv$results)
      write.csv(rv$results$passed_pathways, file, row.names = FALSE)
    }
  )
  
  # Download all IV results CSV
  output$download_iv_csv <- downloadHandler(
    filename = function() {
      paste0("All_IV_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      req(rv$results)
      write.csv(rv$results$iv_table, file, row.names = FALSE)
    }
  )
  
  # Download module assignments (would need writexl package)
  output$download_modules_xlsx <- downloadHandler(
    filename = function() {
      paste0("Module_Assignments_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
    },
    content = function(file) {
      req(rv$module_assignments)
      
      # Create workbook
      library(writexl)
      
      sheets <- list(
        "Microbiome_Modules" = rv$module_assignments$microbiome,
        "Metabolome_Modules" = rv$module_assignments$metabolome
      )
      
      write_xlsx(sheets, file)
    }
  )
  
  # Download individual plots
  output$download_validated_pdf <- downloadHandler(
    filename = function() {
      paste0("Validated_Network_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      req(rv$results, rv$pheno_data)
      
      pdf(file, width = 12, height = 10)
      plot_validated_pathways_network(
        results = rv$results,
        pheno_data = rv$pheno_data,
        layout = "fr",
        show_labels = TRUE,
        output_file = NULL
      )
      dev.off()
    }
  )
  
  output$download_comparison_pdf <- downloadHandler(
    filename = function() {
      paste0("Network_Comparison_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      req(rv$results, rv$pheno_data)
      
      pdf(file, width = 16, height = 8)
      plot_network_comparison(
        results = rv$results,
        pheno_data = rv$pheno_data,
        cor_threshold = 0.3,
        correlation_method = "spearman",
        layout = "fr",
        show_labels = TRUE,
        output_file = NULL
      )
      dev.off()
    }
  )
  
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() {
      paste0("Module_Heatmap_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      req(rv$results)
      
      pdf(file, width = 10, height = 10)
      MEs <- rv$results$metabolome_modules$MEs
      cor_matrix <- cor(MEs, use = "pairwise.complete.obs")
      
      pheatmap::pheatmap(cor_matrix,
                        main = "Metabolome Module Correlation",
                        color = colorRampPalette(c("blue", "white", "red"))(100),
                        breaks = seq(-1, 1, length.out = 101),
                        display_numbers = TRUE,
                        number_format = "%.2f",
                        fontsize = 10)
      dev.off()
    }
  )
  
  output$download_effects_pdf <- downloadHandler(
    filename = function() {
      paste0("Effect_Sizes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      req(rv$results)
      
      if (nrow(rv$results$passed_pathways) == 0) {
        return()
      }
      
      pw_data <- rv$results$passed_pathways %>%
        arrange(desc(abs(beta))) %>%
        head(20) %>%
        mutate(pathway = paste(module, "→", outcome),
               pathway = factor(pathway, levels = pathway))
      
      p <- ggplot(pw_data, aes(x = beta, y = pathway)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
        geom_errorbarh(aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se),
                       height = 0.3, color = "gray40") +
        geom_point(aes(size = weak_F, color = p_value < 0.01), alpha = 0.8) +
        scale_color_manual(values = c("TRUE" = "red", "FALSE" = "orange"),
                          name = "p < 0.01") +
        scale_size_continuous(name = "Weak-F", range = c(3, 10)) +
        labs(title = "Top 20 Pathway Effect Sizes",
             subtitle = "Error bars show 95% CI",
             x = "Effect Size (β)", y = "") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "right")
      
      ggsave(file, p, width = 10, height = 8)
    }
  )
  
  output$download_weakf_pdf <- downloadHandler(
    filename = function() {
      paste0("WeakF_Statistics_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file) {
      req(rv$results)
      
      if (nrow(rv$results$passed_pathways) == 0) {
        return()
      }
      
      pw_data <- rv$results$passed_pathways %>%
        arrange(desc(weak_F)) %>%
        head(30) %>%
        mutate(pathway = paste(module, "→", outcome),
               pathway = factor(pathway, levels = rev(pathway)),
               strength = case_when(
                 weak_F > 20 ~ "Very Strong (F > 20)",
                 weak_F > 10 ~ "Strong (F > 10)",
                 TRUE ~ "Moderate (F < 10)"
               ))
      
      p <- ggplot(pw_data, aes(x = weak_F, y = pathway, fill = strength)) +
        geom_col(alpha = 0.8) +
        geom_vline(xintercept = 10, linetype = "dashed", color = "red", size = 1) +
        scale_fill_manual(values = c("Very Strong (F > 20)" = "#2E7D32",
                                     "Strong (F > 10)" = "#66BB6A",
                                     "Moderate (F < 10)" = "#FFA726")) +
        labs(title = "Instrumental Variable Strength (Weak-F Statistics)",
             subtitle = "Dashed line indicates F = 10 threshold",
             x = "Weak-F Statistic", y = "",
             fill = "Instrument Strength") +
        theme_minimal(base_size = 12) +
        theme(legend.position = "bottom")
      
      ggsave(file, p, width = 10, height = 10)
    }
  )
  
  # Download all plots as ZIP
  output$download_all_plots <- downloadHandler(
    filename = function() {
      paste0("All_Plots_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".zip")
    },
    content = function(file) {
      req(rv$results, rv$pheno_data)
      
      # Create temp directory
      temp_dir <- tempdir()
      plot_files <- c()
      
      # Generate all plots
      tryCatch({
        # Validated network
        f1 <- file.path(temp_dir, "validated_network.pdf")
        pdf(f1, width = 12, height = 10)
        plot_validated_pathways_network(rv$results, rv$pheno_data, "fr", TRUE, NULL)
        dev.off()
        plot_files <- c(plot_files, f1)
        
        # Network comparison
        f2 <- file.path(temp_dir, "network_comparison.pdf")
        pdf(f2, width = 16, height = 8)
        plot_network_comparison(rv$results, rv$pheno_data, 0.3, "spearman", "fr", TRUE, NULL)
        dev.off()
        plot_files <- c(plot_files, f2)
        
        # Heatmap
        f3 <- file.path(temp_dir, "module_heatmap.pdf")
        pdf(f3, width = 10, height = 10)
        MEs <- rv$results$metabolome_modules$MEs
        cor_matrix <- cor(MEs, use = "pairwise.complete.obs")
        pheatmap::pheatmap(cor_matrix, main = "Metabolome Module Correlation",
                          color = colorRampPalette(c("blue", "white", "red"))(100),
                          breaks = seq(-1, 1, length.out = 101),
                          display_numbers = TRUE, number_format = "%.2f", fontsize = 10)
        dev.off()
        plot_files <- c(plot_files, f3)
        
        # Effect sizes (if pathways exist)
        if (nrow(rv$results$passed_pathways) > 0) {
          f4 <- file.path(temp_dir, "effect_sizes.pdf")
          pw_data <- rv$results$passed_pathways %>%
            arrange(desc(abs(beta))) %>% head(20) %>%
            mutate(pathway = paste(module, "→", outcome),
                   pathway = factor(pathway, levels = pathway))
          
          p <- ggplot(pw_data, aes(x = beta, y = pathway)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
            geom_errorbarh(aes(xmin = beta - 1.96*se, xmax = beta + 1.96*se),
                           height = 0.3, color = "gray40") +
            geom_point(aes(size = weak_F, color = p_value < 0.01), alpha = 0.8) +
            scale_color_manual(values = c("TRUE" = "red", "FALSE" = "orange"), name = "p < 0.01") +
            scale_size_continuous(name = "Weak-F", range = c(3, 10)) +
            labs(title = "Top 20 Pathway Effect Sizes", subtitle = "Error bars show 95% CI",
                 x = "Effect Size (β)", y = "") +
            theme_minimal(base_size = 12) + theme(legend.position = "right")
          
          ggsave(f4, p, width = 10, height = 8)
          plot_files <- c(plot_files, f4)
          
          # Weak-F
          f5 <- file.path(temp_dir, "weakf_statistics.pdf")
          pw_data2 <- rv$results$passed_pathways %>%
            arrange(desc(weak_F)) %>% head(30) %>%
            mutate(pathway = paste(module, "→", outcome),
                   pathway = factor(pathway, levels = rev(pathway)),
                   strength = case_when(weak_F > 20 ~ "Very Strong (F > 20)",
                                       weak_F > 10 ~ "Strong (F > 10)",
                                       TRUE ~ "Moderate (F < 10)"))
          
          p2 <- ggplot(pw_data2, aes(x = weak_F, y = pathway, fill = strength)) +
            geom_col(alpha = 0.8) +
            geom_vline(xintercept = 10, linetype = "dashed", color = "red", size = 1) +
            scale_fill_manual(values = c("Very Strong (F > 20)" = "#2E7D32",
                                        "Strong (F > 10)" = "#66BB6A",
                                        "Moderate (F < 10)" = "#FFA726")) +
            labs(title = "Instrumental Variable Strength (Weak-F Statistics)",
                 subtitle = "Dashed line indicates F = 10 threshold",
                 x = "Weak-F Statistic", y = "", fill = "Instrument Strength") +
            theme_minimal(base_size = 12) + theme(legend.position = "bottom")
          
          ggsave(f5, p2, width = 10, height = 10)
          plot_files <- c(plot_files, f5)
        }
        
        # Create ZIP
        zip(file, plot_files, flags = "-j")
        
      }, error = function(e) {
        showNotification(paste("Error creating plots:", e$message), type = "error")
      })
    }
  )
  
  # Generate HTML report (simplified)
  output$download_report <- downloadHandler(
    filename = function() {
      paste0("MLN_SEM_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html")
    },
    content = function(file) {
      req(rv$results)
      
      # Create simple HTML report
      html_content <- sprintf("
<!DOCTYPE html>
<html>
<head>
  <title>MLN-SEM Analysis Report</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 40px; }
    h1 { color: #2E7D32; }
    h2 { color: #1976D2; margin-top: 30px; }
    table { border-collapse: collapse; width: 100%%; margin-top: 20px; }
    th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
    th { background-color: #4CAF50; color: white; }
    .summary { background-color: #f0f0f0; padding: 15px; border-radius: 5px; }
  </style>
</head>
<body>
  <h1>MLN-SEM Analysis Report</h1>
  <p><strong>Generated:</strong> %s</p>
  
  <div class='summary'>
    <h2>Analysis Summary</h2>
    <p><strong>Metabolome Modules:</strong> %d</p>
    <p><strong>Microbiome Modules:</strong> %d</p>
    <p><strong>Pathways Tested:</strong> %d</p>
    <p><strong>Pathways Validated:</strong> %d</p>
  </div>
  
  <h2>Top Validated Pathways</h2>
  <table>
    <tr>
      <th>Module</th>
      <th>Outcome</th>
      <th>Beta</th>
      <th>SE</th>
      <th>P-value</th>
      <th>Weak-F</th>
    </tr>
    %s
  </table>
  
  <h2>Notes</h2>
  <p>This report was automatically generated by the MLN-SEM Shiny Application.</p>
  <p>For detailed visualizations and data files, please download from the Downloads tab.</p>
</body>
</html>
      ",
      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      ncol(rv$results$metabolome_modules$MEs),
      ncol(rv$results$microbiome_modules$MEs),
      nrow(rv$results$iv_table),
      nrow(rv$results$passed_pathways),
      paste(
        apply(head(rv$results$passed_pathways, 10), 1, function(row) {
          sprintf("<tr><td>%s</td><td>%s</td><td>%.3f</td><td>%.3f</td><td>%.4f</td><td>%.2f</td></tr>",
                 row['module'], row['outcome'], as.numeric(row['beta']),
                 as.numeric(row['se']), as.numeric(row['p_value']), as.numeric(row['weak_F']))
        }),
        collapse = "\n"
      )
      )
      
      writeLines(html_content, file)
    }
  )
}

################################################################################
# RUN APP
################################################################################

shinyApp(ui = ui, server = server)
