library(shiny)
library(shinydashboard)
library(DT)

ui <- dashboardPage(
  dashboardHeader(title = "CRISPR Off-target Detector"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "input", icon = icon("edit")),
      menuItem("Results", tabName = "results", icon = icon("table")),
      menuItem("About", tabName = "about", icon = icon("info-circle"))
    )
  ),
  dashboardBody(
    tabItems(
      # Input tab
      tabItem(tabName = "input",
        fluidRow(
          box(
            title = "Input Parameters", status = "primary", solidHeader = TRUE,
            width = 12,
            textInput("spacer", 
                      label = "gRNA Spacer Sequence (20 nt, without PAM):",
                      placeholder = "e.g., TCGCCCAGCGACCCTGCTCC",
                      value = ""),
            helpText("Enter a 20-nucleotide spacer sequence (A, C, G, T only)"),
            
            sliderInput("seed_length",
                       label = "Seed Sequence Length:",
                       min = 8, max = 12, value = 12, step = 1),
            helpText("Seed length: 8 to 12 nucleotides (inclusive)"),
            
            selectInput("pam", 
                       label = "PAM Sequence:",
                       choices = list("NGG" = "NGG", "NRG" = "NRG"),
                       selected = "NGG"),
            helpText("Protospacer Adjacent Motif sequence"),
            
            br(),
            actionButton("run", "Run Analysis", 
                        class = "btn-primary btn-lg",
                        icon = icon("play")),
            br(), br(),
            helpText("Note: Analysis may take several minutes depending on your internet connection and the number of candidates found.")
          )
        )
      ),
      
      # Results tab
      tabItem(tabName = "results",
        fluidRow(
          box(
            title = "Summary Statistics", status = "info", solidHeader = TRUE,
            width = 12,
            verbatimTextOutput("summary_stats")
          )
        ),
        fluidRow(
          box(
            title = "Off-target Candidates", status = "success", solidHeader = TRUE,
            width = 12,
            DTOutput("results_table"),
            br(),
            downloadButton("download_summary", "Download Summary (CSV)", class = "btn-primary"),
            downloadButton("download_full", "Download Full Results (CSV)", class = "btn-primary"),
            downloadButton("download_annotated", "Download Annotated (TSV)", class = "btn-primary")
          )
        )
      ),
      
      # About tab
      tabItem(tabName = "about",
        fluidRow(
          box(
            title = "About CRISPR Off-target Detector", status = "info", solidHeader = TRUE,
            width = 12,
            h3("Software Description"),
            p("This Shiny application provides a web interface for CRISPR-Cas9 off-target detection and annotation."),
            h4("Features:"),
            tags$ul(
              tags$li("GGGenome API integration for off-target candidate search"),
              tags$li("PAM exact match filtering"),
              tags$li("Coordinate-based overlap detection using bedtools"),
              tags$li("Gene annotation (exon/intron mapping)"),
              tags$li("Interactive results table with download options")
            ),
            h4("Usage:"),
            tags$ol(
              tags$li("Enter a 20-nucleotide gRNA spacer sequence (without PAM)"),
              tags$li("Select seed length (8-12 nucleotides)"),
              tags$li("Choose PAM sequence (default: NGG)"),
              tags$li("Click 'Run Analysis' and wait for results"),
              tags$li("View results in the Results tab and download as needed")
            ),
            h4("Requirements:"),
            tags$ul(
              tags$li("Internet connection (for GGGenome API)"),
              tags$li("bedtools (for overlap detection and annotation)"),
              tags$li("R packages: shiny, shinydashboard, DT, httr, readr, dplyr, stringr, GenomicRanges")
            ),
            br(),
            p("For command-line usage, see the main README.md file.")
          )
        )
      )
    )
  )
)
