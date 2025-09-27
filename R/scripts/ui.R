library(shiny)
library(shinydashboard)
library(DT)

# UI
ui <- dashboardPage(
  dashboardHeader(title = "OT Detector"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Input", tabName = "input", icon = icon("edit")),
      menuItem("Results", tabName = "results", icon = icon("table"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              fluidRow(
                box(
                  title = "Input Parameters",
                  textInput("spacer", "gRNA Sequence (20 nt, without PAM):", ""),
                  selectInput("seed_length", "Seed Sequence Length:", choices = c(8, 12)),
                  selectInput("pam", "PAM:", choices = c("NGG", "NRG")),
                  actionButton("run", "Run OT Detection"),
                  width = 12
                )
              )
      ),
      tabItem(tabName = "results",
              tabsetPanel(
                tabPanel("Spacer Match", DTOutput("spacer_match_table")),
                tabPanel("Spacer Mismatch", DTOutput("spacer_mismatch_table")),
                tabPanel("Seed Match", DTOutput("seed_match_table")),
                tabPanel("Seed Mismatch", DTOutput("seed_mismatch_table")),
                tabPanel("All Results", DTOutput("all_results_table"))
              ),
              downloadButton("download_all", "Download All Results")
      )
    )
  )
)