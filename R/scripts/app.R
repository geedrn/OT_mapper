# =============================================================================
# Shiny App Entry Point
# =============================================================================
# This file allows running the Shiny app with: shiny::runApp()
# =============================================================================

library(shiny)

# Source UI and server
source("ui.R")
source("server.R")

# Create Shiny app object
shinyApp(ui = ui, server = server)

