library(shiny)
library(shinydashboard)
library(DT)
library(httr)
library(readr)
library(dplyr)
library(stringr)
library(GenomicRanges)

# Source required scripts
# Get script directory (works when sourced or run from RStudio)
script_dir <- tryCatch({
  # Try to get the directory of the current script
  if (exists("RStudio.Version")) {
    # Running in RStudio
    rstudioapi::getActiveDocumentContext()$path
  } else {
    # Try commandArgs method
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("--file=", cmd_args, value = TRUE)
    if (length(file_arg) > 0) {
      dirname(normalizePath(sub("--file=", "", file_arg)))
    } else {
      # Fallback: assume scripts are in the same directory
      getwd()
    }
  }
}, error = function(e) {
  # Fallback to current working directory
  getwd()
})

# If script_dir is a file path, get its directory
if (file.exists(script_dir) && !dir.exists(script_dir)) {
  script_dir <- dirname(script_dir)
}

# Ensure script_dir points to the R/scripts directory
if (!file.exists(file.path(script_dir, "OT_mapper.R"))) {
  # Try common locations relative to shiny directory
  if (file.exists("../R/scripts/OT_mapper.R")) {
    script_dir <- "../R/scripts"
  } else if (file.exists("R/scripts/OT_mapper.R")) {
    script_dir <- "R/scripts"
  } else if (file.exists("scripts/OT_mapper.R")) {
    script_dir <- "scripts"
  } else {
    stop("Cannot find OT_mapper.R. Please ensure R/scripts directory is accessible.")
  }
}

# Source core functions
source(file.path(script_dir, "OT_mapper.R"), local = TRUE)
source(file.path(script_dir, "OT_annotator.R"), local = TRUE)
source(file.path(script_dir, "primer_generate.R"), local = TRUE)

server <- function(input, output, session) {
  results <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    tryCatch({
      req(input$spacer, input$seed_length, input$pam)
      
      # Clean and validate spacer
      spacer <- gsub("[[:space:]]", "", toupper(input$spacer))
      
      validate(
        need(nchar(spacer) == 20, "Error: Spacer must be exactly 20 nucleotides!"),
        need(grepl("^[ACGT]+$", spacer), "Error: Spacer must contain only A, C, G, T characters!")
      )
      
      # Initialize progress
      progress <- shiny::Progress$new()
      progress$set(message = "Processing", value = 0)
      on.exit(progress$close())
      
      updateProgress <- function(detail) {
        progress$inc(amount = 1/6, detail = detail)
      }
      
      # Step 1: GGGenome search
      updateProgress("Step 1/6: Searching GGGenome API...")
      list <- tryCatch({
        gggenome_to_dataframe(spacer, input$seed_length, input$pam)
      }, error = function(e) {
        stop(paste("Error querying GGGenome API:", e$message))
      })
      
      if (nrow(list$plus_full) == 0 && nrow(list$minus_full) == 0) {
        stop("No candidates found in GGGenome search. Try adjusting parameters.")
      }
      
      # Step 2: PAM filtering
      updateProgress("Step 2/6: Filtering PAM matches...")
      filtered_list <- tryCatch({
        filter_exact_pam(list, input$pam)
      }, error = function(e) {
        stop(paste("Error filtering PAM:", e$message))
      })
      
      # Step 3: Overlap detection
      updateProgress("Step 3/6: Detecting overlaps...")
      overlaps <- tryCatch({
        find_overlaps(filtered_list)
      }, error = function(e) {
        stop(paste("Error detecting overlaps:", e$message))
      })
      
      # Step 4: Combine results
      updateProgress("Step 4/6: Combining results...")
      combined_df <- tryCatch({
        combine_results(overlaps)
      }, error = function(e) {
        stop(paste("Error combining results:", e$message))
      })
      
      if (nrow(combined_df) == 0) {
        stop("No off-target candidates found after filtering.")
      }
      
      # Step 5: Gene annotation (optional, if annotation DBs are available)
      updateProgress("Step 5/6: Annotating with gene information...")
      annotated_df <- NULL
      # Annotation DBs are in R/scripts/data/
      exon_db <- file.path(script_dir, "data/UCSC_exons_modif_canonical.bed")
      intron_db <- file.path(script_dir, "data/UCSC_introns_modif_canonical.bed")
      
      if (file.exists(exon_db) && file.exists(intron_db)) {
        annotated_df <- tryCatch({
          annotate_with_bedtools(combined_df, exon_db, intron_db, output_file = NULL)
        }, error = function(e) {
          warning(paste("Annotation failed:", e$message))
          NULL
        })
      }
      
      # Step 6: Prepare summary table
      updateProgress("Step 6/6: Preparing results...")
      
      # Create summary table
      summary_table <- combined_df %>%
        mutate(
          Chromosome = chrom,
          Start = start,
          End = end,
          Strand = ifelse("strand" %in% colnames(combined_df), strand, "+"),
          Sequence = ifelse("sbjct" %in% colnames(combined_df), sbjct, ""),
          Mismatches = ifelse("mis" %in% colnames(combined_df), mis, NA)
        ) %>%
        select(Chromosome, Start, End, Strand, Sequence, Mismatches) %>%
        arrange(if("Mismatches" %in% colnames(.)) Mismatches else Chromosome, Chromosome, Start)
      
      # Store results
      results(list(
        summary_table = summary_table,
        combined_df = combined_df,
        annotated_df = annotated_df,
        spacer = spacer,
        seed_length = input$seed_length,
        pam = input$pam
      ))
      
      showNotification("Analysis completed successfully!", type = "success")
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
      print(paste("Error:", e$message))
    })
  })
  
  # Render results table
  output$results_table <- renderDT({
    req(results())
    datatable(
      results()$summary_table,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(5, 'asc'))  # Sort by mismatches if available
      ),
      rownames = FALSE
    )
  })
  
  # Download handlers
  output$download_summary <- downloadHandler(
    filename = function() { 
      paste0("OT_summary_", Sys.Date(), ".csv") 
    },
    content = function(file) {
      write.csv(results()$summary_table, file, row.names = FALSE)
    }
  )
  
  output$download_full <- downloadHandler(
    filename = function() { 
      paste0("OT_full_", Sys.Date(), ".csv") 
    },
    content = function(file) {
      write.csv(results()$combined_df, file, row.names = FALSE)
    }
  )
  
  output$download_annotated <- downloadHandler(
    filename = function() { 
      paste0("OT_annotated_", Sys.Date(), ".tsv") 
    },
    content = function(file) {
      if (!is.null(results()$annotated_df) && nrow(results()$annotated_df) > 0) {
        write.table(results()$annotated_df, file, row.names = FALSE, col.names = FALSE, sep = "\t")
      } else {
        write.table(data.frame(), file, row.names = FALSE, col.names = FALSE, sep = "\t")
      }
    }
  )
  
  # Summary statistics
  output$summary_stats <- renderText({
    req(results())
    r <- results()
    total <- nrow(r$summary_table)
    annotated <- if (!is.null(r$annotated_df)) nrow(r$annotated_df) else 0
    
    paste0(
      "Total candidates: ", total, "\n",
      "Annotated: ", annotated, "\n",
      "Spacer: ", r$spacer, "\n",
      "Seed length: ", r$seed_length, "\n",
      "PAM: ", r$pam
    )
  })
}
