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
  analysis_completed <- reactiveVal(FALSE)
  
  # Observe when results are updated and switch to Results tab
  observe({
    if (!is.null(results()) && analysis_completed()) {
      updateTabItems(session, "sidebar", "results")
      analysis_completed(FALSE)  # Reset flag
    }
  })
  
  observeEvent(input$run, {
    analysis_completed(FALSE)  # Reset flag at start
    tryCatch({
      req(input$spacer, input$seed_length, input$pam, input$full_mismatch, input$seed_mismatch)
      
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
        suppressMessages({
          gggenome_to_dataframe(spacer, input$seed_length, input$pam, 
                               full_mismatch = input$full_mismatch, 
                               seed_mismatch = input$seed_mismatch)
        })
      }, error = function(e) {
        stop(paste("Error querying GGGenome API:", e$message), call. = FALSE)
      })
      
      if (nrow(list$plus_full) == 0 && nrow(list$minus_full) == 0) {
        stop("No candidates found in GGGenome search. Try adjusting parameters.", call. = FALSE)
      }
      
      # Step 2: PAM filtering
      updateProgress("Step 2/6: Filtering PAM matches...")
      filtered_list <- tryCatch({
        suppressMessages({
          filter_exact_pam(list, input$pam)
        })
      }, error = function(e) {
        stop(paste("Error filtering PAM:", e$message), call. = FALSE)
      })
      
      # Step 3: Overlap detection
      updateProgress("Step 3/6: Detecting overlaps...")
      overlaps <- tryCatch({
        suppressMessages({
          find_overlaps(filtered_list)
        })
      }, error = function(e) {
        stop(paste("Error detecting overlaps:", e$message), call. = FALSE)
      })
      
      # Step 4: Combine results
      updateProgress("Step 4/6: Combining results...")
      combined_df <- tryCatch({
        suppressMessages({
          combine_results(overlaps)
        })
      }, error = function(e) {
        stop(paste("Error combining results:", e$message), call. = FALSE)
      })
      
      if (nrow(combined_df) == 0) {
        stop("No off-target candidates found after filtering.", call. = FALSE)
      }
      
      # Step 5: Gene annotation (optional, if annotation DBs are available)
      updateProgress("Step 5/6: Annotating with gene information...")
      annotated_df <- NULL
      # Annotation DBs are in R/scripts/data/
      exon_db <- file.path(script_dir, "data/UCSC_exons_modif_canonical.bed")
      intron_db <- file.path(script_dir, "data/UCSC_introns_modif_canonical.bed")
      
      if (file.exists(exon_db) && file.exists(intron_db)) {
        # Check if bedtools is available
        bedtools_available <- system("which bedtools", ignore.stdout = TRUE, ignore.stderr = TRUE) == 0
        
        if (!bedtools_available) {
          showNotification("bedtools not found. Skipping gene annotation. Install bedtools to enable annotation.", 
                          type = "warning", duration = 5)
          print("Warning: bedtools not found. Gene annotation skipped.")
        } else {
          annotated_df <- tryCatch({
            # Suppress print() output in Shiny app (output_file=NULL means no file output)
            suppressMessages({
              annotate_with_bedtools(combined_df, exon_db, intron_db, output_file = NULL)
            })
          }, error = function(e) {
            showNotification(paste("Annotation failed:", e$message), 
                            type = "warning", duration = 5)
            print(paste("Annotation error:", e$message))
            NULL
          })
        }
      } else {
        showNotification("Annotation database files not found. Skipping gene annotation.", 
                        type = "warning", duration = 5)
        print("Warning: Annotation database files not found.")
      }
      
      # Step 6: Prepare summary table and generate Primer-BLAST links
      updateProgress("Step 6/6: Preparing results...")
      
      # Create summary table first
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
      
      # Generate Primer-BLAST links if annotated data is available
      primer_blast_df <- NULL
      if (!is.null(annotated_df) && nrow(annotated_df) > 0) {
        tryCatch({
          # annotated_df format from bedtools intersect -wa -wb:
          # V1-V6: input BED (chrom, start, end, name, score, strand)
          # V7-V12: annotation BED (chrom, start, end, geneSymbol, feature_type, feature_number)
          # We need to use input BED coordinates (V1, V2, V3) for Primer-BLAST
          
          # Extract unique coordinates from annotated_df (input BED columns)
          unique_coords <- annotated_df[, 1:3]
          colnames(unique_coords) <- c("chrom", "start", "end")
          unique_coords <- unique_coords[!duplicated(unique_coords), ]
          
          temp_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", sample(1000:9999, 1))
          temp_primer_input <- file.path(tempdir(), paste0("primer_input_", temp_id, ".tsv"))
          temp_primer_output <- file.path(tempdir(), paste0("primer_output_", temp_id, ".tsv"))
          
          # Write unique coordinates to temporary file
          write.table(unique_coords, temp_primer_input, 
                     row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
          
          # Generate Primer-BLAST links
          generate_primer_blast_links(temp_primer_input, temp_primer_output)
          
          # Read back the results with Primer-BLAST URLs
          if (file.exists(temp_primer_output) && file.size(temp_primer_output) > 0) {
            primer_blast_df <- read.table(temp_primer_output, header = FALSE, sep = "\t", 
                                          stringsAsFactors = FALSE, quote = "")
            
            # Check column count - generate_primer_blast_links outputs:
            # Input: 3 columns (chrom, start, end)
            # The function adds V6 column with URL, so output has 6 columns total
            # But when reading, if V4 and V5 are empty, they may be missing
            print(paste("Primer-BLAST output columns:", ncol(primer_blast_df)))
            print(paste("Primer-BLAST output rows:", nrow(primer_blast_df)))
            
            # Handle different column counts
            if (ncol(primer_blast_df) == 4) {
              # 4 columns: V1=RefSeq_ID, V2=Start, V3=End, V4=Primer_BLAST_URL (V6 was written as V4)
              colnames(primer_blast_df) <- c("RefSeq_ID", "Start", "End", "Primer_BLAST_URL")
            } else if (ncol(primer_blast_df) == 6) {
              # 6 columns: V1=RefSeq_ID, V2=Start, V3=End, V4=NA, V5=NA, V6=Primer_BLAST_URL
              colnames(primer_blast_df) <- c("RefSeq_ID", "Start", "End", "V4", "V5", "Primer_BLAST_URL")
              primer_blast_df <- primer_blast_df[, c("RefSeq_ID", "Start", "End", "Primer_BLAST_URL")]
            } else {
              # Try to extract URL from the last column
              print(paste("Unexpected column count:", ncol(primer_blast_df), "- trying to extract URL from last column"))
              colnames(primer_blast_df)[1:3] <- c("RefSeq_ID", "Start", "End")
              primer_blast_df$Primer_BLAST_URL <- primer_blast_df[, ncol(primer_blast_df)]
              primer_blast_df <- primer_blast_df[, c("RefSeq_ID", "Start", "End", "Primer_BLAST_URL")]
            }
            
            if (ncol(primer_blast_df) == 4 && "Primer_BLAST_URL" %in% colnames(primer_blast_df)) {
              # Store Primer-BLAST URLs in primer_blast_df for later use in download files
              # Don't add to summary_table (not displayed in table)
              matched_count <- nrow(primer_blast_df)
              print(paste("Primer-BLAST links generated:", matched_count, "links"))
            } else {
              print(paste("Error: Primer-BLAST output has unexpected column count:", ncol(primer_blast_df)))
            }
            
            # Clean up temporary files
            unlink(c(temp_primer_input, temp_primer_output))
          } else {
            print("Primer-BLAST output file is empty or does not exist")
          }
        }, error = function(e) {
          print(paste("Primer-BLAST link generation failed:", e$message))
          print(paste("Error details:", toString(e)))
          # Continue without Primer-BLAST links
        })
      } else {
        print("No annotated data available for Primer-BLAST link generation")
      }
      
      # Store results
      results(list(
        summary_table = summary_table,
        combined_df = combined_df,
        annotated_df = annotated_df,
        primer_blast_df = primer_blast_df,
        spacer = spacer,
        seed_length = input$seed_length,
        pam = input$pam
      ))
      
      showNotification("Analysis completed successfully!", type = "message")
      
      # Set flag to trigger tab switch
      analysis_completed(TRUE)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
      print(paste("Error:", e$message))
    })
  })
  
  # Render results table
  output$results_table <- renderDT({
    req(results())
    df <- results()$summary_table
    
    # Remove Primer_BLAST_URL column if it exists (not displayed in table, only in download files)
    if ("Primer_BLAST_URL" %in% colnames(df)) {
      df <- df %>% select(-Primer_BLAST_URL)
    }
    
    datatable(
      df,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        order = list(list(5, 'asc'))  # Sort by mismatches if available
      ),
      rownames = FALSE
    )
  })
  
  # Download handler - All results in ZIP
  output$download_all <- downloadHandler(
    filename = function() { 
      paste0("OT_results_", Sys.Date(), ".zip") 
    },
    content = function(file) {
      # Create temporary directory for files
      temp_dir <- tempdir()
      temp_files <- c()
      
      # Generate filenames
      date_str <- format(Sys.Date(), "%Y%m%d")
      summary_file <- file.path(temp_dir, paste0("OT_summary_", date_str, ".csv"))
      full_file <- file.path(temp_dir, paste0("OT_full_", date_str, ".csv"))
      annotated_file <- file.path(temp_dir, paste0("OT_annotated_", date_str, ".tsv"))
      
      # Write summary CSV (without Primer-BLAST URLs)
      write.csv(results()$summary_table, summary_file, row.names = FALSE)
      temp_files <- c(temp_files, summary_file)
      
      # Write full results CSV (without Primer-BLAST URLs)
      write.csv(results()$combined_df, full_file, row.names = FALSE)
      temp_files <- c(temp_files, full_file)
      
      # Write annotated TSV with Primer-BLAST links (if available)
      if (!is.null(results()$annotated_df) && nrow(results()$annotated_df) > 0) {
        # If primer_blast_df exists, merge with annotated_df
        if (!is.null(results()$primer_blast_df) && nrow(results()$primer_blast_df) > 0) {
          # Merge annotated_df with primer_blast_df by coordinates
          # annotated_df format: V1=chrom, V2=start, V3=end, ...
          annotated_with_primer <- results()$annotated_df
          annotated_with_primer$Primer_BLAST_URL <- NA
          for (i in 1:nrow(annotated_with_primer)) {
            matching_idx <- which(
              results()$primer_blast_df$Start == annotated_with_primer[i, 2] &
              results()$primer_blast_df$End == annotated_with_primer[i, 3]
            )
            if (length(matching_idx) > 0) {
              annotated_with_primer$Primer_BLAST_URL[i] <- results()$primer_blast_df$Primer_BLAST_URL[matching_idx[1]]
            }
          }
          write.table(annotated_with_primer, annotated_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
        } else {
          write.table(results()$annotated_df, annotated_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
        }
        temp_files <- c(temp_files, annotated_file)
      } else {
        # Create empty file if no annotation
        write.table(data.frame(), annotated_file, row.names = FALSE, col.names = FALSE, sep = "\t")
        temp_files <- c(temp_files, annotated_file)
      }
      
      # Create ZIP file
      zip_file <- file.path(temp_dir, paste0("OT_results_", date_str, ".zip"))
      old_wd <- getwd()
      setwd(temp_dir)
      tryCatch({
        utils::zip(zip_file, files = basename(temp_files))
        file.copy(zip_file, file)
      }, finally = {
        setwd(old_wd)
        # Clean up temporary files
        unlink(temp_files)
        unlink(zip_file)
      })
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
