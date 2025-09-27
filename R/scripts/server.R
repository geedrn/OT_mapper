library(shiny)
library(httr)
library(readr)
library(dplyr)
library(stringr)
library(DT)
library(GenomicRanges)

source("gggenome_functions.R")

server <- function(input, output, session) {
  results <- reactiveVal(NULL)
  
  observeEvent(input$run, {
    tryCatch({
      req(input$spacer, input$seed_length, input$pam)
      
      # スペースを除去し、バリデーション
      spacer <- clean_and_validate_input(input$spacer)
      pam <- clean_and_validate_input(input$pam)

      validate(need(nchar(spacer) == 20, "Error: Spacer must be 20 nt!"))
      
      spacer_and_pam <- paste0(spacer, input$pam)
      seed <- substr(spacer_and_pam, ifelse(input$seed_length == 12, 9, 13), nchar(spacer_and_pam))
      
      # Initialize progress
      progress <- shiny::Progress$new()
      progress$set(message = "Processing", value = 0)
      on.exit(progress$close())
      
      updateProgress <- function(detail) {
        progress$inc(amount = 1/6, detail = detail)
      }
      
      spacer_and_pam <- paste0(input$spacer, input$pam)
      print(paste("Spacer and PAM:", spacer_and_pam))
      
      seed <- substr(spacer_and_pam, ifelse(input$seed_length == 12, 9, 13), 23)
      print(paste("Seed:", seed))
      
      updateProgress("Downloading plus strand data")
      
      ot1_plus <- tryCatch({
        process_gggenome(query_gggenome(spacer_and_pam, "+", 3), spacer_and_pam)
      }, error = function(e) {
        stop(paste("Error querying GGGenome for plus strand full sequence:", e$message))
      })
      
      ot2_plus <- tryCatch({
        process_gggenome(query_gggenome(seed, "+", 1), seed)
      }, error = function(e) {
        stop(paste("Error querying GGGenome for plus strand seed sequence:", e$message))
      })
      
      updateProgress("Processing plus strand data")
      
      print("Plus strand data:")
      print(paste("ot1_plus rows:", nrow(ot1_plus)))
      print(paste("ot2_plus rows:", nrow(ot2_plus)))
      
      # Join plus strand results with sequence containment filtering using 'name', 'start', and 'end'
      plus_strand <- tryCatch({
        
        # Ensure that columns 'name', 'start', and 'end' exist in both data frames
        if (!all(c("name", "start", "end") %in% colnames(ot1_plus)) ||
            !all(c("name", "start", "end") %in% colnames(ot2_plus))) {
          stop("Required columns ('name', 'start', 'end') not found in one or both dataframes.")
        }
        
        # Create GRanges objects from ot1_plus and ot2_plus
        gr1 <- GRanges(
          seqnames = ot1_plus$name,
          ranges = IRanges(start = ot1_plus$start, end = ot1_plus$end),
          strand = ot1_plus$strand
        )
        
        gr2 <- GRanges(
          seqnames = ot2_plus$name,
          ranges = IRanges(start = ot2_plus$start, end = ot2_plus$end),
          strand = ot2_plus$strand
        )
        
        # Find overlaps between gr1 and gr2
        overlaps <- findOverlaps(gr1, gr2)
        
        # Extract indices of overlapping ranges
        ot1_indices <- queryHits(overlaps)
        ot2_indices <- subjectHits(overlaps)
        
        # Get overlapping entries from ot1_plus and ot2_plus
        ot1_plus_filtered <- ot1_plus[ot1_indices, ]
        ot2_plus_filtered <- ot2_plus[ot2_indices, ]
        
        # Perform the join on 'name' and apply further filtering if needed
        overlapping_results <- ot1_plus_filtered %>%
          inner_join(ot2_plus_filtered, by = "name") %>%
          filter(start.x <= start.y, end.y <= end.x, strand.x == strand.y)
        
        overlapping_results
        
      }, error = function(e) {
        stop(paste("Error joining filtered plus strand results:", e$message))
      })
      
      print(paste("Filtered plus_strand rows after join:", nrow(plus_strand)))
      
      updateProgress("Downloading minus strand data")
      
      ot1_minus <- tryCatch({
        process_gggenome(query_gggenome(spacer_and_pam, "-", 3), spacer_and_pam)
      }, error = function(e) {
        stop(paste("Error querying GGGenome for minus strand full sequence:", e$message))
      })
      
      ot2_minus <- tryCatch({
        process_gggenome(query_gggenome(seed, "-", 0), seed)
      }, error = function(e) {
        stop(paste("Error querying GGGenome for minus strand seed sequence:", e$message))
      })
      
      updateProgress("Processing minus strand data")
      
      print("Minus strand data:")
      print(paste("ot1_minus rows:", nrow(ot1_minus)))
      print(paste("ot2_minus rows:", nrow(ot2_minus)))
      
      # Join minus strand results with sequence containment filtering using 'name', 'start', and 'end'
      minus_strand <- tryCatch({
        
        # Ensure that columns 'name', 'start', and 'end' exist in both data frames
        if (!all(c("name", "start", "end") %in% colnames(ot1_minus)) ||
            !all(c("name", "start", "end") %in% colnames(ot2_minus))) {
          stop("Required columns ('name', 'start', 'end') not found in one or both dataframes.")
        }
        
        # Load the GenomicRanges package
        if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
          install.packages("BiocManager")
          BiocManager::install("GenomicRanges")
        }
        library(GenomicRanges)
        
        # Create GRanges objects from ot1_minus and ot2_minus
        gr1 <- GRanges(
          seqnames = ot1_minus$name,
          ranges = IRanges(start = ot1_minus$start, end = ot1_minus$end),
          strand = ot1_minus$strand
        )
        
        gr2 <- GRanges(
          seqnames = ot2_minus$name,
          ranges = IRanges(start = ot2_minus$start, end = ot2_minus$end),
          strand = ot2_minus$strand
        )
        
        # Find overlaps between gr1 and gr2
        overlaps <- findOverlaps(gr1, gr2)
        
        # Extract indices of overlapping ranges
        ot1_indices <- queryHits(overlaps)
        ot2_indices <- subjectHits(overlaps)
        
        # Get overlapping entries from ot1_minus and ot2_minus
        ot1_minus_filtered <- ot1_minus[ot1_indices, ]
        ot2_minus_filtered <- ot2_minus[ot2_indices, ]
        
        # Perform the join on 'name' and apply further filtering if needed
        overlapping_results <- ot1_minus_filtered %>%
          inner_join(ot2_minus_filtered, by = "name") %>%
          filter(start.x <= start.y, end.y <= end.x, strand.x == strand.y)
        
        overlapping_results
        
      }, error = function(e) {
        stop(paste("Error joining filtered minus strand results:", e$message))
      })
      
      # Combine plus_strand and minus_strand results
      ot_list <- tryCatch({
        combined_results <- bind_rows(plus_strand, minus_strand)
        combined_results
      }, error = function(e) {
        stop(paste("Error combining results:", e$message))
      })
      
      # After processing the data and before setting the results, add this code:
      
      print("Final Off-Target Sites:")
      if(nrow(ot_list) > 0) {
        summary_table <- ot_list %>%
          mutate(
            Chr = name,
            Start = start.x,
            End = end.x,
            Strand = strand.x,
            Sequence = sbjct.x,
            Mismatches = edit.x
          ) %>%
          select(Chr, Start, End, Strand, Sequence, Mismatches) %>%
          arrange(Mismatches, Chr, Start)
        
        print(summary_table)
        
        print(paste("Total off-target sites found:", nrow(summary_table)))
        print("Mismatches distribution:")
        print(table(summary_table$Mismatches))
      } else {
        print("No off-target sites found.")
      }
      
      # Then, modify the results() call to include this summary:
      
      results(list(
        ot_list_final = ot_list,
        summary_table = summary_table,  # Add this line
        ucsc_list = if(nrow(ot_list) > 0) {
          ot_list %>%
            mutate(
              ucsc = paste0(sbjct.x, ",", name, ":", start.x, "-", end.x)
            ) %>%
            select(ucsc)
        } else {
          data.frame(ucsc = character())
        },
        ot_candidate_bed = if(nrow(ot_list) > 0) {
          ot_list %>%
            select(name, start.x, end.x)
        } else {
          data.frame(name = character(), start = integer(), end = integer())
        }
      ))
      
      # Finally, update the renderDT() function to use the summary table:
      
      output$ot_table <- renderDT({
        req(results())
        datatable(results()$summary_table, options = list(pageLength = 10))
      })
      
      updateProgress("Complete")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      print(paste("Error:", e$message))
    })
  })
  
  output$ot_table <- renderDT({
    req(results())
    datatable(results()$ot_list_final)
  })
  
  output$download_ot <- downloadHandler(
    filename = function() { "OT_list_final.csv" },
    content = function(file) {
      write_csv(results()$ot_list_final, file)
    }
  )
  
  output$download_ucsc <- downloadHandler(
    filename = function() { "UCSC_list_final.csv" },
    content = function(file) {
      write_csv(results()$ucsc_list, file)
    }
  )
  
  output$download_bed <- downloadHandler(
    filename = function() { "OT_candidate.bed" },
    content = function(file) {
      write_tsv(results()$ot_candidate_bed, file, col_names = FALSE)
    }
  )
}