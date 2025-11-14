#!/usr/bin/env Rscript

# =============================================================================
# CRISPR Off-target Analysis Script
# =============================================================================
# Command-line script for complete off-target analysis pipeline
#
# Usage:
#   Rscript analysis.R [options]
#   Rscript analysis.R -s <spacer> -l <seed_length> [options]
#
# Options:
#   -s, --spacer      20nt gRNA spacer sequence (without PAM) [required]
#   -l, --seed-length Seed sequence length: 8 to 12 [default: 12]
#   -p, --pam         PAM sequence [default: NGG]
#   -g, --genome      Genome assembly [default: hg38]
#   -f, --full-mismatch Mismatch tolerance for full sequence [default: 3]
#   -m, --seed-mismatch Mismatch tolerance for seed sequence [default: 1]
#   -e, --exon-db     Exon annotation BED file [default: data/UCSC_exons_modif_canonical.bed]
#   -i, --intron-db   Intron annotation BED file [default: data/UCSC_introns_modif_canonical.bed]
#   -o, --output      Output file path [default: annotated_offtargets.tsv]
#   -h, --help        Show this help message
# =============================================================================

suppressPackageStartupMessages({
  options(warn = 1)  # Show warnings immediately
})

# =============================================================================
# Load Required Packages
# =============================================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package 'optparse' is required. Install with: install.packages('optparse')")
}

library(optparse)

# =============================================================================
# Logging Functions
# =============================================================================

log_info <- function(...) {
  cat("[INFO]", ..., "\n")
}

log_success <- function(...) {
  cat("[SUCCESS]", ..., "\n")
}

log_warning <- function(...) {
  cat("[WARNING]", ..., "\n")
}

log_error <- function(...) {
  cat("[ERROR]", ..., "\n", file = stderr())
}

log_step <- function(...) {
  cat("\n=== ", ..., " ===\n", sep = "")
}

# =============================================================================
# Command-Line Argument Parsing
# =============================================================================

option_list <- list(
  make_option(c("-s", "--spacer"), 
              type = "character", 
              default = NULL,
              help = "20nt gRNA spacer sequence (without PAM) [required if using options]",
              metavar = "SEQUENCE"),
  
  make_option(c("-l", "--seed-length"), 
              type = "integer", 
              default = 12,
              help = "Seed sequence length: 8 to 12 [default: %default]",
              metavar = "LENGTH"),
  
  make_option(c("-p", "--pam"), 
              type = "character", 
              default = "NGG",
              help = "PAM sequence [default: %default]",
              metavar = "PAM"),
  
  make_option(c("-g", "--genome"), 
              type = "character", 
              default = "hg38",
              help = "Genome assembly [default: %default]",
              metavar = "GENOME"),
  
  make_option(c("-f", "--full-mismatch"), 
              type = "integer", 
              default = 3,
              help = "Mismatch tolerance for full sequence [default: %default]",
              metavar = "NUMBER"),
  
  make_option(c("-m", "--seed-mismatch"), 
              type = "integer", 
              default = 1,
              help = "Mismatch tolerance for seed sequence [default: %default]",
              metavar = "NUMBER"),
  
  make_option(c("-e", "--exon-db"), 
              type = "character", 
              default = NULL,
              help = "Exon annotation BED file [default: data/UCSC_exons_modif_canonical.bed]",
              metavar = "FILE"),
  
  make_option(c("-i", "--intron-db"), 
              type = "character", 
              default = NULL,
              help = "Intron annotation BED file [default: data/UCSC_introns_modif_canonical.bed]",
              metavar = "FILE"),
  
  make_option(c("-o", "--output"), 
              type = "character", 
              default = "annotated_offtargets.tsv",
              help = "Output file path [default: %default]",
              metavar = "FILE"),
  
  make_option(c("--skip-primer"), 
              action = "store_true", 
              default = FALSE,
              help = "Skip Primer-BLAST link generation [default: %default]"),
  
  make_option(c("-h", "--help"), 
              action = "store_true", 
              default = FALSE,
              help = "Show this help message")
)

opt_parser <- OptionParser(
  option_list = option_list,
  usage = "usage: %prog [options] OR %prog <spacer> [seed_length] [pam] [output_file]",
  description = "CRISPR Off-target Analysis - Complete pipeline for detecting and annotating off-target candidates\n\nBy default, this script runs the complete pipeline:\n  1. Off-target candidate detection\n  2. Gene annotation mapping (automatic)\n  3. Primer-BLAST link generation (automatic)",
  epilogue = "Examples:\n  Rscript analysis.R -s TCGCCCAGCGACCCTGCTCC -l 12\n  Rscript analysis.R -s TCGCCCAGCGACCCTGCTCC -l 8 -p NGG -g hg38 -o results.tsv\n  Rscript analysis.R -s TCGCCCAGCGACCCTGCTCC -l 12 --skip-primer  # Skip primer generation\n  Rscript analysis.R TCGCCCAGCGACCCTGCTCC 12 NGG results.tsv  # Legacy positional args"
)

opt <- parse_args(opt_parser, positional_arguments = TRUE)
positional_args <- opt$args
opt <- opt$options

# Show help and exit
if (opt$help) {
  print_help(opt_parser)
  quit(status = 0)
}

# =============================================================================
# Error Handling Wrapper
# =============================================================================

try_catch_with_message <- function(expr, error_msg) {
  tryCatch(
    expr,
    error = function(e) {
      log_error(error_msg)
      log_error("Details:", conditionMessage(e))
      stop(e)
    }
  )
}

# =============================================================================
# Input Validation
# =============================================================================

validate_inputs <- function(spacer, seed_length, pam) {
  log_info("Validating inputs...")
  
  # Check spacer
  if (is.null(spacer) || nchar(spacer) == 0) {
    stop("Spacer sequence is required. Use -s/--spacer option or provide as first positional argument.")
  }
  
  if (nchar(spacer) != 20) {
    stop(paste0("Invalid spacer length: ", nchar(spacer), " nucleotides.\n",
                "Expected: 20 nucleotides (standard length for SpCas9 gRNAs)\n",
                "Received: ", nchar(spacer), " nucleotides\n",
                "Please provide a 20-nucleotide spacer sequence without the PAM."))
  }
  
  if (!grepl("^[ACGTacgt]+$", spacer)) {
    stop(paste0("Invalid characters in spacer sequence.\n",
                "Spacer must contain only A, C, G, T (case-insensitive).\n",
                "Received: ", spacer))
  }
  
  spacer <- toupper(spacer)
  
  # Check seed length
  if (is.na(seed_length) || seed_length < 8 || seed_length > 12) {
    stop(paste0("Invalid seed length: ", seed_length, "\n",
                "Expected: Integer between 8 and 12 (inclusive)\n",
                "Received: ", seed_length))
  }
  
  # Check PAM
  if (nchar(pam) == 0) {
    stop("PAM sequence cannot be empty.")
  }
  
  log_success("Input validation passed")
  return(spacer)
}

# =============================================================================
# Main Analysis Function
# =============================================================================

main <- function() {
  log_step("CRISPR Off-target Analysis")
  cat("Version: 1.0\n\n")
  
  # Handle backward compatibility: if spacer not provided via option, check positional args
  if (is.null(opt$spacer) && length(positional_args) > 0) {
    # Positional arguments detected (legacy mode)
    spacer <- positional_args[1]
    seed_length <- if (length(positional_args) >= 2) as.numeric(positional_args[2]) else opt$`seed-length`
    pam <- if (length(positional_args) >= 3) positional_args[3] else opt$pam
    output_file <- if (length(positional_args) >= 4) positional_args[4] else opt$output
    log_info("Using legacy positional arguments (consider using -s, -l, -p, -o options)")
  } else {
    # Use options (new preferred method)
    spacer <- opt$spacer
    seed_length <- opt$`seed-length`
    pam <- opt$pam
    output_file <- opt$output
  }
  
  # Set other options (not available in legacy mode)
  genome <- opt$genome
  full_mismatch <- opt$`full-mismatch`
  seed_mismatch <- opt$`seed-mismatch`
  
  # Validate inputs
  spacer <- validate_inputs(spacer, seed_length, pam)
  
  # Display configuration
  log_info("Configuration:")
  cat("  Spacer:", spacer, "\n")
  cat("  Seed length:", seed_length, "\n")
  cat("  PAM:", pam, "\n")
  cat("  Genome:", genome, "\n")
  cat("  Full sequence mismatch tolerance:", full_mismatch, "\n")
  cat("  Seed sequence mismatch tolerance:", seed_mismatch, "\n")
  cat("  Output file:", output_file, "\n\n")
  
  # Check for required R packages
  log_info("Checking required R packages...")
  required_packages <- c("utils", "httr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop(paste0("Missing required R packages: ", paste(missing_packages, collapse = ", "), "\n",
                "Please install them using: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))"))
  }
  log_success("All required packages are available")
  
  # Check for bedtools
  log_info("Checking for bedtools...")
  bedtools_check <- tryCatch(
    {
      system("which bedtools", ignore.stdout = TRUE, ignore.stderr = TRUE)
    },
    error = function(e) 1
  )
  
  if (bedtools_check != 0) {
    log_warning("bedtools not found in PATH.")
    log_warning("Some functions may not work correctly without bedtools.")
    log_info("Install bedtools: conda install -c bioconda bedtools")
  } else {
    log_success("bedtools is available")
  }
  
  # Source required scripts
  log_step("Loading Scripts")
  script_dir <- dirname(commandArgs(trailingOnly = FALSE)[grep("--file=", commandArgs(trailingOnly = FALSE))][1])
  script_dir <- sub("--file=", "", script_dir)
  
  if (script_dir == "") {
    script_dir <- getwd()
  }
  
  mapper_path <- file.path(script_dir, "OT_mapper.R")
  annotator_path <- file.path(script_dir, "OT_annotator.R")
  
  log_info("Script directory:", script_dir)
  
  if (!file.exists(mapper_path)) {
    stop(paste0("OT_mapper.R not found at: ", mapper_path, "\n",
                "Please ensure the script is in the correct location."))
  }
  
  if (!file.exists(annotator_path)) {
    stop(paste0("OT_annotator.R not found at: ", annotator_path, "\n",
                "Please ensure the script is in the correct location."))
  }
  
  log_info("Loading OT_mapper.R...")
  try_catch_with_message(
    source(mapper_path),
    "Failed to load OT_mapper.R"
  )
  
  log_info("Loading OT_annotator.R...")
  try_catch_with_message(
    source(annotator_path),
    "Failed to load OT_annotator.R"
  )
  
  log_success("Scripts loaded successfully")
  
  # Set annotation database paths
  if (is.null(opt$`exon-db`)) {
    exon_db <- file.path(script_dir, "data/UCSC_exons_modif_canonical.bed")
  } else {
    exon_db <- opt$`exon-db`
  }
  
  if (is.null(opt$`intron-db`)) {
    intron_db <- file.path(script_dir, "data/UCSC_introns_modif_canonical.bed")
  } else {
    intron_db <- opt$`intron-db`
  }
  
  log_info("Annotation databases:")
  cat("  Exon DB:", exon_db, "\n")
  cat("  Intron DB:", intron_db, "\n\n")
  
  if (!file.exists(exon_db)) {
    stop(paste0("Exon database file not found: ", exon_db, "\n",
                "Please ensure the annotation files are in the data/ directory or specify with -e option."))
  }
  
  if (!file.exists(intron_db)) {
    stop(paste0("Intron database file not found: ", intron_db, "\n",
                "Please ensure the annotation files are in the data/ directory or specify with -i option."))
  }
  
  # Step 1: GGGenome search
  log_step("Step 1: GGGenome Search")
  log_info("Searching for off-target candidates...")
  log_info("This may take a few moments depending on your internet connection...")
  
  try_catch_with_message(
    {
      list <- gggenome_to_dataframe(spacer, seed_length, pam)
    },
    "Failed to query GGGenome API. Please check your internet connection and try again."
  )
  
  if (is.null(list) || (nrow(list$plus_full) == 0 && nrow(list$minus_full) == 0)) {
    log_warning("No candidates found in GGGenome search.")
    log_warning("This may indicate a very specific gRNA with no off-targets.")
  } else {
    log_success("GGGenome search completed")
    cat("  Plus strand full sequence:", nrow(list$plus_full), "candidates\n")
    cat("  Minus strand full sequence:", nrow(list$minus_full), "candidates\n")
    cat("  Plus strand seed sequence:", nrow(list$plus_seed), "candidates\n")
    cat("  Minus strand seed sequence:", nrow(list$minus_seed), "candidates\n")
  }
  
  # Step 2: PAM filtering
  log_step("Step 2: PAM Exact Match Filtering")
  log_info("Filtering for exact PAM matches...")
  
  try_catch_with_message(
    {
      filtered_list <- filter_exact_pam(list, pam)
    },
    "Failed to filter PAM matches"
  )
  
  log_success("PAM filtering completed")
  
  # Step 3: Overlap detection
  log_step("Step 3: Overlap Detection")
  log_info("Finding overlaps between full and seed sequences...")
  log_info("This may take a moment for large datasets...")
  
  try_catch_with_message(
    {
      overlaps <- find_overlaps(filtered_list)
    },
    "Failed to detect overlaps. Please check if bedtools is installed and accessible."
  )
  
  plus_overlap_count <- nrow(overlaps$plus_overlaps)
  minus_overlap_count <- nrow(overlaps$minus_overlaps)
  total_overlap_count <- plus_overlap_count + minus_overlap_count
  
  log_success("Overlap detection completed")
  cat("  Plus strand overlaps:", plus_overlap_count, "\n")
  cat("  Minus strand overlaps:", minus_overlap_count, "\n")
  cat("  Total overlaps:", total_overlap_count, "\n")
  
  if (total_overlap_count == 0) {
    log_warning("No overlapping candidates found.")
    log_warning("This may indicate very stringent filtering criteria.")
  }
  
  # Step 4: Combine results
  log_step("Step 4: Combining Results")
  log_info("Combining plus and minus strand results...")
  
  try_catch_with_message(
    {
      combined_df <- combine_results(overlaps)
    },
    "Failed to combine results"
  )
  
  log_success("Results combined successfully")
  cat("  Total off-target candidates:", nrow(combined_df), "\n")
  
  if (nrow(combined_df) == 0) {
    log_warning("No off-target candidates found after filtering.")
    log_info("An empty output file will be created.")
  }
  
  # Step 5: Gene annotation
  log_step("Step 5: Gene Annotation")
  log_info("Mapping candidates to gene annotations...")
  log_info("This may take a moment depending on the number of candidates...")
  
  try_catch_with_message(
    {
      annotated_df <- annotate_with_bedtools(
        combined_df,
        exon_db = exon_db,
        intron_db = intron_db,
        output_file = output_file
      )
    },
    "Failed to annotate candidates with gene information"
  )
  
  # Step 6: Primer-BLAST link generation (if not skipped)
  primer_output_file <- NULL
  if (!opt$`skip-primer`) {
    log_step("Step 6: Primer-BLAST Link Generation")
    log_info("Generating Primer-BLAST URLs...")
    
    if (is.null(annotated_df) || nrow(annotated_df) == 0) {
      log_warning("No annotated candidates found. Skipping primer generation.")
      log_info("Primer-BLAST links require annotated candidates with chromosome coordinates.")
    } else {
      try_catch_with_message(
        {
          # Generate primer output filename
          primer_output_file <- gsub("\\.tsv$", "_with_primer.tsv", output_file)
          
          # Source primer generation function
          primer_script <- file.path(script_dir, "primer_generate.R")
          if (!file.exists(primer_script)) {
            log_warning("primer_generate.R not found. Skipping primer generation.")
          } else {
            # Load primer generation functions
            source(primer_script, local = TRUE)
            
            # Generate primer links
            generate_primer_blast_links(output_file, primer_output_file)
            
            if (file.exists(primer_output_file) && file.size(primer_output_file) > 0) {
              primer_count <- length(readLines(primer_output_file)) - 1  # Subtract header
              log_success(paste("Primer-BLAST links generated for", primer_count, "candidate(s)"))
            } else {
              log_warning("Primer-BLAST output file is empty")
            }
          }
        },
        "Failed during Primer-BLAST link generation"
      )
    }
  } else {
    log_info("Skipping Primer-BLAST link generation (--skip-primer specified)")
  }
  
  # Final summary
  log_step("Pipeline Complete")
  cat("\n")
  
  if (!is.null(annotated_df) && nrow(annotated_df) > 0) {
    log_success("Analysis completed successfully!")
    cat("  Off-target candidates detected:", nrow(combined_df), "\n")
    cat("  Annotated candidates:", nrow(annotated_df), "\n")
    cat("  Results saved to:", output_file, "\n")
    
    if (!is.null(primer_output_file) && file.exists(primer_output_file)) {
      cat("  Primer-BLAST links:", primer_output_file, "\n")
    }
    
    if (nrow(combined_df) > 0 && nrow(annotated_df) < nrow(combined_df)) {
      missing <- nrow(combined_df) - nrow(annotated_df)
      log_warning(paste0(missing, " candidate(s) were not annotated (no overlap with exons or introns)"))
      log_info("This is normal - not all genomic regions contain annotated genes")
    }
  } else {
    log_warning("No annotated candidates found.")
    log_info("This may indicate:")
    cat("  - No candidates overlap with annotated exons or introns\n")
    cat("  - Candidates are in intergenic regions\n")
    cat("  - Annotation database does not cover these regions\n")
    cat("\n")
    log_info("Output file created:", output_file)
    log_info("(File may be empty if no candidates were found)")
  }
  
  cat("\n")
  log_success("Done!")
}

# Run main function
main()
