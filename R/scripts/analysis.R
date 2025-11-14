#!/usr/bin/env Rscript

# CRISPR Off-target Analysis Script
#
# Usage:
# Rscript analysis.R [spacer] [seed_length] [PAM] [output_file]
#
# Example:
# Rscript analysis.R TCGCCCAGCGACCCTGCTCC 8 NGG results.tsv

# Suppress package loading messages
suppressPackageStartupMessages({
  options(warn = 1)  # Show warnings immediately
})

# Logging functions
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

# Error handling wrapper
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

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
spacer <- if (length(args) >= 1) args[1] else "TCGCCCAGCGACCCTGCTCC"
seed_length <- if (length(args) >= 2) as.numeric(args[2]) else 8
pam <- if (length(args) >= 3) args[3] else "NGG"
output_file <- if (length(args) >= 4) args[4] else "annotated_offtargets.tsv"

log_step("CRISPR Off-target Analysis")
cat("Version: 1.0\n\n")

# Validate inputs
log_info("Validating inputs...")

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

if (is.na(seed_length) || seed_length < 8 || seed_length > 12) {
  stop(paste0("Invalid seed length: ", seed_length, "\n",
              "Expected: Integer between 8 and 12 (inclusive)\n",
              "Received: ", seed_length))
}

if (nchar(pam) == 0) {
  stop("PAM sequence cannot be empty.")
}

log_success("Input validation passed")

log_info("Configuration:")
cat("  Spacer:", spacer, "\n")
cat("  Seed length:", seed_length, "\n")
cat("  PAM:", pam, "\n")
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
exon_db <- file.path(script_dir, "data/UCSC_exons_modif_canonical.bed")
intron_db <- file.path(script_dir, "data/UCSC_introns_modif_canonical.bed")

log_info("Annotation databases:")
cat("  Exon DB:", exon_db, "\n")
cat("  Intron DB:", intron_db, "\n\n")

if (!file.exists(exon_db)) {
  stop(paste0("Exon database file not found: ", exon_db, "\n",
              "Please ensure the annotation files are in the data/ directory."))
}

if (!file.exists(intron_db)) {
  stop(paste0("Intron database file not found: ", intron_db, "\n",
              "Please ensure the annotation files are in the data/ directory."))
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

# Final summary
log_step("Analysis Complete")
cat("\n")

if (!is.null(annotated_df) && nrow(annotated_df) > 0) {
  log_success("Analysis completed successfully!")
  cat("  Off-target candidates detected:", nrow(combined_df), "\n")
  cat("  Annotated candidates:", nrow(annotated_df), "\n")
  cat("  Results saved to:", output_file, "\n")
  
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
