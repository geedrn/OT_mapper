#!/bin/bash

# OT Mapper - Main workflow script
# Runs the complete pipeline: OT detection -> Annotation -> Primer generation
#
# Usage:
#   Interactive mode: bash main.sh
#   Command-line mode: bash main.sh -s <spacer> -l <seed_length> [options]

# Color codes for better log messages
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

log_step() {
    echo ""
    echo -e "${BLUE}=== $1 ===${NC}"
}

# Error handling function
error_exit() {
    log_error "$1"
    exit 1
}

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || error_exit "Failed to change to script directory: $SCRIPT_DIR"

log_step "OT Mapper - Complete Pipeline"
echo "This script runs the complete analysis pipeline:"
echo "  1. Off-target candidate detection"
echo "  2. Gene annotation mapping"
echo "  3. Primer-BLAST link generation"
echo ""

# Check if command-line arguments are provided
if [ $# -eq 0 ]; then
    # Interactive mode
    log_info "Starting in interactive mode..."
    log_info "You will be prompted for input parameters"
    echo ""
    bash scripts/OT_detector.sh
    DETECTOR_EXIT_CODE=$?
else
    # Command-line mode - pass all arguments to OT_detector.sh
    log_info "Starting in command-line mode..."
    bash scripts/OT_detector.sh "$@"
    DETECTOR_EXIT_CODE=$?
fi

# Check if OT detection was successful
if [ $DETECTOR_EXIT_CODE -ne 0 ]; then
    error_exit "OT detection failed with exit code $DETECTOR_EXIT_CODE.

Cannot proceed with annotation. Please fix the errors above and try again."
fi

if [ ! -f "analysis/OT/OT_candidate.bed" ]; then
    error_exit "OT detection completed but output file not found: analysis/OT/OT_candidate.bed

This may indicate:
  - OT detection script failed silently
  - Output directory was changed
  - File permissions issue

Please check the OT detection step output above."
fi

# Check if BED file is empty
if [ ! -s "analysis/OT/OT_candidate.bed" ]; then
    log_warning "OT candidate BED file is empty."
    log_warning "No candidates were detected. Pipeline will continue but annotation step may produce empty results."
fi

# Step 2: Map to gene annotations
log_step "Step 2: Gene Annotation"
log_info "Mapping off-target candidates to gene annotations..."

if ! bash scripts/OT_mapper.sh analysis/OT/OT_candidate.bed; then
    error_exit "Gene annotation failed.

Please check:
  - Annotation database files exist in data/ directory
  - bedtools is installed and working
  - Input BED file is valid

Cannot proceed with primer generation."
fi

# Check if annotation was successful
if [ ! -f "analysis/OT_mapper_results/OT_mapped.tsv" ]; then
    error_exit "Gene annotation completed but output file not found: analysis/OT_mapper_results/OT_mapped.tsv

Please check the annotation step output above."
fi

# Check if annotation file is empty
if [ ! -s "analysis/OT_mapper_results/OT_mapped.tsv" ]; then
    log_warning "Annotated TSV file is empty."
    log_warning "No candidates were annotated. This may be normal if candidates are in intergenic regions."
    log_info "Pipeline will continue but primer generation may produce empty results."
fi

# Step 3: Generate Primer-BLAST links
log_step "Step 3: Primer-BLAST Link Generation"
log_info "Generating Primer-BLAST URLs for each candidate..."

# Check if R is available
if ! command -v Rscript &> /dev/null; then
    error_exit "Rscript is not installed or not in PATH.

Please install R to generate Primer-BLAST links:
  - macOS: brew install r
  - Ubuntu/Debian: sudo apt-get install r-base
  - Or download from: https://www.r-project.org/

After installation, verify with: Rscript --version"
fi

# Check if primer_generate.R exists
if [ ! -f "scripts/primer_generate.R" ]; then
    error_exit "primer_generate.R not found: scripts/primer_generate.R

Please ensure the script exists in the scripts/ directory."
fi

# Check if optparse package is available
if ! Rscript -e "library(optparse)" 2>/dev/null; then
    log_warning "R package 'optparse' is not installed."
    log_info "Attempting to install optparse..."
    if Rscript -e "install.packages('optparse', repos='https://cran.rstudio.com/')" 2>/dev/null; then
        log_success "optparse installed successfully"
    else
        error_exit "Failed to install optparse R package.

Please install manually:
  Rscript -e \"install.packages('optparse', repos='https://cran.rstudio.com/')\""
    fi
fi

# Run primer generation
if ! Rscript scripts/primer_generate.R \
  -i analysis/OT_mapper_results/OT_mapped.tsv \
  -o analysis/OT_mapper_results/OT_with_primer.tsv 2>&1; then
    error_exit "Primer-BLAST link generation failed.

Please check:
  - Input TSV file format is correct
  - R and required packages are installed
  - File permissions are correct"
fi

# Verify output file was created
if [ ! -f "analysis/OT_mapper_results/OT_with_primer.tsv" ]; then
    log_warning "Primer-BLAST output file was not created."
    log_warning "This may indicate an error in the primer generation step."
else
    PRIMER_COUNT=$(wc -l < "analysis/OT_mapper_results/OT_with_primer.tsv" 2>/dev/null | tr -d ' ')
    if [ "$PRIMER_COUNT" -gt 0 ]; then
        log_success "Primer-BLAST links generated for $PRIMER_COUNT candidate(s)"
    else
        log_warning "Primer-BLAST output file is empty."
    fi
fi

# Final summary
log_step "Pipeline Complete"
echo ""
log_success "All steps completed successfully!"
echo ""
log_info "Output files:"
echo "  analysis/OT/OT_candidate.bed"
echo "  analysis/OT_mapper_results/OT_mapped.tsv"
echo "  analysis/OT_mapper_results/OT_with_primer.tsv"
echo ""

# Count results
if [ -f "analysis/OT/OT_candidate.bed" ] && [ -s "analysis/OT/OT_candidate.bed" ]; then
    CANDIDATE_COUNT=$(wc -l < "analysis/OT/OT_candidate.bed" | tr -d ' ')
    echo "  Total candidates detected: $CANDIDATE_COUNT"
fi

if [ -f "analysis/OT_mapper_results/OT_mapped.tsv" ] && [ -s "analysis/OT_mapper_results/OT_mapped.tsv" ]; then
    ANNOTATED_COUNT=$(wc -l < "analysis/OT_mapper_results/OT_mapped.tsv" | tr -d ' ')
    echo "  Annotated candidates: $ANNOTATED_COUNT"
fi

if [ -f "analysis/OT_mapper_results/OT_with_primer.tsv" ] && [ -s "analysis/OT_mapper_results/OT_with_primer.tsv" ]; then
    PRIMER_COUNT=$(wc -l < "analysis/OT_mapper_results/OT_with_primer.tsv" | tr -d ' ')
    echo "  Candidates with primer links: $PRIMER_COUNT"
fi

echo ""
log_success "Done!"
