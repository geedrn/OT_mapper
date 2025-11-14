#!/bin/bash

# =============================================================================
# OT Mapper - Map off-target candidates to gene annotations
# =============================================================================
# This script maps off-target candidates to gene annotations (exons/introns)
# using bedtools intersect with GENCODE annotation files.
#
# Usage:
#   bash OT_mapper.sh <input.bed> [options]
#
# Options:
#   -e, --exon-db      Exon annotation BED file [default: data/UCSC_exons_modif_canonical.bed]
#   -i, --intron-db    Intron annotation BED file [default: data/UCSC_introns_modif_canonical.bed]
#   -o, --output-dir   Output directory [default: analysis/OT_mapper_results]
#   -h, --help         Show this help message
# =============================================================================

set -euo pipefail

# =============================================================================
# Load Utilities
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/utils.sh"

# =============================================================================
# Configuration
# =============================================================================

readonly DEFAULT_EXON_DB="data/UCSC_exons_modif_canonical.bed"
readonly DEFAULT_INTRON_DB="data/UCSC_introns_modif_canonical.bed"
readonly DEFAULT_OUTPUT_DIR="analysis/OT_mapper_results"

# =============================================================================
# Global Variables
# =============================================================================

BED_INPUT=""
EXON_DB=${DEFAULT_EXON_DB}
INTRON_DB=${DEFAULT_INTRON_DB}
OUTPUT_DIR=${DEFAULT_OUTPUT_DIR}

# =============================================================================
# Helper Functions
# =============================================================================

show_usage() {
    cat << EOF
OT Mapper - Map off-target candidates to gene annotations

Usage:
  bash OT_mapper.sh <input.bed> [options]

Options:
  -e, --exon-db      Exon annotation BED file [default: ${DEFAULT_EXON_DB}]
  -i, --intron-db    Intron annotation BED file [default: ${DEFAULT_INTRON_DB}]
  -o, --output-dir   Output directory [default: ${DEFAULT_OUTPUT_DIR}]
  -h, --help         Show this help message

Examples:
  bash OT_mapper.sh analysis/OT/OT_candidate.bed
  bash OT_mapper.sh input.bed -e custom_exons.bed -i custom_introns.bed
EOF
}

parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -e|--exon-db)
                EXON_DB="$2"
                shift 2
                ;;
            -i|--intron-db)
                INTRON_DB="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                if [ -z "$BED_INPUT" ]; then
                    BED_INPUT="$1"
                else
                    error_exit "Multiple input files specified. Please provide only one input BED file."
                fi
                shift
                ;;
        esac
    done
}

validate_inputs() {
    # Check input file
    if [ -z "$BED_INPUT" ]; then
        error_exit "No input BED file specified.

Usage: bash OT_mapper.sh <input.bed> [options]
Use -h or --help for more information"
    fi
    
    check_file "$BED_INPUT" || error_exit "Input file validation failed"
    
    # Check if input is empty
    if [ ! -s "$BED_INPUT" ]; then
        log_warning "Input BED file is empty: $BED_INPUT"
        log_info "This may indicate no off-target candidates were found in the previous step."
    fi
    
    # Check annotation files
    check_file "$EXON_DB" || error_exit "Exon database file not found: $EXON_DB

Please check:
  - File path is correct (relative to current directory or absolute path)
  - File exists and is readable
  - Default location: ${DEFAULT_EXON_DB}

You can specify a custom path using: -e <path_to_exon_file>"
    
    check_file "$INTRON_DB" || error_exit "Intron database file not found: $INTRON_DB

Please check:
  - File path is correct (relative to current directory or absolute path)
  - File exists and is readable
  - Default location: ${DEFAULT_INTRON_DB}

You can specify a custom path using: -i <path_to_intron_file>"
    
    # Check if annotation files are empty
    [ ! -s "$EXON_DB" ] && log_warning "Exon database file is empty: $EXON_DB"
    [ ! -s "$INTRON_DB" ] && log_warning "Intron database file is empty: $INTRON_DB"
}

run_annotation() {
    local input_bed="$1"
    local exon_db="$2"
    local intron_db="$3"
    local output_file="$4"
    
    log_info "Mapping candidates to gene annotations using bedtools intersect..."
    log_info "This may take a moment depending on the number of candidates..."
    
    # Validate BED format (basic check)
    if ! head -n 1 "$input_bed" | grep -qE '^[^[:space:]]+[[:space:]]+[0-9]+[[:space:]]+[0-9]+'; then
        log_warning "Input BED file may not be in standard format. bedtools will attempt to process it anyway."
    fi
    
    # Run bedtools intersect
    # Output format: annotation name (col 5), chrom (col 1), start (col 2), end (col 3), gene info (cols 8-9)
    if ! bedtools intersect -a "$input_bed" \
                            -b "$exon_db" "$intron_db" \
                            -wb 2>/dev/null | \
         cut -f 5,1-3,8-9 > "$output_file" 2>/dev/null; then
        error_exit "bedtools intersect failed.

Possible reasons:
  - Invalid BED file format
  - bedtools version incompatibility
  - Corrupted annotation files
  - Insufficient memory

Please check:
  - Input BED file format (should be: chrom, start, end)
  - Annotation files are valid BED format
  - bedtools version: bedtools --version"
    fi
    
    [ -f "$output_file" ] || error_exit "Output file was not created: $output_file"
}

# =============================================================================
# Main Function
# =============================================================================

main() {
    log_step "OT Mapper - Gene Annotation Mapping"
    
    # Check prerequisites
    check_bedtools
    
    # Parse arguments
    parse_arguments "$@"
    
    # Validate inputs
    validate_inputs
    
    # Display configuration
    log_info "Configuration:"
    echo "  Input BED file: $BED_INPUT"
    echo "  Exon database: $EXON_DB"
    echo "  Intron database: $INTRON_DB"
    echo "  Output directory: $OUTPUT_DIR"
    echo ""
    
    # Count input candidates
    local input_count=$(count_lines "$BED_INPUT")
    if [ "$input_count" -gt 0 ]; then
        log_info "Input file contains $input_count candidate(s)"
    else
        log_warning "Input file contains no candidates. Output will be empty."
    fi
    
    # Create output directory
    ensure_directory "$OUTPUT_DIR"
    
    # Run annotation
    local output_file="${OUTPUT_DIR}/OT_mapped.tsv"
    run_annotation "$BED_INPUT" "$EXON_DB" "$INTRON_DB" "$output_file"
    
    # Count results
    local annotated_count=$(count_lines "$output_file")
    
    # Summary
    log_step "Mapping Complete"
    echo ""
    if [ "$annotated_count" -gt 0 ]; then
        log_success "Successfully annotated $annotated_count candidate(s)"
        log_success "Results saved to: $output_file"
        
        if [ "$input_count" -gt 0 ] && [ "$annotated_count" -lt "$input_count" ]; then
            local missing=$((input_count - annotated_count))
            log_warning "$missing candidate(s) were not annotated (no overlap with exons or introns)"
            log_info "This is normal - not all genomic regions contain annotated genes"
        fi
    else
        log_warning "No annotated candidates found."
        log_info "This may indicate:"
        echo "  - No candidates overlap with annotated exons or introns"
        echo "  - Candidates are in intergenic regions"
        echo "  - Annotation database does not cover these regions"
        echo ""
        log_info "An empty output file was created: $output_file"
    fi
    
    echo ""
    log_success "Done!"
}

# Run main function
main "$@"
