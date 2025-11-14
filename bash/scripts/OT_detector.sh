#!/bin/bash

# =============================================================================
# OT Detector - Coordinate-based overlap detection using bedtools
# =============================================================================
# This script detects off-target candidates for CRISPR-Cas9 gRNAs by:
#   1. Querying GGGenome API for full and seed sequence matches
#   2. Filtering for exact PAM matches
#   3. Finding coordinate-based overlaps using bedtools
#   4. Outputting BED files with candidate coordinates
#
# Usage:
#   Interactive mode: bash OT_detector.sh
#   Command-line mode: bash OT_detector.sh -s <spacer> -l <seed_length> [options]
#
# Options:
#   -s, --spacer      20nt gRNA spacer sequence (without PAM) [required]
#   -l, --seed-length Seed sequence length: 8 to 12 [default: 12]
#   -p, --pam         PAM sequence [default: NGG]
#   -g, --genome      Genome assembly [default: hg38]
#   -f, --full-mismatch Mismatch tolerance for full sequence [default: 3]
#   -m, --seed-mismatch Mismatch tolerance for seed sequence [default: 1]
#   -o, --output-dir  Output directory [default: analysis]
#   --skip-annotation Skip gene annotation step [default: false]
#   --skip-primer     Skip Primer-BLAST link generation [default: false]
#   -h, --help        Show this help message
#
# Note: By default, this script runs the complete pipeline:
#   1. Off-target candidate detection
#   2. Gene annotation mapping (automatic)
#   3. Primer-BLAST link generation (automatic)
# =============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# =============================================================================
# Load Utilities
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/utils.sh"

# =============================================================================
# Configuration
# =============================================================================

readonly DEFAULT_SEED_LENGTH=12
readonly DEFAULT_PAM="NGG"
readonly DEFAULT_GENOME="hg38"
readonly DEFAULT_FULL_MISMATCH=3
readonly DEFAULT_SEED_MISMATCH=1
readonly DEFAULT_OUTPUT_DIR="analysis"

# =============================================================================
# Global Variables
# =============================================================================

SPACER=""
SEED_LENGTH=${DEFAULT_SEED_LENGTH}
PAM=${DEFAULT_PAM}
GENOME=${DEFAULT_GENOME}
FULL_MISMATCH=${DEFAULT_FULL_MISMATCH}
SEED_MISMATCH=${DEFAULT_SEED_MISMATCH}
OUTPUT_DIR=${DEFAULT_OUTPUT_DIR}
INTERACTIVE=false
SKIP_ANNOTATION=false
SKIP_PRIMER=false

# =============================================================================
# Helper Functions
# =============================================================================

# Display usage information
show_usage() {
    cat << EOF
OT Detector - Off-target candidate detection for CRISPR-Cas9

Usage:
  Interactive mode: bash OT_detector.sh
  Command-line mode: bash OT_detector.sh -s <spacer> [options]

Options:
  -s, --spacer        20nt gRNA spacer sequence (without PAM) [required]
  -l, --seed-length   Seed sequence length: 8 to 12 [default: ${DEFAULT_SEED_LENGTH}]
  -p, --pam           PAM sequence [default: ${DEFAULT_PAM}]
  -g, --genome        Genome assembly [default: ${DEFAULT_GENOME}]
  -f, --full-mismatch Mismatch tolerance for full sequence [default: ${DEFAULT_FULL_MISMATCH}]
  -m, --seed-mismatch Mismatch tolerance for seed sequence [default: ${DEFAULT_SEED_MISMATCH}]
  -o, --output-dir    Output directory [default: ${DEFAULT_OUTPUT_DIR}]
  --skip-annotation   Skip gene annotation step [default: false]
  --skip-primer       Skip Primer-BLAST link generation [default: false]
  -h, --help          Show this help message

Note: By default, this script runs the complete pipeline:
  1. Off-target candidate detection
  2. Gene annotation mapping (automatic)
  3. Primer-BLAST link generation (automatic)

Examples:
  bash OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12
  bash OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 8 -p NGG -g hg38
  bash OT_detector.sh -s GCTGAAGCACTGCACGCCGT -l 12 --skip-primer  # Skip primer generation
EOF
}

# Parse command-line arguments
parse_arguments() {
    while [[ $# -gt 0 ]]; do
        case $1 in
            -s|--spacer)
                SPACER="$2"
                shift 2
                ;;
            -l|--seed-length)
                SEED_LENGTH="$2"
                shift 2
                ;;
            -p|--pam)
                PAM="$2"
                shift 2
                ;;
            -g|--genome)
                GENOME="$2"
                shift 2
                ;;
            -f|--full-mismatch)
                FULL_MISMATCH="$2"
                shift 2
                ;;
            -m|--seed-mismatch)
                SEED_MISMATCH="$2"
                shift 2
                ;;
            -o|--output-dir)
                OUTPUT_DIR="$2"
                shift 2
                ;;
            --skip-annotation)
                SKIP_ANNOTATION=true
                shift
                ;;
            --skip-primer)
                SKIP_PRIMER=true
                shift
                ;;
            -h|--help)
                show_usage
                exit 0
                ;;
            *)
                log_error "Unknown option: $1"
                echo "Use -h or --help for usage information"
                exit 1
                ;;
        esac
done
}

# Validate all inputs
validate_inputs() {
    # Get spacer interactively if not provided
    if [ -z "$SPACER" ]; then
        INTERACTIVE=true
        echo ""
        read -p 'Input gRNA sequence (DO NOT INCLUDE PAM): ' SPACER
    fi
    
    # Validate spacer
    if ! validate_spacer "$SPACER"; then
        error_exit "Spacer validation failed. Please provide a 20-nucleotide sequence (A, C, G, T only)."
    fi
    
    SPACER=$(normalize_sequence "$SPACER")
    
    # Get seed length interactively if not provided
    if [ "$INTERACTIVE" = true ]; then
        echo ""
        read -p "Choose the length of seed sequence (8 to 12 nt) [default: ${DEFAULT_SEED_LENGTH}]: " input_seed
        if [ -n "$input_seed" ]; then
            SEED_LENGTH="$input_seed"
        fi
    fi
    
    # Validate seed length
    if ! validate_seed_length "$SEED_LENGTH"; then
        error_exit "Seed length validation failed. Please provide a number between 8 and 12."
    fi
}

# Prepare sequences (spacer+PAM and seed+PAM)
prepare_sequences() {
    local full_sequence="${SPACER}${PAM}"
    local seed_sequence
    
    seed_sequence=$(extract_seed_sequence "$SPACER" "$PAM" "$SEED_LENGTH")
    
    if [ -z "$seed_sequence" ] || [ ${#seed_sequence} -ne $((SEED_LENGTH + ${#PAM})) ]; then
        error_exit "Failed to extract seed sequence. This is an internal error."
    fi
    
    # Export for use in other functions
    export FULL_SEQUENCE="$full_sequence"
    export SEED_SEQUENCE="$seed_sequence"
    
    log_info "Sequences prepared:"
    echo "  Full sequence (spacer+PAM): $full_sequence"
    echo "  Seed sequence (seed+PAM): $seed_sequence"
    echo ""
}

# Parse GGGenome CSV file to BED format
# Args: csv_file bed_file strand
parse_gggenome_csv() {
    local csv_file="$1"
    local bed_file="$2"
    local strand="$3"
    
    if [ ! -f "$csv_file" ]; then
        log_error "CSV file not found: $csv_file"
        return 1
    fi
    
    if [ ! -s "$csv_file" ]; then
        log_warning "CSV file is empty: $csv_file"
        touch "$bed_file"
        return 0
    fi
    
    # Skip header (first 5 lines) and parse to BED format
    tail -n +6 "$csv_file" 2>/dev/null | \
    awk -F',' '
    {
        # Remove quotes from all fields
        for (i=1; i<=NF; i++) {
            gsub(/^"/, "", $i)
            gsub(/"$/, "", $i)
        }
        
        chrom = $1
        strand_val = $2
        start = $3
        end = $4
        sbjct = $9
        
        # Create BED format: chrom, start, end, name (sbjct), score (0), strand
        if (chrom != "" && start != "" && end != "" && sbjct != "") {
            print chrom "\t" start "\t" end "\t" sbjct "\t0\t" strand_val
        }
    }' > "$bed_file"
    
    [ $? -eq 0 ] || return 1
    return 0
}

# Filter BED file for PAM matches
# Args: input_bed output_bed strand pam_seq
filter_pam_bed() {
    local input_bed="$1"
    local output_bed="$2"
    local strand="$3"
    local pam_seq="$4"
    local pam_len=${#pam_seq}
    
    [ -f "$input_bed" ] || { log_error "Input BED file not found: $input_bed"; return 1; }
    
    if [ "$strand" = "+" ]; then
        # Plus strand: PAM at end
        if [ "$pam_seq" = "NGG" ]; then
            awk -F'\t' "length(\$4) >= $pam_len && \$4 ~ /.GG\$/ {print}" "$input_bed" > "$output_bed" 2>/dev/null
        else
            local pam_pattern=$(echo "$pam_seq" | sed 's/N/./g')
            awk -F'\t' "length(\$4) >= $pam_len && \$4 ~ /${pam_pattern}\$/ {print}" "$input_bed" > "$output_bed" 2>/dev/null
        fi
    else
        # Minus strand: PAM reverse complement at start
        if [ "$pam_seq" = "NGG" ]; then
            awk -F'\t' "length(\$4) >= $pam_len && \$4 ~ /^CC./ {print}" "$input_bed" > "$output_bed" 2>/dev/null
        else
            local pam_rc=$(echo "$pam_seq" | tr 'ACGTN' 'TGCAN' | rev | sed 's/N/./g')
            awk -F'\t' "length(\$4) >= $pam_len && \$4 ~ /^${pam_rc}/ {print}" "$input_bed" > "$output_bed" 2>/dev/null
        fi
    fi
    
    [ $? -eq 0 ] || return 1
    return 0
}

# Process a single strand (plus or minus)
# Args: strand temp_dir
process_strand() {
    local strand="$1"
    local temp_dir="$2"
    local strand_label
    
    [ "$strand" = "+" ] && strand_label="Plus" || strand_label="Minus"
    
    log_step "Step $([ "$strand" = "+" ] && echo "1" || echo "2"): Downloading ${strand_label} Strand Data"
    
    # Download full sequence
    local full_url="https://gggenome.dbcls.jp/ja/${GENOME}/${FULL_MISMATCH}/${strand}/nogap/${FULL_SEQUENCE}.csv"
    local full_csv="${temp_dir}/${strand}_full.csv"
    
    log_info "Downloading full sequence data (mismatch tolerance: $FULL_MISMATCH)..."
    if ! download_with_retry "$full_url" "$full_csv"; then
        error_exit "Failed to download ${strand_label,,} strand full sequence data from GGGenome.

URL: $full_url

Please check your internet connection and try again."
    fi
    log_success "Downloaded ${strand_label,,} strand full sequence data"
    
    # Download seed sequence
    local seed_url="https://gggenome.dbcls.jp/ja/${GENOME}/${SEED_MISMATCH}/${strand}/nogap/${SEED_SEQUENCE}.csv"
    local seed_csv="${temp_dir}/${strand}_seed.csv"
    
    log_info "Downloading seed sequence data (mismatch tolerance: $SEED_MISMATCH)..."
    if ! download_with_retry "$seed_url" "$seed_csv"; then
        error_exit "Failed to download ${strand_label,,} strand seed sequence data from GGGenome.

URL: $seed_url

Please check your internet connection and try again."
    fi
    log_success "Downloaded ${strand_label,,} strand seed sequence data"
    
    # Parse CSV files to BED
    log_info "Processing ${strand_label,,} strand data..."
    parse_gggenome_csv "$full_csv" "${temp_dir}/${strand}_full.bed" "$strand" || return 1
    parse_gggenome_csv "$seed_csv" "${temp_dir}/${strand}_seed.bed" "$strand" || return 1
    
    local full_count=$(count_lines "${temp_dir}/${strand}_full.bed")
    local seed_count=$(count_lines "${temp_dir}/${strand}_seed.bed")
    log_info "${strand_label} strand candidates: $full_count full sequence, $seed_count seed sequence"
    
    # Filter PAM
    log_info "Filtering for exact PAM match (PAM: $PAM)..."
    filter_pam_bed "${temp_dir}/${strand}_full.bed" "${temp_dir}/${strand}_full_pam.bed" "$strand" "$PAM" || return 1
    filter_pam_bed "${temp_dir}/${strand}_seed.bed" "${temp_dir}/${strand}_seed_pam.bed" "$strand" "$PAM" || return 1
    
    local full_pam_count=$(count_lines "${temp_dir}/${strand}_full_pam.bed")
    local seed_pam_count=$(count_lines "${temp_dir}/${strand}_seed_pam.bed")
    log_info "After PAM filtering: $full_pam_count full sequence, $seed_pam_count seed sequence"
    
    # Find overlaps
    log_info "Finding overlaps using bedtools intersect..."
    if ! bedtools intersect -b "${temp_dir}/${strand}_seed_pam.bed" \
                            -a "${temp_dir}/${strand}_full_pam.bed" \
                            -wa > "${temp_dir}/${strand}_overlaps.bed" 2>/dev/null; then
        error_exit "bedtools intersect failed for ${strand_label,,} strand data."
    fi
    
    local overlap_count=$(count_lines "${temp_dir}/${strand}_overlaps.bed")
    log_success "Found $overlap_count ${strand_label,,} strand overlaps"
    
    echo "$overlap_count"
}

# Extract CSV data for matched coordinates
# Args: csv_file bed_file output_file temp_dir
extract_csv_data() {
    local csv_file="$1"
    local bed_file="$2"
    local output_file="$3"
    local temp_dir="$4"
    
    [ -f "$csv_file" ] && [ -f "$bed_file" ] || {
        log_warning "Skipping CSV extraction: missing input files"
        touch "$output_file"
        return 0
    }
    
    # Create coordinate lookup
    awk -F'\t' '{print $1 ":" $2 "-" $3}' "$bed_file" > "${temp_dir}/coords.txt" 2>/dev/null
    
    # Match and extract
    tail -n +6 "$csv_file" 2>/dev/null | \
    awk -F',' -v coords_file="${temp_dir}/coords.txt" '
    BEGIN {
        while ((getline line < coords_file) > 0) coords[line] = 1
        close(coords_file)
    }
    {
        chrom = $1; start = $3; end = $4
        gsub(/^"|"$/, "", chrom)
        gsub(/^"|"$/, "", start)
        gsub(/^"|"$/, "", end)
        coord = chrom ":" start "-" end
        if (coord in coords) print $0
    }' > "$output_file" 2>/dev/null
    
    return 0
}

# Create output files
# Args: temp_dir output_dir
create_output_files() {
    local temp_dir="$1"
    local output_dir="$2"
    
    log_step "Step 3: Creating Output Files"
    
    ensure_directory "${output_dir}/OT"
    ensure_directory "${output_dir}/intermediate"
    
    # Create OT_candidate.bed
    awk -F'\t' '{print $1 "\t" $2 "\t" $3}' \
        "${temp_dir}/all_overlaps.bed" > "${output_dir}/OT/OT_candidate.bed" 2>/dev/null || \
        error_exit "Failed to create OT_candidate.bed"
    log_success "Created: ${output_dir}/OT/OT_candidate.bed"
    
    # Extract CSV data
    log_info "Extracting full CSV data for overlapping candidates..."
    extract_csv_data "${temp_dir}/plus_full.csv" "${temp_dir}/plus_overlaps.bed" \
                    "${temp_dir}/plus_matched.csv" "$temp_dir"
    extract_csv_data "${temp_dir}/minus_full.csv" "${temp_dir}/minus_overlaps.bed" \
                    "${temp_dir}/minus_matched.csv" "$temp_dir"
    
    # Combine CSV data
    cat "${temp_dir}/plus_matched.csv" "${temp_dir}/minus_matched.csv" > \
        "${output_dir}/OT/OT_list_final.csv" 2>/dev/null || \
        error_exit "Failed to create OT_list_final.csv"
    log_success "Created: ${output_dir}/OT/OT_list_final.csv"
    
    # Create UCSC_list_final.csv
    awk -F',' '{
        for (i=1; i<=NF; i++) { gsub(/^"|"$/, "", $i) }
        print $9 "," $1 ":" $3 "-" $4
    }' "${output_dir}/OT/OT_list_final.csv" > "${output_dir}/OT/UCSC_list_final.csv" 2>/dev/null || \
        error_exit "Failed to create UCSC_list_final.csv"
    log_success "Created: ${output_dir}/OT/UCSC_list_final.csv"
    
    # Move intermediate files
    log_info "Moving intermediate files..."
    mv "${temp_dir}"/*.csv "${output_dir}/intermediate/" 2>/dev/null || true
    mv "${temp_dir}"/*.bed "${output_dir}/intermediate/" 2>/dev/null || true
}

# =============================================================================
# Main Function
# =============================================================================

main() {
    log_step "OT Detector - Off-target Candidate Detection"
    echo "Version: ${SCRIPT_VERSION} (created by ${SCRIPT_AUTHOR})"
    echo ""
    
    # Check prerequisites
    log_info "Checking prerequisites..."
    check_bedtools
    check_download_tool
    check_internet || true  # Warning only, don't fail
    
    # Parse arguments
    parse_arguments "$@"
    
    # Validate inputs
    validate_inputs
    
    log_success "Input validation passed"
    echo ""
    log_info "Configuration:"
    echo "  Spacer: $SPACER"
    echo "  PAM: $PAM"
    echo "  Seed length: $SEED_LENGTH"
    echo "  Genome: $GENOME"
    echo "  Full sequence mismatch tolerance: $FULL_MISMATCH"
    echo "  Seed sequence mismatch tolerance: $SEED_MISMATCH"
    echo "  Output directory: $OUTPUT_DIR"
    echo ""
    
    # Prepare sequences
    prepare_sequences
    
    # Create temporary directory
    TEMP_DIR=$(mktemp -d) || error_exit "Failed to create temporary directory"
    trap "rm -rf $TEMP_DIR" EXIT
    log_info "Temporary directory created: $TEMP_DIR"
    
    # Process both strands
    PLUS_COUNT=$(process_strand "+" "$TEMP_DIR")
    MINUS_COUNT=$(process_strand "-" "$TEMP_DIR")
    
    # Combine results
    log_step "Step 3: Combining Results"
    cat "${TEMP_DIR}/plus_overlaps.bed" "${TEMP_DIR}/minus_overlaps.bed" > \
        "${TEMP_DIR}/all_overlaps.bed" 2>/dev/null || \
        error_exit "Failed to combine results"
    
    TOTAL_COUNT=$((PLUS_COUNT + MINUS_COUNT))
    log_success "Total off-target candidates detected: $TOTAL_COUNT"
    echo "  Plus strand: $PLUS_COUNT candidates"
    echo "  Minus strand: $MINUS_COUNT candidates"
    
    if [ "$TOTAL_COUNT" -eq 0 ]; then
        log_warning "No off-target candidates found."
        log_info "This may indicate a very specific gRNA with no off-targets."
    fi
    
    # Create output files
    create_output_files "$TEMP_DIR" "$OUTPUT_DIR"
    
    # Step 2: Gene annotation (if not skipped)
    if [ "$SKIP_ANNOTATION" = false ]; then
        log_step "Step 4: Gene Annotation"
        
        # Get script directory
        local script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
        local mapper_script="${script_dir}/OT_mapper.sh"
        local candidate_bed="${OUTPUT_DIR}/OT/OT_candidate.bed"
        
        if [ ! -f "$mapper_script" ]; then
            log_warning "OT_mapper.sh not found. Skipping annotation step."
            log_info "You can run annotation manually: bash $mapper_script $candidate_bed"
        elif [ ! -f "$candidate_bed" ] || [ ! -s "$candidate_bed" ]; then
            log_warning "OT candidate BED file is empty or missing. Skipping annotation step."
        else
            log_info "Running gene annotation mapping..."
            if bash "$mapper_script" "$candidate_bed" -o "${OUTPUT_DIR}/OT_mapper_results"; then
                log_success "Gene annotation completed"
                
                # Step 3: Primer-BLAST link generation (if not skipped)
                if [ "$SKIP_PRIMER" = false ]; then
                    log_step "Step 5: Primer-BLAST Link Generation"
                    
                    local mapped_tsv="${OUTPUT_DIR}/OT_mapper_results/OT_mapped.tsv"
                    local primer_tsv="${OUTPUT_DIR}/OT_mapper_results/OT_with_primer.tsv"
                    local primer_script="${script_dir}/primer_generate.R"
                    
                    if [ ! -f "$primer_script" ]; then
                        log_warning "primer_generate.R not found. Skipping primer generation."
                    elif [ ! -f "$mapped_tsv" ] || [ ! -s "$mapped_tsv" ]; then
                        log_warning "Annotated TSV file is empty or missing. Skipping primer generation."
                    else
                        # Check if R is available
                        if ! check_command Rscript; then
                            log_warning "Rscript not found. Skipping primer generation."
                            log_info "Install R to generate Primer-BLAST links"
                        else
                            log_info "Generating Primer-BLAST URLs..."
                            if Rscript "$primer_script" -i "$mapped_tsv" -o "$primer_tsv" 2>/dev/null; then
                                local primer_count=$(count_lines "$primer_tsv")
                                if [ "$primer_count" -gt 0 ]; then
                                    log_success "Primer-BLAST links generated for $primer_count candidate(s)"
                                else
                                    log_warning "Primer-BLAST output file is empty"
                                fi
                            else
                                log_warning "Primer generation failed. Check R and optparse package installation."
                            fi
                        fi
                    fi
                else
                    log_info "Skipping Primer-BLAST link generation (--skip-primer specified)"
                fi
            else
                log_warning "Gene annotation failed. Check error messages above."
            fi
        fi
    else
        log_info "Skipping gene annotation step (--skip-annotation specified)"
    fi
    
    # Final summary
    log_step "Pipeline Complete"
    echo ""
    log_success "Results saved to:"
    echo "  ${OUTPUT_DIR}/OT/OT_candidate.bed ($TOTAL_COUNT candidates)"
    echo "  ${OUTPUT_DIR}/OT/OT_list_final.csv"
    echo "  ${OUTPUT_DIR}/OT/UCSC_list_final.csv"
    
    if [ "$SKIP_ANNOTATION" = false ] && [ -f "${OUTPUT_DIR}/OT_mapper_results/OT_mapped.tsv" ]; then
        local annotated_count=$(count_lines "${OUTPUT_DIR}/OT_mapper_results/OT_mapped.tsv")
        echo "  ${OUTPUT_DIR}/OT_mapper_results/OT_mapped.tsv ($annotated_count annotated)"
    fi
    
    if [ "$SKIP_PRIMER" = false ] && [ "$SKIP_ANNOTATION" = false ] && [ -f "${OUTPUT_DIR}/OT_mapper_results/OT_with_primer.tsv" ]; then
        local primer_count=$(count_lines "${OUTPUT_DIR}/OT_mapper_results/OT_with_primer.tsv")
        echo "  ${OUTPUT_DIR}/OT_mapper_results/OT_with_primer.tsv ($primer_count with primer links)"
    fi
    
    echo ""
    log_success "Done!"
}

# Run main function
main "$@"
