#!/bin/bash

# =============================================================================
# OT Mapper - Common Utilities Library
# =============================================================================
# This file contains shared functions and constants used across bash scripts
# Source this file in other scripts: source "$(dirname "$0")/utils.sh"
# =============================================================================

# =============================================================================
# Constants
# =============================================================================

readonly SCRIPT_VERSION="1.0"
readonly SCRIPT_AUTHOR="RN"

# Color codes for terminal output
readonly COLOR_RED='\033[0;31m'
readonly COLOR_GREEN='\033[0;32m'
readonly COLOR_YELLOW='\033[1;33m'
readonly COLOR_BLUE='\033[0;34m'
readonly COLOR_NC='\033[0m' # No Color

# =============================================================================
# Logging Functions
# =============================================================================

log_info() {
    echo -e "${COLOR_BLUE}[INFO]${COLOR_NC} $*"
}

log_success() {
    echo -e "${COLOR_GREEN}[SUCCESS]${COLOR_NC} $*"
}

log_warning() {
    echo -e "${COLOR_YELLOW}[WARNING]${COLOR_NC} $*"
}

log_error() {
    echo -e "${COLOR_RED}[ERROR]${COLOR_NC} $*" >&2
}

log_step() {
    echo ""
    echo -e "${COLOR_BLUE}=== $* ===${COLOR_NC}"
}

# =============================================================================
# Error Handling
# =============================================================================

error_exit() {
    log_error "$1"
    exit 1
}

# =============================================================================
# Validation Functions
# =============================================================================

# Validate spacer sequence
# Args: spacer sequence
# Returns: 0 if valid, 1 if invalid
validate_spacer() {
    local spacer="$1"
    
    if [ -z "$spacer" ]; then
        log_error "Spacer sequence is empty"
        return 1
    fi
    
    if [ ${#spacer} -ne 20 ]; then
        log_error "Invalid spacer length: ${#spacer} nucleotides (expected: 20)"
        return 1
    fi
    
    if ! [[ "$spacer" =~ ^[ACGTacgt]+$ ]]; then
        log_error "Invalid characters in spacer (must be A, C, G, T only)"
        return 1
    fi
    
    return 0
}

# Validate seed length
# Args: seed_length
# Returns: 0 if valid, 1 if invalid
validate_seed_length() {
    local seed_length="$1"
    
    if ! [[ "$seed_length" =~ ^[0-9]+$ ]]; then
        log_error "Seed length must be a number"
        return 1
    fi
    
    if [ "$seed_length" -lt 8 ] || [ "$seed_length" -gt 12 ]; then
        log_error "Invalid seed length: $seed_length (expected: 8-12)"
        return 1
    fi
    
    return 0
}

# =============================================================================
# Tool Checking Functions
# =============================================================================

# Check if a command exists
# Args: command_name
# Returns: 0 if exists, 1 if not
check_command() {
    command -v "$1" &> /dev/null
}

# Check required tools
# Args: tool1 tool2 ...
# Returns: 0 if all exist, exits with error if any missing
check_required_tools() {
    local missing_tools=()
    
    for tool in "$@"; do
        if ! check_command "$tool"; then
            missing_tools+=("$tool")
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        log_error "Missing required tools: ${missing_tools[*]}"
        return 1
    fi
    
    return 0
}

# Check bedtools availability
check_bedtools() {
    if ! check_command bedtools; then
        error_exit "bedtools is not installed or not in PATH.

Please install bedtools using one of the following methods:
  - conda: conda install -c bioconda bedtools
  - brew (macOS): brew install bedtools
  - apt (Ubuntu/Debian): sudo apt-get install bedtools
  - yum (CentOS/RHEL): sudo yum install bedtools

After installation, verify with: bedtools --version"
    fi
}

# Check download tool (wget or curl)
check_download_tool() {
    if check_command wget || check_command curl; then
        return 0
    else
        error_exit "Neither wget nor curl is available. Please install one of them to download data from GGGenome API."
    fi
}

# Check internet connectivity
check_internet() {
    log_info "Checking internet connectivity..."
    
    if check_command wget; then
        if wget -q --spider https://gggenome.dbcls.jp 2>/dev/null; then
            log_success "Internet connection OK"
            return 0
        fi
    elif check_command curl; then
        if curl -s --head https://gggenome.dbcls.jp > /dev/null 2>&1; then
            log_success "Internet connection OK"
            return 0
        fi
    fi
    
    log_warning "Cannot reach GGGenome API. Please check your internet connection."
    log_warning "The script will continue but may fail during data download."
    return 1
}

# =============================================================================
# File Operations
# =============================================================================

# Check if file exists and is readable
# Args: file_path
# Returns: 0 if valid, 1 if invalid
check_file() {
    local file="$1"
    
    if [ ! -f "$file" ]; then
        log_error "File not found: $file"
        return 1
    fi
    
    if [ ! -r "$file" ]; then
        log_error "File is not readable: $file"
        return 1
    fi
    
    return 0
}

# Create directory if it doesn't exist
# Args: directory_path
# Returns: 0 if successful, exits with error if fails
ensure_directory() {
    local dir="$1"
    
    if [ ! -d "$dir" ]; then
        if ! mkdir -p "$dir" 2>/dev/null; then
            error_exit "Failed to create directory: $dir"
        fi
        log_info "Created directory: $dir"
    fi
}

# =============================================================================
# Sequence Processing Functions
# =============================================================================

# Extract seed sequence from spacer+PAM
# Args: spacer pam seed_length
# Returns: seed sequence (echoed to stdout)
extract_seed_sequence() {
    local spacer="$1"
    local pam="$2"
    local seed_length="$3"
    local full_sequence="${spacer}${pam}"
    local total_length=${#full_sequence}
    
    # Calculate start position: (total_length - seed_length - pam_length) + 1
    # For 20nt spacer + 3nt PAM = 23nt total
    # Seed starts at: 23 - seed_length - 3 + 1 = 21 - seed_length
    local seed_start=$((total_length - seed_length - ${#pam} + 1))
    local seed_end=$total_length
    
    echo "$full_sequence" | cut -c ${seed_start}-${seed_end}
}

# Normalize sequence to uppercase
# Args: sequence
# Returns: uppercase sequence (echoed to stdout)
normalize_sequence() {
    echo "$1" | tr '[:lower:]' '[:upper:]'
}

# =============================================================================
# Download Functions
# =============================================================================

# Download from URL with retry
# Args: url output_file [max_retries]
# Returns: 0 if successful, 1 if failed
download_with_retry() {
    local url="$1"
    local output_file="$2"
    local max_retries="${3:-3}"
    local retry_count=0
    
    while [ $retry_count -lt $max_retries ]; do
        if check_command wget; then
            if wget -q --timeout=30 "$url" -O "$output_file" 2>/dev/null && [ -s "$output_file" ]; then
                return 0
            fi
        elif check_command curl; then
            if curl -s --max-time 30 "$url" -o "$output_file" 2>/dev/null && [ -s "$output_file" ]; then
                return 0
            fi
        fi
        
        retry_count=$((retry_count + 1))
        if [ $retry_count -lt $max_retries ]; then
            log_warning "Download failed, retrying ($retry_count/$max_retries)..."
            sleep 2
        fi
    done
    
    return 1
}

# =============================================================================
# BED File Operations
# =============================================================================

# Count lines in file (handles empty files)
# Args: file_path
# Returns: line count (echoed to stdout)
count_lines() {
    local file="$1"
    if [ -f "$file" ] && [ -s "$file" ]; then
        wc -l < "$file" | tr -d ' '
    else
        echo "0"
    fi
}

