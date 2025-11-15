#!/bin/bash
# =============================================================================
# OT Mapper Shiny App Launcher
# =============================================================================
# This script helps users launch the Shiny app easily, even if they're not
# familiar with R. It checks for dependencies and provides helpful error messages.
#
# Usage:
#   bash run_app.sh
#   or
#   chmod +x run_app.sh
#   ./run_app.sh
# =============================================================================

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

echo "=========================================="
echo "OT Mapper Shiny App Launcher"
echo "=========================================="
echo ""

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo -e "${RED}Error: R is not installed or not in PATH${NC}"
    echo ""
    echo "Please install R first:"
    echo "  - macOS: Download from https://cran.r-project.org/"
    echo "  - Linux: sudo apt-get install r-base (Ubuntu/Debian)"
    echo "  - Or use conda: conda install -c conda-forge r-base"
    exit 1
fi

echo -e "${GREEN}✓${NC} R is installed: $(Rscript --version | head -n 1)"
echo ""

# Check if bedtools is installed
if ! command -v bedtools &> /dev/null; then
    echo -e "${YELLOW}Warning: bedtools is not installed${NC}"
    echo ""
    echo "bedtools is required for gene annotation. The app will work but"
    echo "annotation features will be disabled."
    echo ""
    echo "To install bedtools:"
    echo "  - macOS: brew install bedtools"
    echo "  - Linux: sudo apt-get install bedtools"
    echo "  - Or use conda: conda install -c bioconda bedtools"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    echo -e "${GREEN}✓${NC} bedtools is installed: $(bedtools --version | head -n 1)"
fi
echo ""

# Check for required R packages
echo "Checking required R packages..."
Rscript -e "
required_packages <- c('shiny', 'shinydashboard', 'DT', 'httr', 'readr', 'dplyr', 'stringr', 'GenomicRanges')
missing_packages <- c()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat('Missing packages:', paste(missing_packages, collapse=', '), '\n')
  cat('Installing missing packages...\n')
  
  # Install from CRAN
  cran_packages <- missing_packages[missing_packages != 'GenomicRanges']
  if (length(cran_packages) > 0) {
    install.packages(cran_packages, repos='https://cran.rstudio.com/', quiet=TRUE)
  }
  
  # Install GenomicRanges from Bioconductor if needed
  if ('GenomicRanges' %in% missing_packages) {
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
      install.packages('BiocManager', repos='https://cran.rstudio.com/', quiet=TRUE)
    }
    BiocManager::install('GenomicRanges', quiet=TRUE)
  }
  
  cat('Packages installed successfully!\n')
} else {
  cat('All required packages are installed.\n')
}
" || {
    echo -e "${RED}Error: Failed to check/install R packages${NC}"
    exit 1
}

echo ""
echo -e "${GREEN}✓${NC} All dependencies are ready!"
echo ""

# Check if annotation data files exist
if [ ! -f "../R/scripts/data/UCSC_exons_modif_canonical.bed" ] || \
   [ ! -f "../R/scripts/data/UCSC_introns_modif_canonical.bed" ]; then
    echo -e "${YELLOW}Warning: Annotation database files not found${NC}"
    echo ""
    echo "The app will work, but gene annotation features will be disabled."
    echo "Annotation files should be in: R/scripts/data/"
    echo ""
    read -p "Continue anyway? (y/n) " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        exit 1
    fi
else
    echo -e "${GREEN}✓${NC} Annotation database files found"
fi

echo ""
echo "=========================================="
echo "Starting Shiny App..."
echo "=========================================="
echo ""
echo "The app will open in your default web browser."
echo "If it doesn't open automatically, go to: http://localhost:3838"
echo ""
echo "Press Ctrl+C to stop the app"
echo ""

# Run the Shiny app
Rscript -e "shiny::runApp(port=3838, launch.browser=TRUE)"

