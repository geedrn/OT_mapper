#!/bin/bash
# =============================================================================
# OT Mapper Shiny App - Double-Click Launcher (macOS)
# =============================================================================
# This file can be double-clicked in Finder to launch the Shiny app.
# Make sure this file has execute permissions (chmod +x)
# =============================================================================

# Get the directory where this script is located (project root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR/shiny"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Clear screen
clear

echo "=========================================="
echo "  OT Mapper Shiny App"
echo "=========================================="
echo ""
echo -e "${BLUE}Starting application...${NC}"
echo ""

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo -e "${RED}✗ Error: R is not installed${NC}"
    echo ""
    
    # Check if Homebrew is available
    if command -v brew &> /dev/null; then
        echo "Homebrew is available. Attempting to install R via Homebrew..."
        echo ""
        if brew install --cask r; then
            echo -e "${GREEN}✓${NC} R installed successfully via Homebrew"
            echo ""
            echo "Please restart this launcher after R installation completes."
            echo ""
            echo "Press any key to exit..."
            read -n 1
            exit 0
        else
            echo -e "${YELLOW}⚠ Failed to install R via Homebrew${NC}"
            echo ""
        fi
    fi
    
    echo "Please install R manually:"
    echo "  1. Visit: https://cran.r-project.org/"
    echo "  2. Download and install R for macOS"
    echo "  3. Make sure to add R to your PATH"
    echo ""
    echo "Or install Homebrew first, then R will be installed automatically:"
    echo "  /bin/bash -c \"\$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)\""
    echo ""
    echo "Press any key to exit..."
    read -n 1
    exit 1
fi

echo -e "${GREEN}✓${NC} R is installed"
echo ""

# Check if Homebrew is installed (for bedtools installation later)
if ! command -v brew &> /dev/null; then
    echo -e "${YELLOW}⚠ Homebrew not found${NC}"
    echo "  Homebrew is useful for installing bedtools automatically."
    echo ""
    read -p "  Would you like to install Homebrew now? (y/n) " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        echo ""
        echo "Installing Homebrew..."
        echo "This may take a few minutes..."
        echo ""
        /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
        
        # Add Homebrew to PATH (for Apple Silicon Macs)
        if [[ -f "/opt/homebrew/bin/brew" ]]; then
            eval "$(/opt/homebrew/bin/brew shellenv)"
        elif [[ -f "/usr/local/bin/brew" ]]; then
            eval "$(/usr/local/bin/brew shellenv)"
        fi
        
        if command -v brew &> /dev/null; then
            echo -e "${GREEN}✓${NC} Homebrew installed successfully"
        else
            echo -e "${YELLOW}⚠ Homebrew installation may require a terminal restart${NC}"
            echo "  Please restart this launcher after installation completes."
            echo ""
            echo "Press any key to exit..."
            read -n 1
            exit 0
        fi
    fi
    echo ""
else
    echo -e "${GREEN}✓${NC} Homebrew is installed"
fi
echo ""

# Check if bedtools is installed (optional)
if ! command -v bedtools &> /dev/null; then
    echo -e "${YELLOW}⚠ Warning: bedtools not found${NC}"
    echo "  Gene annotation features will be disabled without bedtools."
    echo ""
    
    # Try to install bedtools automatically
    if command -v brew &> /dev/null; then
        echo "  Attempting to install bedtools via Homebrew..."
        if brew install bedtools; then
            echo -e "${GREEN}✓${NC} bedtools installed successfully"
        else
            echo -e "${YELLOW}⚠ Failed to install bedtools automatically${NC}"
            echo "  You can install it manually later: brew install bedtools"
        fi
    elif command -v conda &> /dev/null; then
        echo "  Attempting to install bedtools via Conda..."
        if conda install -c bioconda bedtools -y; then
            echo -e "${GREEN}✓${NC} bedtools installed successfully"
        else
            echo -e "${YELLOW}⚠ Failed to install bedtools automatically${NC}"
            echo "  You can install it manually later: conda install -c bioconda bedtools"
        fi
    else
        echo "  To install bedtools:"
        echo "    - macOS: brew install bedtools"
        echo "    - Linux: sudo apt-get install bedtools (Ubuntu/Debian)"
        echo "    - Or use conda: conda install -c bioconda bedtools"
        echo ""
        read -p "  Continue without bedtools? (y/n) " -n 1 -r
        echo ""
        if [[ ! $REPLY =~ ^[Yy]$ ]]; then
            exit 1
        fi
    fi
    echo ""
else
    echo -e "${GREEN}✓${NC} bedtools is installed"
fi
echo ""

# Check and install R packages
echo "Checking R packages..."
Rscript -e "
required_packages <- c('shiny', 'shinydashboard', 'DT', 'httr', 'readr', 'dplyr', 'stringr', 'GenomicRanges')
missing_packages <- c()

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat('Installing missing packages (this may take a few minutes)...\n')
  
  cran_packages <- missing_packages[missing_packages != 'GenomicRanges']
  if (length(cran_packages) > 0) {
    install.packages(cran_packages, repos='https://cran.rstudio.com/', quiet=TRUE)
  }
  
  if ('GenomicRanges' %in% missing_packages) {
    if (!requireNamespace('BiocManager', quietly = TRUE)) {
      install.packages('BiocManager', repos='https://cran.rstudio.com/', quiet=TRUE)
    }
    BiocManager::install('GenomicRanges', quiet=TRUE)
  }
  
  cat('✓ Packages installed successfully!\n')
} else {
  cat('✓ All packages are ready\n')
}
" || {
    echo -e "${RED}✗ Error: Failed to install R packages${NC}"
    echo "Press any key to exit..."
    read -n 1
    exit 1
}

echo ""
echo "=========================================="
echo -e "${GREEN}Launching Shiny App...${NC}"
echo "=========================================="
echo ""

# Function to find an available port
find_available_port() {
    local start_port=$1
    local port=$start_port
    local max_port=$((start_port + 10))
    
    while [ $port -le $max_port ]; do
        # Check if port is available (macOS/Linux)
        if command -v lsof &> /dev/null; then
            if ! lsof -Pi :$port -sTCP:LISTEN -t >/dev/null 2>&1; then
                echo $port
                return 0
            fi
        elif command -v netstat &> /dev/null; then
            if ! netstat -an | grep -q ":$port.*LISTEN"; then
                echo $port
                return 0
            fi
        else
            # If we can't check, just try the port
            echo $port
            return 0
        fi
        port=$((port + 1))
    done
    
    # Fallback: return original port
    echo $start_port
    return 1
}

# Find an available port starting from 3838
APP_PORT=$(find_available_port 3838)

if [ "$APP_PORT" != "3838" ]; then
    echo -e "${YELLOW}⚠ Port 3838 is already in use${NC}"
    echo "   Using port $APP_PORT instead"
    echo ""
fi

echo "The app will open in your web browser automatically."
echo "If it doesn't open, go to: http://localhost:$APP_PORT"
echo ""
echo -e "${YELLOW}To stop the app, close this window or press Ctrl+C${NC}"
echo ""

# Run the Shiny app
Rscript -e "shiny::runApp(port=$APP_PORT, launch.browser=TRUE)" || {
    echo ""
    echo -e "${RED}✗ Error: Failed to start the app${NC}"
    echo ""
    echo "Possible solutions:"
    echo "  1. Close any other Shiny apps that might be running"
    echo "  2. Restart your computer"
    echo "  3. Try a different port manually:"
    echo "     cd $(pwd)"
    echo "     Rscript -e \"shiny::runApp(port=3839, launch.browser=TRUE)\""
    echo ""
    echo "Press any key to exit..."
    read -n 1
    exit 1
}

